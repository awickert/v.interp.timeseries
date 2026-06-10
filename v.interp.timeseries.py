#!/usr/bin/env python3
############################################################################
#
# MODULE:       v.interp.timeseries
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Interpolate station time series to raster grids or sample points/areas
#
# COPYRIGHT:    (c) 2026 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################

#%module
#% description: Interpolate station time series to raster grids or sample points/areas
#% keyword: vector
#% keyword: raster
#% keyword: temporal
#% keyword: interpolation
#% keyword: hydrology
#% keyword: climate
#% keyword: precipitation
#%end

#%option G_OPT_V_INPUT
#%  key: input
#%  label: Input vector map of station locations (from v.in.ghcn or similar)
#%  required: yes
#%end

#%option
#%  key: table
#%  type: string
#%  label: Time series table name (default: {input}_timeseries)
#%  required: no
#%end

#%option
#%  key: element
#%  type: string
#%  label: Climate element to interpolate (e.g. PRCP, TMAX, TMIN)
#%  answer: PRCP
#%  required: yes
#%end

#%option
#%  key: start_date
#%  type: string
#%  label: Start date (YYYY-MM-DD; omit for first available record)
#%  required: no
#%end

#%option
#%  key: end_date
#%  type: string
#%  label: End date (YYYY-MM-DD; omit for last available record)
#%  required: no
#%end

#%option
#%  key: method
#%  type: string
#%  label: Interpolation method
#%  options: idw,spline,regression
#%  answer: idw
#%  required: yes
#%end

#%option
#%  key: covariates
#%  type: string
#%  label: Covariate raster maps for regression (comma-separated, e.g. dem,slope)
#%  required: no
#%end

#%option
#%  key: power
#%  type: double
#%  label: IDW power parameter
#%  answer: 2.0
#%  required: no
#%end

#%option
#%  key: npoints
#%  type: integer
#%  label: Nearest stations to use for IDW (0 = use all stations)
#%  answer: 0
#%  required: no
#%end

#%option G_OPT_R_OUTPUT
#%  key: output
#%  label: Output raster basename (one map per time step)
#%  required: no
#%end

#%option G_OPT_V_INPUT
#%  key: sample
#%  label: Sample points or area polygons for extracted time series
#%  required: no
#%end

#%option
#%  key: min_stations
#%  type: integer
#%  label: Minimum stations with data required to interpolate a time step
#%  answer: 4
#%  required: no
#%end

#%flag
#%  key: t
#%  description: Register output rasters as a space-time raster dataset (strds)
#%end

#%flag
#%  key: f
#%  description: Fill stations with missing data from nearest neighbor with data
#%end

#%flag
#%  key: m
#%  description: Mask output to the convex hull of active stations
#%end

import os
import sys
import math
import sqlite3
import tempfile
import atexit
import datetime as dt_mod

import numpy as np

import grass.script as gs

_TMPFILES = []


def cleanup():
    for f in _TMPFILES:
        try:
            os.unlink(f)
        except OSError:
            pass


def _tmpfile(suffix=''):
    path = tempfile.mktemp(suffix=suffix)
    _TMPFILES.append(path)
    return path


def require_package(pkg):
    try:
        __import__(pkg)
    except ImportError:
        gs.fatal(
            "Python package '{}' is required. "
            "Install with: pip3 install --break-system-packages {}".format(pkg, pkg)
        )


def get_mapset_db():
    env = gs.gisenv()
    return os.path.join(
        env['GISDBASE'], env['LOCATION_NAME'], env['MAPSET'],
        'sqlite', 'sqlite.db'
    )


def get_point_coords(vector_map, type='point'):
    """Return dict cat -> (x, y) for points or centroids."""
    coords = {}
    proc = gs.pipe_command(
        'v.out.ascii', input=vector_map, format='point',
        type=type, separator='pipe'
    )
    for line in proc.stdout:
        line = line.decode().strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('|')
        if len(parts) < 3:
            continue
        try:
            x, y, cat = float(parts[0]), float(parts[1]), int(parts[2])
            coords[cat] = (x, y)
        except ValueError:
            continue
    proc.wait()
    return coords


def build_grid_xy(region):
    """Return (rows*cols, 2) array of (x, y) cell centres, north-to-south order."""
    rows, cols = region['rows'], region['cols']
    nsres, ewres = region['nsres'], region['ewres']
    y_ctrs = np.linspace(region['n'] - nsres / 2.0, region['s'] + nsres / 2.0, rows)
    x_ctrs = np.linspace(region['w'] + ewres / 2.0, region['e'] - ewres / 2.0, cols)
    xx, yy = np.meshgrid(x_ctrs, y_ctrs)   # shape (rows, cols); yy decreasing
    return np.column_stack([xx.ravel(), yy.ravel()])


def read_raster(map_name, region, null_val=-9999.0):
    """Read a GRASS float raster into a float32 numpy array (north→south)."""
    tmp = _tmpfile('.bin')
    gs.run_command('r.out.bin', flags='f', input=map_name,
                   output=tmp, bytes=4, null=null_val)
    arr = np.fromfile(tmp, dtype=np.float32).reshape(region['rows'], region['cols'])
    arr[arr == np.float32(null_val)] = np.nan
    return arr


def write_raster(arr, map_name, region, null_val=-9999.0):
    """Write a float32 numpy array (north→south) as a GRASS FCELL raster."""
    out = arr.astype(np.float32).copy()
    out[np.isnan(out)] = np.float32(null_val)
    tmp = _tmpfile('.bin')
    out.tofile(tmp)
    gs.run_command(
        'r.in.bin', flags='f', input=tmp, output=map_name,
        bytes=4, anull=null_val,
        north=region['n'], south=region['s'],
        east=region['e'], west=region['w'],
        rows=region['rows'], cols=region['cols'],
        overwrite=True
    )


def rasterize_to_cats(sample_map, region):
    """Rasterize polygon map to a numpy int32 array of cat values (0 = outside)."""
    tmp_cell = 'v_interp_ts_zones_{}'.format(os.getpid())
    tmp_fcell = 'v_interp_ts_zones_f_{}'.format(os.getpid())
    try:
        gs.run_command('v.to.rast', input=sample_map, output=tmp_cell,
                       use='cat', overwrite=True)
        gs.run_command('r.mapcalc',
                       expression='{f} = float({c})'.format(f=tmp_fcell, c=tmp_cell),
                       overwrite=True)
        tmp = _tmpfile('.bin')
        gs.run_command('r.out.bin', flags='f', input=tmp_fcell,
                       output=tmp, bytes=4, null=-1)
        raw = np.fromfile(tmp, dtype=np.float32).reshape(region['rows'], region['cols'])
        zone_arr = np.where(raw < 0, 0, np.round(raw).astype(np.int32))
    finally:
        for name in (tmp_cell, tmp_fcell):
            if gs.find_file(name, element='raster')['name']:
                gs.run_command('g.remove', type='raster', name=name, flags='f')
    return zone_arr


def sample_covariates_at_xy(cov_flat_list, xy, region):
    """Nearest-cell sample of covariate rasters at given (x, y) coordinates.
    Returns array of shape (n_points, n_covariates)."""
    rows, cols = region['rows'], region['cols']
    col_idx = np.clip(
        ((xy[:, 0] - region['w']) / region['ewres']).astype(int), 0, cols - 1)
    row_idx = np.clip(
        ((region['n'] - xy[:, 1]) / region['nsres']).astype(int), 0, rows - 1)
    results = []
    for cov_flat in cov_flat_list:
        cov_2d = cov_flat.reshape(rows, cols)
        results.append(cov_2d[row_idx, col_idx])
    return np.column_stack(results) if results else np.empty((len(xy), 0))


def interp_idw(stn_xy, vals, query_xy, power, npoints):
    from scipy.spatial import cKDTree
    k = min(npoints if npoints > 0 else len(stn_xy), len(stn_xy))
    tree = cKDTree(stn_xy)
    dists, idxs = tree.query(query_xy, k=k)
    if k == 1:
        return vals[idxs.ravel()]
    dists = np.where(dists == 0, 1e-12, dists)
    weights = 1.0 / dists**power
    weights /= weights.sum(axis=1, keepdims=True)
    return (weights * vals[idxs]).sum(axis=1)


def interp_spline(stn_xy, vals, query_xy):
    from scipy.interpolate import RBFInterpolator
    rbf = RBFInterpolator(stn_xy, vals, kernel='thin_plate_spline')
    return rbf(query_xy)


def interp_regression(stn_xy, vals, cov_st, query_xy, cov_q):
    """Linear regression on covariates + thin-plate-spline on residuals."""
    from scipy.interpolate import RBFInterpolator
    A_st = np.column_stack([np.ones(len(vals)), cov_st])
    coeffs, _, _, _ = np.linalg.lstsq(A_st, vals, rcond=None)
    trend_st = A_st @ coeffs
    resid = vals - trend_st
    rbf = RBFInterpolator(stn_xy, resid, kernel='thin_plate_spline')
    A_q = np.column_stack([np.ones(len(query_xy)), cov_q])
    return A_q @ coeffs + rbf(query_xy)


def interpolate(method, stn_xy, vals, query_xy, power, npoints,
                cov_st=None, cov_q=None):
    if method == 'idw':
        return interp_idw(stn_xy, vals, query_xy, power, npoints)
    elif method == 'spline':
        return interp_spline(stn_xy, vals, query_xy)
    elif method == 'regression':
        if cov_st is not None and cov_q is not None and cov_st.shape[1] > 0:
            return interp_regression(stn_xy, vals, cov_st, query_xy, cov_q)
        return interp_spline(stn_xy, vals, query_xy)
    raise ValueError("Unknown method: {}".format(method))


def loo_rmse(method, stn_xy, vals, power, npoints, cov_st=None, max_n=50):
    """Leave-one-out cross-validation RMSE. Returns NaN if skipped."""
    n = len(vals)
    if n < 2:
        return np.nan
    if method == 'spline' and n > max_n:
        return np.nan    # too expensive; caller can note this
    sq_errs = []
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        xy_tr = stn_xy[mask]
        v_tr = vals[mask]
        cov_tr = cov_st[mask] if cov_st is not None else None
        cov_i = cov_st[[i]] if cov_st is not None else None
        try:
            pred = interpolate(method, xy_tr, v_tr, stn_xy[[i]],
                               power, npoints, cov_tr, cov_i)[0]
        except Exception:
            continue
        sq_errs.append((float(pred) - vals[i]) ** 2)
    return math.sqrt(np.mean(sq_errs)) if sq_errs else np.nan


def neighbor_fill(all_cats, cat_to_xy, cat_to_val, active_cats):
    """Return (filled_val_dict, filled_cats_set) with gaps borrowed from nearest neighbor."""
    from scipy.spatial import cKDTree
    active = sorted(active_cats)
    active_xy = np.array([cat_to_xy[c] for c in active])
    active_vals = np.array([cat_to_val[c] for c in active])
    tree = cKDTree(active_xy)
    filled_val = dict(cat_to_val)
    filled_cats = set(active_cats)
    for cat in all_cats:
        if cat in filled_cats:
            continue
        _, idx = tree.query([cat_to_xy[cat]], k=1)
        filled_val[cat] = float(active_vals[idx[0]])
        filled_cats.add(cat)
    return filled_val, filled_cats


def register_strds(output_base, element, map_date_pairs):
    """Create a strds and register all output maps."""
    gs.run_command('t.create', type='strds', temporaltype='absolute',
                   output=output_base,
                   title='{} {}'.format(element, output_base),
                   description='Interpolated {} from station time series'.format(element),
                   overwrite=True)
    reg_file = _tmpfile('.txt')
    days_set = {d[8:] for _, d in map_date_pairs if len(d) >= 10}
    monthly = (days_set == {'01'})
    with open(reg_file, 'w') as f:
        for map_name, date_str in map_date_pairs:
            if monthly:
                yr, mo = int(date_str[:4]), int(date_str[5:7])
                if mo == 12:
                    end_str = '{:04d}-01-01'.format(yr + 1)
                else:
                    end_str = '{:04d}-{:02d}-01'.format(yr, mo + 1)
            else:
                d = dt_mod.date.fromisoformat(date_str)
                end_str = (d + dt_mod.timedelta(days=1)).isoformat()
            f.write('{}|{}|{}\n'.format(map_name, date_str, end_str))
    gs.run_command('t.register', input=output_base, file=reg_file, overwrite=True)
    gs.message("Registered {} maps in strds '{}'.".format(len(map_date_pairs), output_base))


def main():
    options, flags = gs.parser()
    atexit.register(cleanup)

    input_map   = options['input']
    table_name  = options['table'] or '{}_timeseries'.format(input_map.split('@')[0])
    element     = options['element'].upper()
    start_date  = options['start_date'] or None
    end_date    = options['end_date'] or None
    method      = options['method']
    cov_str     = options['covariates'] or ''
    power       = float(options['power'])
    npoints     = int(options['npoints'])
    output_base = options['output'] or None
    sample_map  = options['sample'] or None
    min_stn     = int(options['min_stations'])
    flag_t      = flags['t']
    flag_f      = flags['f']
    flag_m      = flags['m']

    if not output_base and not sample_map:
        gs.fatal("Specify output= and/or sample=.")

    require_package('scipy')

    # --- station coordinates ---
    cat_to_xy = get_point_coords(input_map, type='point')
    all_cats  = sorted(cat_to_xy.keys())
    if not all_cats:
        gs.fatal("No point features found in '{}'.".format(input_map))
    gs.message("Loaded {} station locations.".format(len(all_cats)))

    # --- SQLite: open and gather dates ---
    db_path = get_mapset_db()
    conn = sqlite3.connect(db_path)
    cur  = conn.cursor()

    date_sql  = 'SELECT DISTINCT datetime FROM "{}" WHERE element=?'.format(table_name)
    date_args = [element]
    if start_date:
        date_sql += ' AND datetime >= ?'
        date_args.append(start_date)
    if end_date:
        date_sql += ' AND datetime <= ?'
        date_args.append(end_date)
    date_sql += ' ORDER BY datetime'

    try:
        cur.execute(date_sql, date_args)
    except sqlite3.OperationalError as e:
        gs.fatal("Cannot query table '{}': {}".format(table_name, e))

    dates = [r[0] for r in cur.fetchall()]
    if not dates:
        gs.fatal("No '{}' records in '{}' for the requested period.".format(element, table_name))
    gs.message("Found {:,} time steps for element '{}'.".format(len(dates), element))

    # --- grid setup ---
    region   = gs.region()
    grid_xy  = None
    cov_flats = []          # list of flat float32 arrays, one per covariate
    cov_at_grid = None

    if output_base:
        grid_xy = build_grid_xy(region)
        covariates = [c.strip() for c in cov_str.split(',') if c.strip()]
        if method == 'regression' and not covariates:
            gs.warning("method=regression requires covariates=; falling back to spline.")
            method = 'spline'
        for cov_name in covariates:
            arr = read_raster(cov_name, region)
            cov_flats.append(arr.ravel().astype(np.float32))
        if cov_flats:
            cov_at_grid = sample_covariates_at_xy(cov_flats, grid_xy, region)

    # --- sample setup ---
    sample_is_area   = False
    zone_arr         = None   # (rows, cols) int32 cat array for area samples
    sample_pt_cats   = {}     # cat -> (x, y) for point samples

    if sample_map:
        vinfo = gs.parse_command('v.info', map=sample_map, flags='t')
        if int(vinfo.get('areas', 0)) > 0:
            sample_is_area = True
            if not output_base:
                gs.warning(
                    "Area sample= without output=: interpolating at polygon centroids."
                )
            else:
                gs.message("Rasterizing sample polygons...")
                zone_arr = rasterize_to_cats(sample_map, region)
        else:
            sample_pt_cats = get_point_coords(sample_map, type='point')

        if sample_is_area and not output_base:
            # Centroid fallback: get centroid coords
            sample_pt_cats = get_point_coords(sample_map, type='centroid')

    # --- create output tables ---
    error_table  = None
    sample_table = None

    if output_base:
        error_table = '{}_errors'.format(output_base)
        cur.execute('DROP TABLE IF EXISTS "{}"'.format(error_table))
        cur.execute(
            'CREATE TABLE "{}" '
            '(datetime TEXT, element TEXT, rmse REAL, n_stations INTEGER)'.format(error_table)
        )

    if sample_map:
        sample_table = '{}_timeseries'.format(sample_map.split('@')[0])
        cur.execute('DROP TABLE IF EXISTS "{}"'.format(sample_table))
        cur.execute(
            'CREATE TABLE "{}" '
            '(cat INTEGER, datetime TEXT, element TEXT, value REAL)'.format(sample_table)
        )

    conn.commit()

    # --- main interpolation loop ---
    output_maps = []   # (map_name, date_str) for strds registration
    n_skipped   = 0

    for i, date_str in enumerate(dates):
        gs.percent(i, len(dates), 1)

        # query station values for this date
        cur.execute(
            'SELECT cat, value FROM "{}" WHERE datetime=? AND element=?'.format(table_name),
            (date_str, element)
        )
        cat_to_val  = {row[0]: row[1] for row in cur.fetchall()}
        active_cats = set(cat_to_val) & set(all_cats)

        if len(active_cats) < min_stn:
            gs.verbose("Skip {} — {} station(s) < min_stations={}.".format(
                date_str, len(active_cats), min_stn))
            n_skipped += 1
            continue

        active_list = sorted(active_cats)
        stn_xy = np.array([cat_to_xy[c] for c in active_list], dtype=np.float64)
        vals   = np.array([cat_to_val[c] for c in active_list], dtype=np.float64)

        # covariate values at (unfilled) station locations
        cov_at_stn = None
        if cov_flats and method in ('regression', 'spline'):
            cov_at_stn = sample_covariates_at_xy(cov_flats, stn_xy, region)

        # LOO error metric — computed before fill so it reflects real data density
        if error_table:
            rmse = loo_rmse(method, stn_xy, vals, power, npoints, cov_at_stn)
            cur.execute(
                'INSERT INTO "{}" VALUES (?, ?, ?, ?)'.format(error_table),
                (date_str, element,
                 None if math.isnan(rmse) else rmse,
                 len(active_list))
            )

        # neighbor fill
        if flag_f:
            cat_to_val, active_cats = neighbor_fill(
                all_cats, cat_to_xy, cat_to_val, active_cats)
            active_list = sorted(active_cats)
            stn_xy = np.array([cat_to_xy[c] for c in active_list], dtype=np.float64)
            vals   = np.array([cat_to_val[c] for c in active_list], dtype=np.float64)
            if cov_flats and method in ('regression',):
                cov_at_stn = sample_covariates_at_xy(cov_flats, stn_xy, region)

        # convex hull mask (lazily built per date since station set varies)
        hull = None
        if flag_m and output_base and len(stn_xy) >= 3:
            from scipy.spatial import Delaunay
            hull = Delaunay(stn_xy)

        # --- grid output ---
        if output_base and grid_xy is not None:
            interp_flat = interpolate(
                method, stn_xy, vals, grid_xy,
                power, npoints, cov_at_stn, cov_at_grid
            )
            if hull is not None:
                outside = hull.find_simplex(grid_xy) < 0
                interp_flat[outside] = np.nan

            map_name = '{}_{}'.format(output_base, date_str.replace('-', ''))
            write_raster(interp_flat.reshape(region['rows'], region['cols']),
                         map_name, region)
            output_maps.append((map_name, date_str))

            # area-weighted extraction from interpolated grid
            if sample_map and sample_is_area and zone_arr is not None:
                rows_buf = []
                for cat_val in np.unique(zone_arr):
                    if cat_val == 0:
                        continue
                    cell_vals = interp_flat[zone_arr.ravel() == cat_val]
                    cell_vals = cell_vals[~np.isnan(cell_vals)]
                    if len(cell_vals):
                        rows_buf.append((int(cat_val), date_str, element,
                                         float(np.mean(cell_vals))))
                if rows_buf:
                    cur.executemany(
                        'INSERT INTO "{}" VALUES (?, ?, ?, ?)'.format(sample_table),
                        rows_buf
                    )

            # point extraction from grid
            elif sample_map and sample_pt_cats and not sample_is_area:
                rows, cols = region['rows'], region['cols']
                interp_2d = interp_flat.reshape(rows, cols)
                rows_buf = []
                for cat_val, (sx, sy) in sample_pt_cats.items():
                    ci = int((sx - region['w']) / region['ewres'])
                    ri = int((region['n'] - sy) / region['nsres'])
                    if 0 <= ri < rows and 0 <= ci < cols:
                        v = interp_2d[ri, ci]
                        if not math.isnan(v):
                            rows_buf.append((cat_val, date_str, element, float(v)))
                if rows_buf:
                    cur.executemany(
                        'INSERT INTO "{}" VALUES (?, ?, ?, ?)'.format(sample_table),
                        rows_buf
                    )

        # --- sample-only (no grid output) ---
        elif sample_map and sample_pt_cats:
            pt_list = sorted(sample_pt_cats)
            pt_xy   = np.array([sample_pt_cats[c] for c in pt_list], dtype=np.float64)
            cov_pt  = sample_covariates_at_xy(cov_flats, pt_xy, region) if cov_flats else None
            preds   = interpolate(method, stn_xy, vals, pt_xy,
                                  power, npoints, cov_at_stn, cov_pt)
            rows_buf = [
                (c, date_str, element, float(v))
                for c, v in zip(pt_list, preds)
                if not math.isnan(v)
            ]
            if rows_buf:
                cur.executemany(
                    'INSERT INTO "{}" VALUES (?, ?, ?, ?)'.format(sample_table),
                    rows_buf
                )

        conn.commit()

    gs.percent(len(dates), len(dates), 1)
    conn.close()

    # --- strds registration ---
    if flag_t and output_maps:
        register_strds(output_base, element, output_maps)

    # --- summary ---
    n_done = len(dates) - n_skipped
    gs.message(
        "Done: {:,} time steps interpolated, {:,} skipped "
        "(< {} stations).".format(n_done, n_skipped, min_stn)
    )
    if output_base and output_maps:
        gs.message("Output rasters: {}_<date> ({} maps).".format(
            output_base, len(output_maps)))
        if error_table:
            gs.message("Error metrics: table '{}' in mapset SQLite db.".format(error_table))
    if sample_map and n_done > 0:
        gs.message("Sample time series: table '{}' in mapset SQLite db.".format(sample_table))
        gs.message("  db.select sql=\"SELECT * FROM {} LIMIT 10\"".format(sample_table))


if __name__ == '__main__':
    main()
