# v.interp.timeseries — Claude Code context

## What this module does

GRASS GIS addon that interpolates a station-based climate time series (written by
v.in.ghcn) to raster grids, to point samples, or to polygon area-weighted averages.

**Primary inputs:**
- Station vector map with point geometry (from v.in.ghcn)
- Companion SQLite table `{input}_timeseries` with columns `cat, datetime, element, value`

**Outputs (at least one required):**
- `output=`: one FCELL raster per time step named `{output}_YYYYMMDD` (or `_YYYYMM01` for monthly); optionally registered as a strds with `-t`
- `sample=`: SQLite table `{sample}_timeseries (cat, datetime, element, value)` for point or polygon extraction
- `{output}_errors` table: LOO cross-validation RMSE per time step

## Interpolation methods

| method= | Algorithm | Notes |
|---------|-----------|-------|
| idw | Inverse-distance weighting | scipy.spatial.cKDTree; power= and npoints= |
| spline | Thin-plate spline | scipy.interpolate.RBFInterpolator(kernel='thin_plate_spline'); exact interpolator |
| regression | OLS on covariates= + spline residuals | numpy.linalg.lstsq; analogous to PRISM detrending |

## Key implementation details

- `get_point_coords()`: `v.out.ascii format=point type=point|centroid separator=pipe` → cat→(x,y) dict
- `build_grid_xy()`: meshgrid of cell centres from `gs.region()`; north→south row order
- `read_raster()`: `r.out.bin -f bytes=4 null=-9999` → float32 numpy array; NaN for null
- `write_raster()`: float32 binary → `r.in.bin -f bytes=4 anull=-9999 north= south= east= west= rows= cols=`
- `rasterize_to_cats()`: `v.to.rast use=cat` → CELL raster → `r.mapcalc float()` → `r.out.bin -f` → int32 numpy
- `sample_covariates_at_xy()`: nearest-cell sampling of flat covariate arrays
- `loo_rmse()`: LOO cross-validation; skips if `method=spline` and `n > 50` (O(n³) per fold)
- `neighbor_fill()`: cKDTree nearest-neighbor fill for stations missing on a given date; fills before grid write, after LOO error metric
- Zone averaging: `zone_arr.ravel() == cat_val` mask on flat interpolated array; equal-weight cell mean
- Progress: `gs.percent(i, n_dates, 1)` inside date loop

## Convex hull mask (-m)

Uses `scipy.spatial.Delaunay.find_simplex(grid_xy) < 0` to identify points outside the convex hull of active stations. Recomputed per date (station set varies). Applied after interpolation, before raster write.

## Raster I/O byte order

- Writes/reads IEEE 754 float32 in native (little-endian on x86-64) byte order
- GRASS null → -9999.0 sentinel in binary; `anull=-9999` or `null=-9999` in GRASS commands

## strds registration (-t)

- `t.create type=strds temporaltype=absolute output={output_base} ...`
- `t.register input={output_base} file={reg_file}` where reg_file lists `mapname|start|end` per line
- Daily: end = start + 1 day; monthly: end = first of next month

## Error table

`{output}_errors (datetime TEXT, element TEXT, rmse REAL, n_stations INTEGER)`
- Created fresh each run (DROP TABLE IF EXISTS + CREATE TABLE)
- LOO RMSE computed before neighbor fill; NULL stored for spline with n > 50

## Sample table

`{sample}_timeseries (cat INTEGER, datetime TEXT, element TEXT, value REAL)`
- Created fresh each run
- Area samples: equal-weight mean of grid cells within polygon (requires output=)
- Point samples: nearest-cell from grid (or direct interpolation at coordinates if no output=)
- Area-without-output fallback: interpolates at polygon centroids via `type=centroid` in v.out.ascii

## Repo and status

- GitHub: `https://github.com/awickert/v.interp.timeseries`
- Local repo: `/home/awickert/dataanalysis/v.interp.timeseries`
- **Not yet submitted to GRASS addons** (as of June 2026, initial implementation)

## Broader context

Third module in the GIS-native hydrological forcings layer, after v.in.ghcn (station download)
and r.in.gridmet/r.in.prism/r.gpm.imerg (gridded products). Intended workflow:

```
g.region watershed=my_watershed
v.in.ghcn output=gauges elements=PRCP,TMAX,TMIN start_date=1980-01-01
v.interp.timeseries input=gauges element=PRCP method=regression \
    covariates=dem output=prcp sample=watersheds -t -m
```

The `sample=` output feeds directly into GSFLOW-PRMS as daily forcing per HRU.
