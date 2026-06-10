MODULE_TOPDIR = $(shell grass --config path)

PGM = v.interp.timeseries

include $(MODULE_TOPDIR)/include/Make/Script.make

default: script
