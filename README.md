[![Build Status](https://travis-ci.org/scivision/FClimDex.svg?branch=master)](https://travis-ci.org/scivision/FClimDex)

# FClimDex
http://etccdi.pacificclimate.org/software.shtml

Plots NetCDF4 climate data on animated map.
Processes climate text files with Python calling Fortran FClimDex on all files in a directory. 

## Install

### Prereqs
Python >= 3.6

```sh
apt install gfortran
```

### Build
```sh
python -m pip install -e .
```

## Programs

```sh
python RunFClimdex.py data/
```
processes each text `*.dat` file in directory `data/` into 33 output files, in a unique directory for each input file

```sh
python index2nc.py data/ out.nc
```
Reads the text data files from the Fortran program FClimdex (as part of RunFClimdex.py) and converts them to a single NetCDF4 file (3-D: time, lat, lon).

```sh
python PlotPrecip.py
```
Takes NetCDF4 input files and plots a movie of precipitation from them, overlaid on a colorful matplotlib map.


## Notes

If you get error

> PermissionError: [Errno 13] Permission denied: './FClimDex'

then as usual
```sh
chmod +x FClimDex
```
