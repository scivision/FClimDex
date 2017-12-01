# FClimDex
http://etccdi.pacificclimate.org/software.shtml

Plots NetCDF4 climate data on animated map.
Processes climate text files with Python calling Fortran FClimDex on all files in a directory. 

## Prereqs
Python >= 3.6

```sh
apt install gfortran libgeos-dev
pip install https://github.com/matplotlib/basemap/archive/v1.1.0.tar.gz
```

## Programs

```sh
python RunPrecip.py data
```
will process each text `*.dat` file into 33 output files, in a unique directory for each input file

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
