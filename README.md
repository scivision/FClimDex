# FClimDex
http://etccdi.pacificclimate.org/software.shtml


## Prereqs
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
