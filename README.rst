.. image:: https://travis-ci.org/scivision/FClimDex.svg?branch=master)
   :target: https://travis-ci.org/scivision/FClimDex

========
FClimDex
========
http://etccdi.pacificclimate.org/software.shtml

Plots NetCDF4 climate data on animated map.
Processes climate text files with Python calling Fortran FClimDex on all files in a directory. 

Install
=======

Prereqs
-------
Python >= 3.6

::

    apt install gfortran


Build
-----

::

    python -m pip install -e .


Usage
=====
Process each text `*.dat` file in directory `data/` into 33 output files, in a unique directory for each input file::

    python RunFClimdex.py data/

Reads the text data files from the Fortran program FClimdex (as part of RunFClimdex.py) and converts them to a single NetCDF4 file (3-D: time, lat, lon)::

    python index2nc.py data/ out.nc

Takes NetCDF4 input files and plots a movie of precipitation from them, overlaid on a colorful map::

    python PlotPrecip.py



Notes
=====

If you get error

    PermissionError: [Errno 13] Permission denied: './FClimDex'

then as usual

.. code:: bash

    chmod +x FClimDex

