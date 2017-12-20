#!/usr/bin/env python
"""convert Fortran FClimdex output text files indices to single NetCDF4 file"""
from pathlib import Path
import numpy as np
import pandas
from netCDF4 import Dataset
import xarray
import datetime

def index2nc(path:Path, glob:str, ofn:Path, cfn:Path):
    cfn = Path(cfn).expanduser()
    ofn = Path(ofn).expanduser()
# %% initial pass on input text files
    path = Path(path).expanduser()
    flist = path.glob(glob)
    flist = [f for f in flist if f.is_file()] # eliminate directories
# %% load coordinates from NetCDF4 file
    with Dataset(str(cfn),'r') as f:
        lon  = f['X'][:]
        lat  = f['Y'][:]
# %% read time from first file--assumes all files have same time span!
    dat = pandas.read_csv(flist[0], sep='\s+', index_col=0)
    time = [datetime.date(year=i,month=12,day=31) for i in dat.index]

    nc = xarray.DataArray(data=np.empty((dat.shape[0],lat.size,lon.size)),
                          coords={'time':time,'lat':lat,'lon':lon},
                          dims=['time','lat','lon'])

# %% wrangle input text files

    for f in flist:
        if f.name.endswith('_CDD'):
            cdd


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('path',help='path containing files to convert')
    p.add_argument('coordfile',help='NetCDF4 file containing coordinates to READ')
    p.add_argument('ofn',help='NetCDF4 filename to write')
    p.add_argument('-g','--glob',help='glob pattern to convert',default='*')
    p = p.parse_args()

    index2nc(p.path, p.glob, p.ofn, p.coordfile)