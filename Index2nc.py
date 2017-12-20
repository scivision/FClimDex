#!/usr/bin/env python
"""convert Fortran FClimdex output text files indices to single NetCDF4 file

currently handles:
CDD, CSDI, CWD, RX5day

python Index2nc.py data/CDD data/data.nc test.nc

python Index2nc.py data/RX5day data/data.nc test.nc
"""
from pathlib import Path
import numpy as np
import pandas
from netCDF4 import Dataset
import xarray
from datetime import datetime
import calendar

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
# selecting Ghana
    ilat = (lat >= 4.5 ) & (lat <= 11.5)
    ilon = (lon >= -3.5) & (lon <= 1.5)

    lat = lat[ilat]
    lon = lon[ilon]
# %% read time from first file--assumes all files have same time span!
    dat = _getdat(flist[0],True)

    time = _gettimes(flist[0], dat)
# %% setup output array
    nc = xarray.DataArray(data=np.empty((len(time),lat.size,lon.size)),
                          coords={'time':time,'lat':lat,'lon':lon},
                          dims=['time','lat','lon'])
# %% wrangle input text files
    # TODO: this unraveling is a total guess. Which file cooresponds to which coordinate?
    iu = np.unravel_index(range(len(flist)),
                          (lat.size, lon.size),
                          order='F')

    for i,f in enumerate(flist):
        dat = _getdat(f)
        nc[:,iu[0][i],iu[1][i]] = dat
# %% write NetCDF4 output
    nc.to_netcdf(ofn)


def _getdat(fn:Path, init:bool=False):
    tail = fn.name.split('_')[-1]

    if tail in ('CDD','CSDI','CWD'):
        dat = pandas.read_csv(fn, sep='\s+', index_col=0).squeeze()
    elif tail in ('DTR','RX5day'): # FIXME: leaving off "annual" for now
        dat = pandas.read_csv(fn, sep='\s+',
                              index_col=0, usecols=range(13))
        if not init:
            dat = dat.values.ravel()  # vectorize in time order
    else:
        raise ValueError(f'unknown filetype {fn}')

    return dat

def _gettimes(fn:Path, dat):
    tail = fn.name.split('_')[-1]

    if tail in ('CDD','CSDI','CWD'):
        time = [datetime(year=i,month=12,day=31) for i in dat.index]
    elif tail in ('DTR','RX5day'):
        years = dat.index
        time = []
        for y in years:
            for m in range(1,13):
                time.append(datetime(year=y,month=m,day=calendar.mdays[m]))
    else:
        raise ValueError(f'unknown filetype {fn}')

    return time


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('path',help='path containing files to convert')
    p.add_argument('coordfile',help='NetCDF4 file containing coordinates to READ')
    p.add_argument('ofn',help='NetCDF4 filename to write')
    p.add_argument('-g','--glob',help='glob pattern to convert',default='*')
    p = p.parse_args()

    index2nc(p.path, p.glob, p.ofn, p.coordfile)