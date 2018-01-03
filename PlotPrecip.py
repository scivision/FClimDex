#!/usr/bin/env python
"""
E.g. get coordinates .nc donor file from:
http://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.daily-improved/.global/.0p25/.prcp/X/%28-20%29/%2815%29/RANGE/Y/%280%29/%2830%29/RANGE/T/%281995%29/%282015%29/RANGEEDGES/datafiles.html

Example: movie plot over West Africa
python PlotPrecip.py data/data.nc
"""
from pathlib import Path
from datetime import datetime
import numpy as np
from netCDF4 import Dataset
from matplotlib.pyplot import figure,draw,pause
import cartopy

PROJ = cartopy.crs.PlateCarree()  # arbitrary

def plotprecip(fn:Path):
    fn = Path(fn).expanduser()
# %% Read in netcdf data
    with Dataset(str(fn), 'r') as f:
        lon  = f['X'][:]
        lat  = f['Y'][:]
    # %% Figure setting
    fig = figure(figsize=(10,8))
    ax = fig.gca(projection=PROJ)

    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')

    ax.set_extent((-10,5,0,15))

#    ax.legend(loc='upper center', shadow='True')
    ht = ax.set_title('')
    clvl= np.arange(0, 6, 0.05)
    #ppt = m.contourf(a,b,prec[0,:],clevs)
    # %% Selection of an area of interest
    ilat = (lat >= 4.5 ) & (lat <= 11.5)
    ilon = (lon >= -3.5) & (lon <= 1.5)
    lat = lat[ilat]
    lon = lon[ilon]

    lon,lat = np.meshgrid(lon, lat)
    # %% precipitation for selected area above
    with Dataset(str(fn),'r') as f:
        prec = f['prcp'][:,ilat, ilon]
        jtime = f['T'][:].astype(int)  # data is in one day increments

    # TODO improvised date conversion. Is it correct?
    time = []
    for t in jtime:
        y = f'{int(t/365.25)-4713:04d}'
        j = f'{int(t%365.25)+1:03d}'
        time.append(datetime.strptime(y+j,'%Y%j'))

    # %% animation
    for p,t in zip(prec,time):
        hc = ax.contourf(lon, lat, p, clvl,
                         transform=PROJ,
                         zorder=2,
                         cmap='cubehelix_r',
                         alpha=0.5)
        ht.set_text(str(t))
        draw()
        pause(1)
        for h in hc.collections:
            h.remove()  # else it gets slow and uses tons of RAM


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('fn',help='netcdf file to plot')
    p = p.parse_args()

    plotprecip(p.fn)
