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
#from matplotlib import gridspec
from mpl_toolkits.basemap import Basemap

def plotprecip(fn:Path):
    fn = Path(fn).expanduser()
# %% Read in netcdf data
    with Dataset(str(fn), 'r') as f:
        lon  = f['X'][:]
        lat  = f['Y'][:]
    # %% Figure setting
    fig = figure(figsize=(10,8))
    ax = fig.gca()

    spc=0.5
    lat_beg=0
    lat_end=15
    lon_beg=-10
    lon_end=5

    m = Basemap(projection='merc',llcrnrlon=lon_beg,llcrnrlat=lat_beg,urcrnrlon=lon_end,urcrnrlat=lat_end,resolution='l')
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # can use HTML names or codes for colors

    parallels = np.arange(lat_beg,lat_end,spc) # make latitude lines ever 5 degrees from 30N-50N
    meridians = np.arange(lon_beg,lon_end,spc) # make longitude lines every 5 degrees from 95W to 70W
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize='small')
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize='small')

    ax.legend(loc='upper center', shadow='True', fontsize='normal')
    ht = ax.set_title('')
    clvl= np.arange(0, 6, 0.05)
    #ppt = m.contourf(a,b,prec[0,:],clevs)
    # %% Selection of an area of interest
    ilat = (lat >= 4.5 ) & (lat <= 11.5)
    ilon = (lon >= -3.5) & (lon <= 1.5)
    lat = lat[ilat]
    lon = lon[ilon]

    x,y= np.meshgrid(lon, lat)
    a,b = m(x,y)
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
        hc=m.contourf(a,b,p,clvl)
        ht.set_text(str(t))
        draw()
        pause(0.05)
        for h in hc.collections:
            h.remove()  # else it gets slow and uses tons of RAM


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('fn',help='netcdf file to plot')
    p = p.parse_args()

    plotprecip(p.fn)
