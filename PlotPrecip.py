#!/usr/bin/env python
"""
http://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.daily-improved/.global/.0p25/.prcp/X/%28-20%29/%2815%29/RANGE/Y/%280%29/%2830%29/RANGE/T/%281995%29/%282015%29/RANGEEDGES/datafiles.htmlsudo

apt install libgeos-dev
pip install https://github.com/matplotlib/basemap/archive/v1.1.0.tar.gz
"""
from datetime import datetime
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.basemap import Basemap

in_file='data/data.nc'

# %% Read in netcdf data
with Dataset(in_file,'r') as f:
    lon  = f['X'][:]
    lat  = f['Y'][:]
# %% Figure setting
fig = plt.figure(figsize=(10,8))
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
with Dataset(in_file,'r') as f:
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
    plt.draw()
    plt.pause(0.05)
    for h in hc.collections:
        h.remove()  # else it gets slow and uses tons of RAM


raise SystemExit
#x,y= np.meshgrid(lo1,la1 )

#a,b = m(x,y)


#ppt = m.contourf(a,b,prec[0,:],clevs)

####selection of a particular grid point
####creating lists to append all precipitaion values from grids
#sizes = (7670, 560)
#list_grids = np.zeros(sizes)


#### Looping through each grid to produce 560 grid boxes of precipittaion data
#rains=[]
#for i in range(lat.size):
#	for j in range(lon.size):
#		 rains.append(prec[1:,i,j])


####Creating the year, month, day data for the file

yr = []
mon = []
day = []


beg = datetime(1981, 1, 1)
end = datetime(2015, 12, 31)
t_step = datetime.timedelta(days=1)

###################################################Mke this work

while beg < end:
    yr.append(beg.strftime('%Y'))
    mon.append(beg.strftime('%m'))
    day.append(beg.strftime('%d'))
    beg += t_step

######## creation of arrays for dates, minT, maxT data
dates = np.array([yr,mon,day]) ####dates in array format
minT  = [-99.9]*np.size(dates[0]) ###an array of min tempertature
maxT  = [-99.9]*np.size(dates[0]) ###an array of max tempertature

######## 35 years of data
'''
###### Saving files ready for R code on various indices
for i in range(len(rains)):
	print('Done with RR_GH_'+str(i+1))

	np.savetxt('Out_files/RR_GH_'+str(i),np.c_[dates[0],dates[1],dates[2],rains[i],minT,maxT],fmt='%s')

print('Completed for all grids')

'''

############## Reading file one after the other to be executed in R
'''
path=os.getcwd()
infiles=glob.glob(os.path.join('Out_files/RR_GH_0'))
out_file=open('1.txt','w')
for filename in infiles:


	with open(filename, 'r') as f:
		for line in f:
			line=line.split()

			#a = format(line[3],'.2f')
			a = round(float(line[3]),1)


			out_file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\r\n')


'''
#####rclimdex1.1_131115.r


#np.savetxt('Out/out_'+basename(file),out,fmt=['%4.4d','%2.2d','%2.2d','%.1f'],delimiter=',')




		#list_grids.append(m)

		#plt.plot(m,'r',linewidth=2)

#plt.show()



'''
with open('testdata.txt', 'r') as input:
    for (counter,line) in enumerate(input):
        with open('filename{0}.txt'.format(counter), 'w') as output:
            output.write(line)
'''
'''
x,y= np.meshgrid(lo1,la1 )

a,b = m(x,y)


ppt = m.contourf(a,b,m,clevs)
'''


'''
i =0

while (i <len(lo1)):



count=0
for lo in xlo1:
	for la in la1:
		count +=1
		#print lo,la,count

		a=prec[:,2,5]

'''
