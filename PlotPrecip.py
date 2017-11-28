########## Import necessary libraries
import numpy as np
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.basemap import Basemap
import datetime
import os
import io


in_file='data.nc'
	
#### Reading in the variables of the netcdf data
nc_fil = Dataset(in_file,'r')
lons= nc_fil.variables['X'][:]
lats= nc_fil.variables['Y'][:]
prec = nc_fil.variables['prcp'][:]
time = nc_fil.variables['T'][:]

#### Figure setting 
fig = plt.figure(figsize=(10,8))

#lon_range=np.arange(-3.5,1.5,0.5)

#lat_range=np.arange(4.5,11.5,0.5)

'''
###command to give long lat ranges
lon_beg = input('Please enter your lowest longitude: ')
lon_end = input('Please enter your highest longitude: ')

lat_beg = input('Please enter your lowest latitude: ')
lat_end = input('Please enter your highest latitude: ')

'''

###indicate spacing of interest for meridians and parallels
#spc = input('Please indicate the grid spacing of your interest: ')

#print('This will only take a second')
	
#m = Basemap(projection='merc',llcrnrlon=lon_beg,llcrnrlat=lat_beg,urcrnrlon=lon_end,urcrnrlat=lat_end,resolution='l')

spc=0.5
lat_beg=0
lat_end=30
lon_beg=-20
lon_end=15

m = Basemap(projection='merc',llcrnrlon=lon_beg,llcrnrlat=lat_beg,urcrnrlon=lon_end,urcrnrlat=lat_end,resolution='l')
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # can use HTML names or codes for colors
#m.drawcounties()
parallels = np.arange(lat_beg,lat_end,spc) # make latitude lines ever 5 degrees from 30N-50N
meridians = np.arange(lon_beg,lon_end,spc) # make longitude lines every 5 degrees from 95W to 70W
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=5.5)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=5.5)

x,y= np.meshgrid(lons,lats)

a,b = m(x,y)

plt.legend(loc='upper center', shadow='True', fontsize='14')

clevs = np.arange(0,6,0.05)


#ppt = m.contourf(a,b,prec[0,:],clevs)

####Selection of an area of interest 


idx_latarea = np.where((lats >= 4.5 ) & (lats <= 11.5))[0]
idx_lonarea = np.where((lons >= -3.5) & (lons <= 1.5))[0]
la1 = lats[idx_latarea[0]:idx_latarea[-1]+1]
lo1 = lons[idx_lonarea[0]:idx_lonarea[-1]+1]

##############precipitation for selected area above

prec = prec[:,idx_latarea[0]:idx_latarea[-1]+1, idx_lonarea[0]:idx_lonarea[-1]+1]

#x,y= np.meshgrid(lo1,la1 )

#a,b = m(x,y)


#ppt = m.contourf(a,b,prec[0,:],clevs)

####selection of a particular grid point
####creating lists to append all precipitaion values from grids
#sizes = (7670, 560)
#list_grids = np.zeros(sizes)


#### Looping through each grid to produce 560 grid boxes of precipittaion data
rains=[] 
for i in range(len(la1)):
	for j in range(len(lo1)):
		 rains.append(prec[1:,i,j])


####Creating the year, month, day data for the file

yr = []
mon = []
day = []


beg = datetime.datetime(1981, 01, 01)
end = datetime.datetime(2015, 12, 31)
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
