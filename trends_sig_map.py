#!/usr/bin/env python
import numpy as np
from scipy import stats
import cartopy
import matplotlib.pyplot as plt
#import iris
#from eofs.iris import Eof
#import cm

# %% puting files in an ascending order
flist = [f'Indice_com_files/Indices/PRCPTOT/RR_GH_{k}_PRCPTOT' for k in range(1,561)]

'''
prec = []
for f in flist:
	date=np.transpose(sp.genfromtxt(str(f),usecols=[0],skip_header=1))
	data=np.transpose(sp.genfromtxt(str(f),usecols=[1,2,3,4,5,6,7,8,9,10,11,12],skip_header=1))
	DJF=np.transpose(sp.genfromtxt(str(f),usecols=[9,10,11],skip_header=1))
	djf = np.average(DJF,axis=0)
	#col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12= data
	sum_all=np.sum(data,axis=0)
	prec.append(sum_all)
	pp=np.transpose(prec)
'''
	#np.savetxt('rx5day.txt',np.c_[np.transpose(prec)])

#flist = glob.glob(path+'/Try/Indices/PRCPTOT/'+'*')
# %% annual indices

prec = []
for f in flist:
	data = np.genfromtxt(str(f),usecols=[0,1],skip_header=1).T
	date,ind = data



	prec.append(ind)
	#ppt=np.round("%.2f" % prec)
#	pp=np.transpose(prec)

	#print prec

#sst = iris.load(pp)
#solver = Eof(sst, weights='coslat')

#sst = iris.load_cube(prec)

#dat= np.loadtxt('CDD.txt')




'''
lon = np.arange(-3.5,2,0.5)
lat = np.arange(4.5,12,0.5)
'''


m = Basemap(projection='merc',llcrnrlon=-3.5,llcrnrlat=4.5,urcrnrlon=1.5,urcrnrlat=11.5,resolution='i')
m.drawcoastlines(linewidth=2)
m.drawstates()
m.drawcountries(linewidth=2)
m.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # can use HTML names or codes for colors
#m.drawcounties()
parallels = np.arange(4.5,11.5,1) # make latitude lines ever 5 degrees from 30N-50N
meridians = np.arange(-3.5,1.5,1) # make longitude lines every 5 degrees from 95W to 70W
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12,linewidth=2.5)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12,linewidth=2.5)

lon = np.arange(-3.375, 1.625, 0.25)
lat = np.arange(4.625,11.625,0.25)



x,y= np.meshgrid(lon,lat)
a,b = m(x,y)

n = np.size(lon)
mn = np.size(lat)
t = np.size(date)

X = np.reshape(pp,(t,mn,n))
X2 = np.average(X,axis=0)


x = np.arange(1,36,1)

#np.savetxt('ind.txt',X)
#plt.show()

slopes = []
r_values =[]
p_values= []
stdev=[]
intecep = []
####Statistics computation
for i in range(len(prec)):




	slope,intercept,r_value,p_value,std_err=stats.linregress(x,prec[i])
	#t_stat, p_val = stats.ttest_ind(x, prec[i], nan_policy='omit')
	r_values.append(r_value)
	p_values.append(p_value)
	stdev.append(std_err)
	intecep.append(intercept)
	slopes.append(slope)


#sl,inte,r,p,s=stats.linregress(slopes[:])


clevs = np.arange(-0.8,0.8,0.15)
#clevs = np.arange(-1,1.5,0.5)

slop = np.reshape(slopes,(mn,n))

lim = np.arange(1,561)

#slope,intercept,r_value,p_value,std_err=stats.linregress(lim,slopes)

sslop=np.str_(slopes)
ppval=np.str_(p_values)
np.savetxt('slop_pvalue.txt',np.c_[np.float_(slopes),np.float_(p_values)])

sig=[]
non_sig =[]

non_sig_pos=[]
non_sig_neg=[]

sig_pos=[]
sig_neg=[]




data = np.loadtxt('slop_pvalue.txt')
sloppe = data[:,0]
pval=data[:,1]


for i in range(len(sloppe)):

	if pval[i]<=0.1:
		sig.append(sloppe[i])
		#sloppe[i]=sloppe[i]

		ppsig=(len(sig)/560.0)*100

		if sloppe[i]>0:
			sig_pos.append(sloppe[i])
			psigpos=float(len(sig_pos))/float(len(sig))*ppsig

		if sloppe[i]<0:
			sig_neg.append(sloppe[i])
			psigneg=float(len(sig_neg))/float(len(sig))*ppsig


	if pval[i]>0.1:
		#sloppe[i]=np.nan

		non_sig.append(sloppe[i])
		ppnonsig=(len(non_sig)/560.0)*100

		if sloppe[i]>0:
			non_sig_pos.append(sloppe[i])
			pnonsigpos=(float(len(non_sig_pos))/float(len(non_sig)))*ppnonsig





cc = np.reshape(sloppe,(mn,n))


#print psigpos,(ppsig-psigpos),(ppnonsig-pnonsigpos),len(sig),len(non_sig)
print(psigpos,(ppsig-psigpos),pnonsigpos,(ppnonsig-pnonsigpos),len(sig),len(non_sig))






slop = np.reshape(slopes,(mn,n))
slop[slop==0.0]=np.nan


#s_masked = np.ma.array(slop,mask=np.isnan(slop))

#cmap = plt.cm.OrRd
#cmap=plt.cm.RdGy
#cmap=plt.cm.Accent
#cmap.set_over(color='black')
cmap='gist_rainbow'
#cmap='gist_ncar'



ppt = m.contourf(a,b,cc,cmap=cmap)
#plt.colorbar(ppt)

cb = m.colorbar(ppt,"right", size="10%", pad="7%")


#plt.title('()',fontsize=20)

plt.tight_layout()
plt.savefig('Figures/trends/trends_90_95/90/cwd.png',dpi=300,facecolor='w')
plt.show()



####check sig_pvalues


#sig_point.append(p_)'''

