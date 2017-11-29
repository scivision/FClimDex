#!/usr/bin/env python

def old():
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
