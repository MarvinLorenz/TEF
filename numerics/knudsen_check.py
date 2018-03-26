# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 14:35:57 2018

@author: lorenz
"""
import numpy as np
import sys

def convert_to_knudsen_check_format(salt,velocity,height,ds,out_path,save_name):
	
	##open the file we want to write in
	outfile = open(out_path+save_name,'w')
	outfile.write(str(np.shape(salt)[0]) + '\t' + str(np.shape(salt)[1]*np.shape(salt)[2]) + '\t' + str(np.float64(np.min(salt))) + '\t' + str(np.float64(np.max(salt))) + '\n')	
	print(np.shape(salt))
	t=0
	while t < np.shape(salt)[0]:
		counter = 1
		k=0
		while k<np.shape(salt)[1]:
			i=0
			while i < np.shape(salt)[2]:
				outfile.write(str(t+1) + '\t' + str(counter) + '\t' + str(np.float64(velocity[t,k,i]))+ '\t' + str(np.float64(salt[t,k,i])) + '\t' + str(np.float64(height[t,k,i]*ds[i])) + '\n')
				i+=1
				counter +=1
			k+=1
		t+=1
	print('Convertion done')
	return('jop')

"""
!     This program calculates exchange flow bulk properties 
!     (the inflow and outflow volume fluxes and salinities) 
!     across hydrographic transects across estuarine channels.
!
!     As input, data for inflow velocity normal to the transect, u,
!     salinity, s, and area a of discrete interfaces along the transect
!     need to be given.
!
!     The format of the input file is as follows:
!     First line is a header line giving # of time steps, 
!     # of spatial increments per time step, minimum salinity and 
!     maximum salinity.
!     All other lines include the following 5 entries:
!        1. number of time step
!        2. number of spatial increment 
!        3. inward velocity normal to the transect [m/s]
!        4. salinity [g/kg]
!        5. area of grid cell [m**2]
!
!     First, the volume and salinity fluxes per salinity class are calculated 
!     for a  maximum of salinity classes (here: nnmax = 2**23). Then, these
!     volume and salinity fluxes per salinity class are subsequently subsampled 
!     into fewer and larger salinity classes.
!
!     As output, the inflow and outflow volume fluxes and salinities are
!     given to help estimating an optimal number of salinity classes 
!     for the calculation of the Total Exchange Flow (TEF) and the
!     derived bulk properties.
"""
def knudsen_check(path_to_file,path_to_savefile):

	t,i,u,s,a = np.loadtxt(path_to_file, unpack=True, skiprows=1,dtype=np.float64) #load all data and split them into one variable for each column
	kmax,imax,smin,smax = np.genfromtxt(path_to_file, unpack=True, skip_footer =len(t),dtype=np.float64)# read the header and store in 4 variables
	mmax = 23
	nmax = 2**mmax
	ds=np.float64((smax-smin)/float(nmax))
	qv = np.zeros(shape=(nmax+1,),dtype=np.float64) # Initialise the variable qv with length of salinity classes with zeros;+1 so that the first class is zero and the indexing is the same as in the paper
	qs = np.zeros(shape=(nmax+1,),dtype=np.float64) # Initialise the variable qs with length of salinity classes with zeros 
	
	for tt in range (0,int(kmax)): #time loop
		for ii in range (0,int(imax)): #cross section loop
			if not ii+1 == i[ii+tt*imax]:
				sys.exit('wrong i-indexing in input file')
			if not tt+1 == t[ii+tt*imax]:
				sys.exit('wrong k-indexing in input file')
			if s[tt*imax+ii] > smax:
				sys.exit('s is greater than smax')
			if s[tt*imax+ii] < smin:
				sys.exit('s is smaller than smin')
			idx=int((s[tt*imax+ii]-smin)/(smax-smin)*(nmax))+1 #computation of the index
			if s[tt*imax+ii] == smax:
				idx = nmax
			qv[idx]=qv[idx]+u[tt*imax+ii]*a[tt*imax+ii]/kmax/ds # Adding to high-res volume flux
			qs[idx]=qs[idx]+u[tt*imax+ii]*s[tt*imax+ii]*a[tt*imax+ii]/kmax/ds # Adding high-res salinity flux
	
	"""
	      ! At this stage, qv and qs contain volume and salinity fluxes 
	      ! at very high resolution with nmax=2**mmax salinity classes
	      ! In the next step we apply coarse-graining to check the
	      ! dependence of Qin, Qout, sin and sout 
	      ! on the number of salinity classes.
	"""
	qv = qv[1:] #slice the first class away since it was only needed for correct indexing
	qs = qs[1:]
	print('#classes    inflow sal.    outflow sal.    inflow vol.    outflow vol.')
	output=open(path_to_savefile,'w')
	output.write('#classes \t inflow sal \t outflow sal. \t inflow vol. \t outflow vol. \n')
	for mm in range (1,mmax+1):#number of experiments mm
		nmax=2**mm
		jmax=2**mmax/nmax # number of low-res per high-res salinity classes
		dss = np.float64((smax-smin)/float(nmax))
		qvl=np.zeros(shape=(nmax,),dtype=np.float64) #define low res variables
		qsl=np.zeros(shape=(nmax,),dtype=np.float64)
		for nn in range (0,nmax): #loop over salinity classes of the low res
			qvl[nn] = np.float64(np.sum(ds/dss*qv[nn*jmax:(nn+1)*jmax])) #put the sum of high-res into the lower-res variables
			qsl[nn] = np.float64(np.sum(ds/dss*qs[nn*jmax:(nn+1)*jmax]))
		Qin = np.float64(dss*np.sum((qvl).clip(min=0),dtype=np.float64)) #calculate the bulk values
		Qout = np.float64(dss*np.sum((qvl).clip(max=0),dtype=np.float64))
		Sin = np.float64(dss*np.sum((qsl).clip(min=0),dtype=np.float64))
		Sout = np.float64(dss*np.sum((qsl).clip(max=0),dtype=np.float64))
		sin = np.float64(Sin / Qin)
		sout = np.float64(Sout / Qout)
		#Output of # of salinity classes, Knudsen bulk values
		print "%8d" % nmax,"%12.10f" % sin,"%12.10f" % sout,"%15.8f" % Qin,"%15.8f" % Qout 
		output.write(str(nmax) + '\t' + str(sin) +'\t' +str(sout) + '\t' + str(Qin) + '\t' + str(Qout) +'\n')
	return('done')

def knudsen_check_2(path_to_file,path_to_savefile):#different approach to calculate the bulk values!

	t,i,u,s,a = np.loadtxt(path_to_file, unpack=True, skiprows=1,dtype=np.float64) #load all data and split them into one variable for each column
	kmax,imax,smin,smax = np.genfromtxt(path_to_file, unpack=True, skip_footer =len(t),dtype=np.float64)# read the header and store in 4 variables
	mmax = 23
	nmax = 2**mmax
	ds=np.float64((smax-smin)/float(nmax))
	qv = np.zeros(shape=(nmax+1,),dtype=np.float64) # Initialise the variable qv with length of salinity classes with zeros;+1 so that the first class is zero and the indexing is the same as in the paper
	qs = np.zeros(shape=(nmax+1,),dtype=np.float64) # Initialise the variable qs with length of salinity classes with zeros 
	
	for tt in range (0,int(kmax)): #time loop
		for ii in range (0,int(imax)): #cross section loop
			if not ii+1 == i[ii+tt*imax]:
				sys.exit('wrong i-indexing in input file')
			if not tt+1 == t[ii+tt*imax]:
				sys.exit('wrong k-indexing in input file')
			if s[tt*imax+ii] > smax:
				sys.exit('s is greater than smax')
			if s[tt*imax+ii] < smin:
				sys.exit('s is smaller than smin')
			idx=int((s[tt*imax+ii]-smin)/(smax-smin)*(nmax))+1 #computation of the index
			if s[tt*imax+ii] == smax:
				idx = nmax
			qv[idx]=qv[idx]+u[tt*imax+ii]*a[tt*imax+ii]/kmax/ds # Adding to high-res volume flux
			qs[idx]=qs[idx]+u[tt*imax+ii]*s[tt*imax+ii]*a[tt*imax+ii]/kmax/ds # Adding high-res salinity flux
	
	"""
	      ! At this stage, qv and qs contain volume and salinity fluxes 
	      ! at very high resolution with nmax=2**mmax salinity classes
	      ! In the next step we apply coarse-graining to check the
	      ! dependence of Qin, Qout, sin and sout 
	      ! on the number of salinity classes.
	"""
	qv = qv[1:] #slice the first class away since it was only needed for correct indexing
	qs = qs[1:]
	print('#classes    inflow sal.    outflow sal.    inflow vol.    outflow vol.')
	output=open(path_to_savefile,'w')
	output.write('#classes \t inflow sal \t outflow sal. \t inflow vol. \t outflow vol. \n')
	for mm in range (1,mmax+1):#number of experiments mm
		nmax=2**mm
		jmax=2**mmax/nmax # number of low-res per high-res salinity classes
		dss = np.float64((smax-smin)/float(nmax))
		qvl=np.zeros(shape=(nmax,),dtype=np.float64) #define low res variables
		qsl=np.zeros(shape=(nmax,),dtype=np.float64)
		for nn in range (0,nmax): #loop over salinity classes of the low res
			qvl[nn] = np.float64(np.sum(ds/dss*qv[nn*jmax:(nn+1)*jmax])) #put the sum of high-res into the lower-res variables
			qsl[nn] = np.float64(np.sum(ds/dss*qs[nn*jmax:(nn+1)*jmax]))
		Qvl=np.zeros(shape=(nmax,),dtype=np.float64) #define low res Q(s) variables
		Qsl=np.zeros(shape=(nmax,),dtype=np.float64) #Q_s(s)
		for i in range(0,int(len(Qvl))): #calculate Q(s) and Q_s(s)
			Qvl[i]=np.sum(qvl[i:]*dss)
			Qsl[i]=np.sum(qsl[i:]*dss)
		#Here we take advantage that the maximum of Q corresponds to Qin 		
		#do a check for inverse estuary:
		if np.abs(np.max(Qvl))> np.abs(np.min(Qvl)):#classic estuary
			Qin = np.max(Qvl)
			Qout = Qvl[0] - Qin #if no river discharge
			Sin = np.max(Qsl)
			Sout = Qsl[0]-Sin
			sin = np.float64(Sin / Qin)
			sout = np.float64(Sout / Qout)
		elif np.abs(np.max(Qvl)) < np.abs(np.min(Qvl)):#inverse
			Qout = np.min(Qvl)
			Qin = Qvl[0]-Qout
			Sout = np.min(Qsl)
			Sin = Qsl[0]-Sout
			sin = np.float64(Sin / Qin)
			sout = np.float64(Sout / Qout)
		else:
			sys.exit('min(Q)=max(Q): strange estuary')
		#Output of # of salinity classes, Knudsen bulk values
		print "%8d" % nmax,"%12.10f" % sin,"%12.10f" % sout,"%15.8f" % Qin,"%15.8f" % Qout 
		output.write(str(nmax) + '\t' + str(sin) +'\t' +str(sout) + '\t' + str(Qin) + '\t' + str(Qout) +'\n')
	return('done')

