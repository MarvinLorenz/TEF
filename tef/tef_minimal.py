# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:53:23 2018

@author: lorenz
"""

"""
TEF function which calculates q and q_s for all timesteps and  
"""

import numpy as np
import sys


def find_extrema(x,comp):
	
	indices = []
	minmax = []
	i = 0
	while i < np.shape(x)[0]:
		if i-comp < 0:
			a = 0
		else:
			a=i-comp
		if i+comp>=len(x):
			b=None
		else:
			b=i+comp+1
		if x[i] == np.max(x[a:b]):
			indices.append(i)
			minmax.append('max')
		elif x[i] == np.min(x[a:b]):
			indices.append(i)
			minmax.append('min')
		i+=1

	return indices,minmax

def TEF_divsal(salt,vel,h,dx,N,salinity_array='auto'):
	
		#do a check if all variables have the same dimensions
	if np.shape(salt)==np.shape(vel) and np.shape(salt) == np.shape(h):
		print('input is good')
	else:
		sys.exit('salt, vel and h dont have the same dimension!')
		
	#check if dx is int or 1d array
	if type(dx) is int or type(dx) is float:
		dx = np.full((np.shape(salt)[-1],),dx,dtype=np.float64)
	elif type(dx) is np.ndarray:
		print('dx is already an array')
	elif type(dx) is np.ma.core.MaskedArray:
		print('dx is already an array')
	else:
		print(type(dx))
		sys.exit('dx format not known')
	if salinity_array == 'auto':
		#calculate number of salinity classes and create an array accordingly
		s_min = np.float64(np.min(salt))
		s_max = np.float64(np.max(salt))
		#print(s_min,s_max)
		Nmax = 2**N
		DeltaS = np.float64((s_max-s_min)/Nmax)
		s = np.arange(s_min,s_max,DeltaS,dtype=np.float64)	
	else:
		s_min=salinity_array[0]
		s_max=salinity_array[-1]
		Nmax = len(salinity_array)
		DeltaS = np.float64((s_max-s_min)/Nmax)
		s = np.arange(s_min,s_max,DeltaS,dtype=np.float64)	
#	if np.max(salt)>s_max:
#		sys.exit('a salinity greater than s_max is found')
#	if np.min(salt)<s_min:
#		sys.exit('a salinity smaller than s_min is found')	
	
	#define qv and qs:
	qv = np.zeros(shape=(len(s)+1,),dtype=np.float64) 
	qs = np.zeros(shape=(len(s)+1,),dtype=np.float64)	
	
	#check if there is a time dim and do the calculation accordingly:
	if len(np.shape(salt)) == 3:
		tmax = np.float64(np.shape(salt)[0])
		kmax = np.float64(np.shape(salt)[1])
		imax = np.float64(np.shape(salt)[2])
		tt =0
		while tt < int(tmax):
			kk=0
			while kk < int(kmax):
				ii=0
				while ii < int(imax):	
					if not np.ma.is_masked(salt[tt,kk,ii]) and salt[tt,kk,ii]>=s_min:
						idx = int((salt[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax)+1 #find the suiting salinity class
						if salt[tt,kk,ii] == s_max:
							idx = Nmax
						qv[idx] = qv[idx] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]/DeltaS/tmax
						qs[idx] = qs[idx] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]*salt[tt,kk,ii]/DeltaS/tmax
					ii+=1
				kk+=1
			tt+=1
	else:
		tmax=1
		kmax = np.shape(salt)[0]
		imax = np.shape(salt)[1]
		kk=0
		while kk < kmax:
			ii=0
			while ii < imax:	
				if not np.ma.is_masked(salt[kk,ii]) and salt[kk,ii]>=s_min:				
					idx = int((salt[kk,ii]-s_min)/(s_max-s_min)*Nmax)+1
					if salt[kk,ii] == s_max:
							idx = Nmax
					qv[idx] = qv[idx] + vel[kk,ii]*h[kk,ii]*dx[ii]/DeltaS/tmax
					qs[idx] = qs[idx] + vel[kk,ii]*h[kk,ii]*dx[ii]*salt[kk,ii]/DeltaS/tmax
				ii+=1
			kk+=1
	qv = qv[1:]
	qs = qs[1:]
	#Calculation of Qv and Qs:
	Qv = np.zeros(shape=np.shape(qv))
	Qs = np.zeros(shape=np.shape(qs))
	for i in range(0,int(len(qv))): #calculate Q(s) and Q_s(s)
		Qv[i]=np.sum(qv[i:]*DeltaS)
		Qs[i]=np.sum(qs[i:]*DeltaS)
	
	#now find the Maxima and Minima	
	if N > 6:
		teiler = 50
	elif N == 1:
		teiler = Nmax+1
	else:
		teiler = 1
	ind,minmax = find_extrema(Qv,int(Nmax/teiler))# find the extrema
	div_sal = []
	for i in range(0,len(ind)):
		div_sal.append(s_min+DeltaS*ind[i]) #compute dividing salinities
	print(Nmax,ind,minmax,div_sal)
	index=[]
	ii=0
	while ii < len(ind):
		#print(ii)
		if minmax[ii] == minmax[ii-1]:
			index.append(ii-1)
		ii+=1
	minmax = np.asarray(minmax)
	ind = np.asarray(ind)
	if len(ind)>2:
		ind = np.delete(ind, index)
		minmax = np.delete(minmax, index)
		div_sal = np.delete(div_sal, index)
	for i in range(0,len(ind)):
		print(s_min+DeltaS*ind[i])
	print(Nmax,ind,minmax,div_sal)
	ind = ind[:-1]
	#compute the bulk values
	transports = []
	saltfluxes = []
	salinities = []
	print(ind)

	for i in range(0,len(ind)-1):
		transports.append(Qv[ind[i]]-Qv[ind[i+1]])
		saltfluxes.append(Qs[ind[i]]-Qs[ind[i+1]])
		salinities.append(saltfluxes[i]/transports[i])
	transports.append(Qv[ind[-1]])
	saltfluxes.append(Qs[ind[-1]])
	salinities.append(saltfluxes[-1]/transports[-1])

	return(qv, Qv, s, div_sal, transports, salinities)

def TEF_q_timeseries(salt,vel,h,dx,salinity_array,time,dt):
	
	if np.shape(salt)==np.shape(vel) and np.shape(salt) == np.shape(h) and np.shape(salt)[0] == np.shape(time)[0]:
			print('input is good')
	else:
		sys.exit('shapes dont match')
	#q_ar=np.zeros(shape=(int(len(time)/dt),len(salinity_array)))	
	print(np.shape(salt))
	q_ar=[]	
	time_new = []
	print(int(len(time)/dt))	
	t=0
	while t<int(len(time)/dt):
		print((t+1)*dt, len(time))
		if (t+1)*dt<len(time):
			qv,qs,s,bulk_values = TEF_minimal(salt[t*dt:(t+1)*dt,:,:],vel[t*dt:(t+1)*dt,:,:],h[t*dt:(t+1)*dt,:,:],dx,10,salinity_array=salinity_array)
			q_ar.append(qv)
			time_new.append(time[t*dt])
		t+=1
	
	return(np.asarray(q_ar), np.asarray(time_new))

def TEF_entrainment(salt_1,vel_1,h_1,dx_1,dy_1,salt_2,vel_2,h_2,dx_2,dy_2,s,N,dt):#following Wang et al. 2017
	print(np.shape(salt_1),np.shape(vel_1),np.shape(h_1),np.shape(salt_2),np.shape(vel_2),np.shape(h_2))
	#need two timepoints to calculate dV/dt, but for TEF-profile only the later of the two timesteps is analyzed	
	print('doing standard TEF for left and right boundary')
	qv_1,qs_1,s_1,bulk_values_1 = TEF_minimal(salt_1[-1,:,:],vel_1[-1,:,:],h_1[-1,:,:],dx_1,N,salinity_array=s)
	print('TEF2')
	qv_2,qs_2,s_2,bulk_values_2 = TEF_minimal(salt_2[-1,:,:],vel_2[-1,:,:],h_2[-1,:,:],dx_2,N,salinity_array=s)
	print('calculating dV/ds')
	V_left = np.zeros(shape=(len(s_1)+1,2),dtype=np.float64) #shape=(salt,time) len(time)=2
	V_right = np.zeros(shape=(len(s_2)+1,2),dtype=np.float64)
	s_min=s[0]
	s_max=s[-1]
	Nmax = len(s)
	DeltaS = np.float64((s_max-s_min)/Nmax)
	
	kmax = np.shape(salt_1)[1]
	imax = np.shape(salt_1)[2]
	tt=0
	while tt <2:
		kk=0
		while kk < kmax:
			ii=0
			while ii < imax:	
				if not np.ma.is_masked(salt_1[tt,kk,ii]):	
					idx = int((salt_1[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax)+1
					if salt_1[tt,kk,ii] == s_max:
							idx = Nmax
					V_left[idx,tt] = V_left[idx,tt] + 0.5*h_1[tt,kk,ii]*dx_1[ii]*dy_1[ii]/DeltaS #dV ds, 0.5 because of T point position in C grid
				if not np.ma.is_masked(salt_2[tt,kk,ii]):
					idy = int((salt_2[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax)+1
					if salt_2[tt,kk,ii] == s_max:
							idy = Nmax
					V_right[idy,tt] = V_right[idy,tt] + 0.5*h_2[tt,kk,ii]*dx_2[ii]*dy_2[ii]/DeltaS
				ii+=1
			kk+=1
		tt+=1
	DVdt = (V_left[:,1]-V_left[:,0]+V_right[:,1]-V_right[:,0])/dt #d^2V dsdt
	DVdt = DVdt[1:]	
	
	print(np.shape(s), np.shape(DVdt))
	
	#Calculate Q_e(s)		
	#A:
	A = qv_1-qv_2+DVdt
	
	print(np.sum(A*DeltaS))
	#B:
	i=0
	B=np.zeros(shape=np.shape(s))
	while i < len(s):
		B[i]=np.sum(A[i:]*DeltaS)
		i+=1
	#C: molecular flux:
	C=np.zeros(shape=np.shape(s))
	i=0
	while i < len(s):	
		C[i]= np.sum(B[i:]*DeltaS)
		i+=1
	#D: dB/ds
	#D=(B[1:]-B[:-1])/DeltaS
	return(qv_1,qv_2,DVdt,A,B,C)
	
	


def TEF_minimal(salt,vel,h,dx,N,salinity_array='auto'):

	#do a check if all variables have the same dimensions
	if np.shape(salt)==np.shape(vel) and np.shape(salt) == np.shape(h):
		print('input is good')
	else:
		sys.exit('salt, vel and h dont have the same dimension!')
		
	#check if dx is int or 1d array
	if type(dx) is int or type(dx) is float:
		dx = np.full((np.shape(salt)[-1],),dx,dtype=np.float64)
	elif type(dx) is np.ndarray:
		print('dx is already an array')
	elif type(dx) is np.ma.core.MaskedArray:
		print('dx is already an array')
	else:
		print(type(dx))
		sys.exit('dx format not known')
	if salinity_array == 'auto':
		#calculate number of salinity classes and create an array accordingly
		s_min = np.float64(np.min(salt))
		s_max = np.float64(np.max(salt))
		#print(s_min,s_max)
		Nmax = 2**N
		DeltaS = np.float64((s_max-s_min)/Nmax)
		s = np.arange(s_min,s_max,DeltaS,dtype=np.float64)	
	else:
		s_min=salinity_array[0]
		s_max=salinity_array[-1]
		Nmax = len(salinity_array)
		DeltaS = np.float64((s_max-s_min)/Nmax)
		s = np.arange(s_min,s_max,DeltaS,dtype=np.float64)	
#	if np.max(salt)>s_max:
#		sys.exit('a salinity greater than s_max is found')
#	if np.min(salt)<s_min:
#		sys.exit('a salinity smaller than s_min is found')	
	
	#define qv and qs:
	qv = np.zeros(shape=(len(s)+1,),dtype=np.float64) 
	qs = np.zeros(shape=(len(s)+1,),dtype=np.float64)	
	
	#check if there is a time dim and do the calculation accordingly:
	if len(np.shape(salt)) == 3:
		tmax = np.float64(np.shape(salt)[0])
		kmax = np.float64(np.shape(salt)[1])
		imax = np.float64(np.shape(salt)[2])
		tt =0
		while tt < int(tmax):
			kk=0
			while kk < int(kmax):
				ii=0
				while ii < int(imax):	
					if not np.ma.is_masked(salt[tt,kk,ii]) and salt[tt,kk,ii]>=s_min:
						idx = int((salt[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax)+1 #find the suiting salinity class
						if salt[tt,kk,ii] == s_max:
							idx = Nmax
						qv[idx] = qv[idx] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]/DeltaS/tmax
						qs[idx] = qs[idx] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]*salt[tt,kk,ii]/DeltaS/tmax
					ii+=1
				kk+=1
			tt+=1
	else:
		tmax=1
		kmax = np.shape(salt)[0]
		imax = np.shape(salt)[1]
		kk=0
		while kk < kmax:
			ii=0
			while ii < imax:	
				if not np.ma.is_masked(salt[kk,ii]) and salt[kk,ii]>=s_min:				
					idx = int((salt[kk,ii]-s_min)/(s_max-s_min)*Nmax)+1
					if salt[kk,ii] == s_max:
							idx = Nmax
					qv[idx] = qv[idx] + vel[kk,ii]*h[kk,ii]*dx[ii]/DeltaS/tmax
					qs[idx] = qs[idx] + vel[kk,ii]*h[kk,ii]*dx[ii]*salt[kk,ii]/DeltaS/tmax
				ii+=1
			kk+=1
	qv = qv[1:]
	qs = qs[1:]
	
	Qin = np.float64(np.sum((DeltaS*qv).clip(min=0),dtype=np.float64)) #calculate the bulk values
	Qout = np.float64(np.sum((DeltaS*qv).clip(max=0),dtype=np.float64))
	Sin = np.float64(np.sum((DeltaS*qs).clip(min=0),dtype=np.float64))
	Sout = np.float64(np.sum((DeltaS*qs).clip(max=0),dtype=np.float64))
	sin = np.float64(Sin / Qin)
	sout = np.float64(Sout / Qout)

	bulk_values = [Qin, Qout, Sin, Sout, sin, sout]
	
	return(qv, qs, s, bulk_values)
	
def TEF_minimal_2d(salt,temp,vel,h,dx,N):
	
	#do a check if all variables have the same dimensions
	if np.shape(salt)==np.shape(vel) and np.shape(salt) == np.shape(h) and np.shape(salt)==np.shape(temp):
		print('input is good')
	else:
		sys.exit('salt, temp, vel and h dont have the same dimension!')
		
	#check if ds is int or 1d array
	if type(dx) is int or type(dx) is float:
		dx = np.full((np.shape(salt)[-1],),dx,dtype=np.float64)
	elif type(dx) is np.ndarray:
		print('dx is already an array')
	else:
		print(type(dx))
		sys.exit('dx format not known')
	
	#calculate number of salinity classes and create an array accordingly
	s_min = np.float64(np.min(salt))
	s_max = np.float64(np.max(salt))
	#print(s_min,s_max)
	t_min = np.float64(np.min(temp))
	t_max = np.float64(np.max(temp))
	#print(t_min,t_max)
	Nmax = 2**N
	DeltaS = np.float64((s_max-s_min)/Nmax)
	s = np.arange(s_min,s_max,DeltaS,dtype=np.float64)	
	DeltaT = np.float64((t_max-t_min)/Nmax)
	t = np.arange(t_min,t_max,DeltaT,dtype=np.float64)	
	
	
	#define qv and qs:
	qv = np.zeros(shape=(len(s)+1,len(t)+1),dtype=np.float64) 
	qs = np.zeros(shape=(len(s)+1,len(t)+1),dtype=np.float64)	
	qt = np.zeros(shape=(len(s)+1,len(t)+1),dtype=np.float64)
	
	#check if there is a time dim and do the calculation accordingly:
	#print(len(np.shape(salt)))	
	if len(np.shape(salt)) == 3:
		tmax = np.float64(np.shape(salt)[0])
		kmax = np.float64(np.shape(salt)[1])
		imax = np.float64(np.shape(salt)[2])
		tt =0
		while tt < int(tmax):
			kk=0
			while kk < int(kmax):
				ii=0
				while ii < int(imax):						
					idx = int((salt[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax)+1 #find the suiting salinity class
					#print((salt[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax, 'salt')
					if salt[tt,kk,ii] == s_max:
						#print('found smax')
						idx = Nmax
					idy = int((temp[tt,kk,ii]-t_min)/(t_max-t_min)*Nmax)+1 #find the suiting salinity class
					#print((temp[tt,kk,ii]-t_min)/(t_max-t_min)*Nmax, 'temp')					
					if temp[tt,kk,ii] == t_max:
						#print('found tmax')
						idy = Nmax
#					qv[idx,idy] = qv[idx,idy] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]/DeltaS/tmax
#					qs[idx,idy] = qs[idx,idy] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]*salt[tt,kk,ii]/DeltaS/tmax
					qv[idx,idy] = qv[idx,idy] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]/DeltaS/tmax/DeltaT
					qs[idx,idy] = qs[idx,idy] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]*salt[tt,kk,ii]/DeltaS/tmax/DeltaT
					qt[idx,idy] = qt[idx,idy] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]*temp[tt,kk,ii]/DeltaS/tmax/DeltaT
					ii+=1
				kk+=1
			tt+=1
	else:
		tmax=1
		kmax = np.shape(salt)[0]
		imax = np.shape(salt)[1]
		kk=0
		while kk < kmax:
			ii=0
			while ii < imax:						
				idx = int((salt[kk,ii]-s_min)/(s_max-s_min)*Nmax)+1
				if salt[kk,ii] == s_max:
						idx = Nmax
				idy = int((temp[kk,ii]-t_min)/(t_max-t_min)*Nmax)+1
				if temp[kk,ii] == t_max:
						idy = Nmax
				qv[idx,idy] = qv[idx,idy] + vel[kk,ii]*h[kk,ii]*dx[ii]/DeltaS/tmax/DeltaT
				qs[idx,idy] = qs[idx,idy] + vel[kk,ii]*h[kk,ii]*dx[ii]*salt[kk,ii]/DeltaS/tmax/DeltaT
				qt[idx,idy] = qs[idx,idy] + vel[kk,ii]*h[kk,ii]*dx[ii]*temp[kk,ii]/DeltaS/tmax/DeltaT
				ii+=1
			kk+=1
	#print(qv[,:])
	qv = qv[1:,1:]
	qs = qs[1:,1:]
	qt = qt[1:,1:]
	
	Qin = np.float64(np.sum((DeltaS*DeltaT*qv).clip(min=0),dtype=np.float64)) #calculate the bulk values
	Qout = np.float64(np.sum((DeltaS*DeltaT*qv).clip(max=0),dtype=np.float64))
	Sin = np.float64(np.sum((DeltaS*DeltaT*qs).clip(min=0),dtype=np.float64))
	Sout = np.float64(np.sum((DeltaS*DeltaT*qs).clip(max=0),dtype=np.float64))
	Tin = np.float64(np.sum((DeltaS*DeltaT*qt).clip(min=0),dtype=np.float64))
	Tout = np.float64(np.sum((DeltaS*DeltaT*qt).clip(max=0),dtype=np.float64))
	sin = np.float64(Sin / Qin)
	sout = np.float64(Sout / Qout)
	tin = np.float64(Tin / Qin)
	tout = np.float64(Tout / Qout)

	bulk_values = [Qin, Qout, Sin, Sout, sin, sout, Tin, Tout, tin, tout]
	
	#values for salt TEF by integrating over all temperatures
	qvl = np.sum(qv, axis=1)
	qsl = np.sum(qs, axis=1)
	qtl = np.sum(qt, axis=1)
	Qin = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(min=0),dtype=np.float64)) #calculate the bulk values
	Qout = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(max=0),dtype=np.float64))
	Sin = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(min=0),dtype=np.float64))
	Sout = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(max=0),dtype=np.float64))
	Tin = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(min=0),dtype=np.float64))
	Tout = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(max=0),dtype=np.float64))
	sin = np.float64(Sin / Qin)
	sout = np.float64(Sout / Qout)
	tin = np.float64(Tin / Qin)
	tout = np.float64(Tout / Qout)
	
	bulk_val_sal = [Qin, Qout, Sin, Sout, sin, sout, Tin, Tout, tin, tout]
	
	#values for temp TEF by integrating over all salinities
	qvl = np.sum(qv, axis=0)
	qsl = np.sum(qs, axis=0)
	qtl = np.sum(qt, axis=0)
	Qin = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(min=0),dtype=np.float64)) #calculate the bulk values
	Qout = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(max=0),dtype=np.float64))
	Sin = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(min=0),dtype=np.float64))
	Sout = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(max=0),dtype=np.float64))
	Tin = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(min=0),dtype=np.float64))
	Tout = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(max=0),dtype=np.float64))
	sin = np.float64(Sin / Qin)
	sout = np.float64(Sout / Qout)
	tin = np.float64(Tin / Qin)
	tout = np.float64(Tout / Qout)
	
	bulk_val_temp = [Qin, Qout, Sin, Sout, sin, sout, Tin, Tout, tin, tout]
	
	return(qv, qs, qt, s, t, bulk_val_sal)
	
def TEF_minimal_2d_k(salt,temp,vel,A,N):
	
	#do a check if all variables have the same dimensions
	if np.shape(salt)==np.shape(vel) and np.shape(salt)==np.shape(temp):
		print('input is good')
	else:
		sys.exit('salt, temp and vel dont have the same dimension!')
	
	
	#calculate number of salinity classes and create an array accordingly
	s_min = np.float64(np.min(salt))
	s_max = np.float64(np.max(salt))
	#print(s_min,s_max)
	t_min = np.float64(np.min(temp))
	t_max = np.float64(np.max(temp))
	#print(t_min,t_max)
	Nmax = 2**N
	DeltaS = np.float64((s_max-s_min)/Nmax)
	s = np.arange(s_min,s_max,DeltaS,dtype=np.float64)	
	DeltaT = np.float64((t_max-t_min)/Nmax)
	t = np.arange(t_min,t_max,DeltaT,dtype=np.float64)	
	
	
	#define qv and qs:
	qv = np.zeros(shape=(len(s)+1,len(t)+1),dtype=np.float64) 
	qs = np.zeros(shape=(len(s)+1,len(t)+1),dtype=np.float64)	
	qt = np.zeros(shape=(len(s)+1,len(t)+1),dtype=np.float64)
	
	#check if there is a time dim and do the calculation accordingly:
	#print(len(np.shape(salt)))	
	if len(np.shape(salt)) == 3:
		tmax = np.float64(np.shape(salt)[0])
		kmax = np.float64(np.shape(salt)[1])
		imax = np.float64(np.shape(salt)[2])
		tt =0
		while tt < int(tmax):
			kk=0
			while kk < int(kmax):
				ii=0
				while ii < int(imax):						
					if np.ma.is_masked(salt[tt,kk,ii]):
						idx=0
						idy=0
					else:
						idx = int((salt[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax)+1 #find the suiting salinity class
					#print((salt[tt,kk,ii]-s_min)/(s_max-s_min)*Nmax, 'salt')
						if salt[tt,kk,ii] == s_max:
						#print('found smax')
							idx = Nmax
						idy = int((temp[tt,kk,ii]-t_min)/(t_max-t_min)*Nmax)+1 #find the suiting salinity class
					#print((temp[tt,kk,ii]-t_min)/(t_max-t_min)*Nmax, 'temp')					
						if temp[tt,kk,ii] == t_max:
						#print('found tmax')
							idy = Nmax
#					qv[idx,idy] = qv[idx,idy] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]/DeltaS/tmax
#					qs[idx,idy] = qs[idx,idy] + vel[tt,kk,ii]*h[tt,kk,ii]*dx[ii]*salt[tt,kk,ii]/DeltaS/tmax
					qv[idx,idy] = qv[idx,idy] + vel[tt,kk,ii]*A[kk,ii]/DeltaS/tmax/DeltaT
					qs[idx,idy] = qs[idx,idy] + vel[tt,kk,ii]*A[kk,ii]*salt[tt,kk,ii]/DeltaS/tmax/DeltaT
					qt[idx,idy] = qt[idx,idy] + vel[tt,kk,ii]*A[kk,ii]*temp[tt,kk,ii]/DeltaS/tmax/DeltaT
					ii+=1
				kk+=1
			tt+=1
	else:
		tmax=1
		kmax = np.shape(salt)[0]
		imax = np.shape(salt)[1]
		kk=0
		while kk < kmax:
			ii=0
			while ii < imax:		
				if np.ma.is_masked(salt[kk,ii]):
					idx=0
					idy=0
				else:				
					idx = int((salt[kk,ii]-s_min)/(s_max-s_min)*Nmax)+1
					if salt[kk,ii] == s_max:
							idx = Nmax
					idy = int((temp[kk,ii]-t_min)/(t_max-t_min)*Nmax)+1
					if temp[kk,ii] == t_max:
							idy = Nmax
				qv[idx,idy] = qv[idx,idy] + vel[kk,ii]*A[kk,ii]/DeltaS/tmax/DeltaT
				qs[idx,idy] = qs[idx,idy] + vel[kk,ii]*A[kk,ii]*salt[kk,ii]/DeltaS/tmax/DeltaT
				qt[idx,idy] = qs[idx,idy] + vel[kk,ii]*A[kk,ii]*temp[kk,ii]/DeltaS/tmax/DeltaT
				ii+=1
			kk+=1
	#print(qv[,:])
	qv = qv[1:,1:]
	qs = qs[1:,1:]
	qt = qt[1:,1:]
	
	Qin = np.float64(np.sum((DeltaS*DeltaT*qv).clip(min=0),dtype=np.float64)) #calculate the bulk values
	Qout = np.float64(np.sum((DeltaS*DeltaT*qv).clip(max=0),dtype=np.float64))
	Sin = np.float64(np.sum((DeltaS*DeltaT*qs).clip(min=0),dtype=np.float64))
	Sout = np.float64(np.sum((DeltaS*DeltaT*qs).clip(max=0),dtype=np.float64))
	Tin = np.float64(np.sum((DeltaS*DeltaT*qt).clip(min=0),dtype=np.float64))
	Tout = np.float64(np.sum((DeltaS*DeltaT*qt).clip(max=0),dtype=np.float64))
	sin = np.float64(Sin / Qin)
	sout = np.float64(Sout / Qout)
	tin = np.float64(Tin / Qin)
	tout = np.float64(Tout / Qout)

	bulk_values = [Qin, Qout, Sin, Sout, sin, sout, Tin, Tout, tin, tout]
	
	#values for salt TEF by integrating over all temperatures
	qvl = np.sum(qv, axis=1)
	qsl = np.sum(qs, axis=1)
	qtl = np.sum(qt, axis=1)
	Qin = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(min=0),dtype=np.float64)) #calculate the bulk values
	Qout = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(max=0),dtype=np.float64))
	Sin = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(min=0),dtype=np.float64))
	Sout = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(max=0),dtype=np.float64))
	Tin = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(min=0),dtype=np.float64))
	Tout = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(max=0),dtype=np.float64))
	sin = np.float64(Sin / Qin)
	sout = np.float64(Sout / Qout)
	tin = np.float64(Tin / Qin)
	tout = np.float64(Tout / Qout)
	
	bulk_val_sal = [Qin, Qout, Sin, Sout, sin, sout, Tin, Tout, tin, tout]
	
	#values for temp TEF by integrating over all salinities
	qvl = np.sum(qv, axis=0)
	qsl = np.sum(qs, axis=0)
	qtl = np.sum(qt, axis=0)
	Qin = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(min=0),dtype=np.float64)) #calculate the bulk values
	Qout = np.float64(np.sum((DeltaS*DeltaT*qvl).clip(max=0),dtype=np.float64))
	Sin = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(min=0),dtype=np.float64))
	Sout = np.float64(np.sum((DeltaS*DeltaT*qsl).clip(max=0),dtype=np.float64))
	Tin = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(min=0),dtype=np.float64))
	Tout = np.float64(np.sum((DeltaS*DeltaT*qtl).clip(max=0),dtype=np.float64))
	sin = np.float64(Sin / Qin)
	sout = np.float64(Sout / Qout)
	tin = np.float64(Tin / Qin)
	tout = np.float64(Tout / Qout)
	
	bulk_val_temp = [Qin, Qout, Sin, Sout, sin, sout, Tin, Tout, tin, tout]
	
	return(qv, qs, qt, s, t, bulk_values, bulk_val_sal, bulk_val_temp)