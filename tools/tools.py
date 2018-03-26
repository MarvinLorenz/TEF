# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:44:05 2018

@author: lorenz
"""

import numpy as np
import sys
sys.path.append('/fast/lorenz/Programs/seawater')
import gsw as gsw
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import netCDF4

#module to do a rotation of an arbitrary transect contructed with nctransect, uu and vv need to have the dimensions
#time,k,distance, where k is the vertical coordinate and distance the horizontal one lonc,latc and distance are only 
#dependent on distance. time is only optional, k and distance are necessary!
def tef_q_timeseries_show(q_ar,s,time,path,filename,time_since=0):
	
	print(np.shape(q_ar), np.shape(time), np.shape(s))
	if np.shape(q_ar)[0]==np.shape(time)[0] and np.shape(q_ar)[1] == np.shape(s)[0]:
			print('input is good')
	else:
		sys.exit('shapes dont match')
	#convert imte axis im time_since is given
	if time_since !=0:
		time_date = netCDF4.num2date(time,time_since)
		x_label=[time_date[0].strftime('%b')+' '+str(time_date[0].year)]
		ind=[time[0]]
		i=1
		while i < len(time):
			if time_date[i].month != time_date[i-1].month:
				x_label.append(time_date[i].strftime('%b')+' '+str(time_date[i].year))
				ind.append(time[i])
			i+=1
	
	q_ar = q_ar.transpose()
	#q_ar = np.ma.masked_where(q_ar==0, q_ar)
	fs=20
	plt.figure(figsize=(15,5))
	ax1 = plt.subplot(111) #left plot of q
	plt.gca().invert_yaxis()
	pc1=ax1.pcolor(time,s,q_ar,cmap = 'bwr', vmin=-np.max(q_ar)/10, vmax = np.max(q_ar)/10)
	cbar1 = plt.colorbar(pc1, format='%d')
	cbar1.set_label('$-\partial Q(s) / \partial s $ [m$^3$s$^{-1}$ / (g/kg)]', fontsize = fs)
	ax1.set_ylabel('salinity [g/kg]', fontsize = fs)
	ax1.set_xlabel('time', fontsize = fs)
	ax1.tick_params('both', colors='black', labelsize=fs)
	plt.ylim([np.max(s),np.min(s)])
	plt.xlim([np.min(time),np.max(time)])
	if time_since!=0:
		plt.xticks(ind, x_label, rotation = 45, fontsize = fs)
	plt.gcf().subplots_adjust(bottom=0.15)
	plt.grid()
	#fig.subplots_adjust(hspace=0)
	print('saving png...')
	#plt.savefig(directory_TEF+ 'Q_and_q_example'+save_name+'.pdf',format = 'pdf',bbox_inches='tight')
	plt.savefig(path+filename,format = 'png',bbox_inches='tight')
	#plt.show
	plt.close()
	return('done')
	
def tef_2d_show(p,s,t,path,filename,plot_dens=True):
	
	fs=20
	p=p.transpose() #transpose q from tef_2d output, to be able to plot it
	if plot_dens:
		#calculate density lines in plot
		ydim=np.shape(t)[0]
		xdim=np.shape(s)[0]
	
		dens = np.zeros(shape=(xdim,ydim))
		for i in range(0,int(xdim)):
			#print(i)	
			dens[:,i]=gsw.rho(s[i],t,0)
		dens = dens -1000.0
	vmax = np.max(p)
	#do the plot:
	fig = plt.figure(figsize=(10,10))
	fig.set_size_inches(10,10)
	ax = plt.subplot()
	p4 = ax.pcolor(s,t, p, cmap='seismic',vmin=-vmax,vmax=vmax)
	#p4=ax.scatter(ss, tt, c=qq/1000000.0, marker='o', s = 4**2,cmap='seismic',vmin=-vmax,vmax=vmax,edgecolors='face')
	#levels = np.linspace(-vmax,vmax,64)	
	#p4=ax.contourf(s_new2, t_new2, qv/1000000.0, levels=levels, cmap='seismic',vmin=-vmax,vmax=vmax)
	CS = ax.contour(s,t,dens, linestyles='dashed', colors='k')
	plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
	ax.tick_params('both', colors='black', labelsize=fs)
	ax.set_xlabel('salinity [g/kg]', fontsize = fs)
	ax.set_ylabel('temperature [$^\circ$C]', fontsize = fs)
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	v = np.linspace(-vmax, vmax, 10, endpoint=True)
	cbar=plt.colorbar(p4,ticks = v, format = "%.2f")	
	cbar.set_label('p [m$^3$s$^{-1}$(g/kg)$^{-1}$K$^{-1}$]', fontsize = fs)
	cbar.ax.tick_params(labelsize=fs)
	plt.ylim([t[0],t[-1]])
	plt.xlim([s[0],s[-1]])
	plt.gcf().subplots_adjust(bottom=0.15)
	print('saving png...')
	plt.savefig(path+filename,format = 'png',bbox_inches='tight')
	plt.close()	
	
	return('done')

def tef_show(q,Q,s,path,filename,hlines=[]):
	
	fs=20
	plt.figure(figsize=(15,10))
	ax1 = plt.subplot(121) #left plot of q
	plt.gca().invert_yaxis()
	ax1.plot(q, s, color = "black")
	#p2 = ax.plot(Qv/1000000.0, s_new2, color = "blue")
	ax1.axvline(0,ls='-',color = 'black')
	if len(hlines)!=0:
		for i in range(0,len(hlines)):
			ax1.axhline(hlines[i],ls='-',color = 'black')

	ax1.tick_params('both', colors='black', labelsize=fs)
	#ax.set_xlabel('$Q \: $[Sv]', fontsize = fs)
	ax1.set_ylabel('salinity [g/kg]', fontsize = fs)
	xfmt=ScalarFormatter(useMathText=True)
	xfmt.set_powerlimits((-4,7))
	ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax1.ticklabel_format(style='sci', scilimits=(1,4), axis='x')	
	ax1.xaxis.major.formatter._useMathText = True
	#ax1.xaxis.set_major_formatter(xfmt)
	#ax1.xaxis.set_major_formatter(xfmt.set_powerlimits((-4,7)))
	#plt.ticklabel_format(style='sci',axis='x')
	ax1.set_xlabel('$-\partial Q(s) / \partial s $ [m$^3$s$^{-1}$ / (g/kg)]', fontsize = fs)
#	plt.xlim([np.min(q)*1.1,np.max(q)*1.1])
	#plt.ylim([35,29])
	
	ax2 = plt.subplot(122) #right plot of Q
	plt.gca().invert_yaxis()
	ax2.plot(Q, s, color = "black")
	#p2 = ax.plot(Qv/1000000.0, s_new2, color = "blue")
	ax2.axvline(0,ls='-',color = 'black')
	if len(hlines)!=0:
		for i in range(0,len(hlines)):
			ax2.axhline(hlines[i],ls='-',color = 'black')

	ax2.tick_params('both', colors='black', labelsize=fs)
	#ax.set_xlabel('$Q \: $[Sv]', fontsize = fs)
	ax2.set_ylabel('salinity [g/kg]', fontsize = fs)
	ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax2.ticklabel_format(style='sci', scilimits=(1,4), axis='x')	
	ax2.xaxis.major.formatter._useMathText = True
	ax2.set_xlabel('$Q(s)$ [m$^3$s$^{-1}$]', fontsize = fs)
#	plt.xlim([np.min(Q)*1.1,np.max(Q)*1.1])
#	plt.ylim([np.max(s),np.min(s)])
	#plt.ylim([35,29])
	plt.gcf().subplots_adjust(bottom=0.15)
	#fig.subplots_adjust(hspace=0)
	print('saving png...')
	#plt.savefig(directory_TEF+ 'Q_and_q_example'+save_name+'.pdf',format = 'pdf',bbox_inches='tight')
	plt.savefig(path+filename,format = 'png',bbox_inches='tight')
	#plt.show
	plt.close()
	
	return('done')

def running_mean(data_1d, step): #data should have one discrete timestep, will give back the mean of values of step/2 before and after initial values
	data_return = np.zeros(shape=(np.shape(data_1d)))
	i=0
	while i < len(data_1d):
		if i < step/2.0:
			data_return[i] = np.mean(data_1d[:step])
		elif len(data_1d)-i < step/2.0:
			data_return[i] = np.mean(data_1d[len(data_1d)-step:])	
		else:
			data_return[i] = np.mean(data_1d[i-step/2:i+step/2])
		i+=1
	return(data_return)

def interpolate_to_z_1d(z_levels, bathymetry, h, var_to_interpolate):

	z_old = bathymetry-np.cumsum(h)
	func = interpolate.interp1d(z_old, var_to_interpolate)
	var_out=func(z_levels)	
	#print(var_out)
	
	return (var_out)
	
def interpolate_to_z_2d(z_levels, bathymetry, h, var_to_interpolate):
	
	var_out=np.zeros(shape=(len(z_levels),np.shape(h)[1]))
	print(np.shape(bathymetry))
	print(np.shape(h))
	z_old = bathymetry-np.cumsum(h,axis=0)
	print(np.shape(z_old))
	i=0
	while i < np.shape(var_out)[1]:
		var_out[:,i] = np.interp(-z_levels,-z_old[:,i], var_to_interpolate[:,i])
		i+=1
	#var_out=func(z_levels)	
	#print(var_out)
	
	return (var_out)
	
def interpolate_to_z_4d(z_levels, bathymetry, h, var_to_interpolate):
	
	z_levels=np.asarray(z_levels)
	z_old= np.zeros(shape=(np.shape(h)[0],np.shape(h)[1],np.shape(h)[2],np.shape(h)[3]))		
	t=0
	while t < np.shape(z_old)[0]:
		#print(t, 'of', np.shape(z_old)[0]+1)
		z_old[t,:,:,:] = bathymetry-np.cumsum(h[t,:,:,:],axis=0)+0.5*h[t,:,:,:] #to be in the middle of the layer
		t+=1
	print('z_old done')
	#print(z_old[0,:,15,15])
	#print(var_to_interpolate[0,:,15,15])
	var_out=	np.zeros(shape=(np.shape(h)[0],len(z_levels),np.shape(h)[2],np.shape(h)[3]))	
	t=0
	while t < np.shape(var_out)[0]:
		print(t, 'of', np.shape(var_out)[0])
		j=0
		while j <  np.shape(var_out)[2]:
			i=0
			while i < np.shape(var_out)[3]:
#				print(bathymetry[j,i])
#				print(z_levels)
#				print(z_old[t,:,j,i])
				var_out[t,:,j,i] = np.interp(-z_levels, -z_old[t,:,j,i],var_to_interpolate[t,:,j,i])
				k=0
				while k < np.shape(var_out)[1]:
					if not bathymetry[j,i] > z_levels[k]:
						var_out[t,k,j,i] = -9999.0
					k+=1
				i+=1
			j+=1
		t+=1
	#print(np.interp([0.5,1.5],[0,1,2,3,4], [4,3,2,1,0]))
	#print(np.interp(z_levels, z_old[0,:,15,15],var_to_interpolate[0,:,15,15]))
	var_out = np.ma.masked_where(var_out == -9999.0, var_out)
	#print(z_levels)
	#print(var_out[0,:,15,15])
	return (var_out, z_old)
	
def convert_to_scatter(s,t,q):
		
	#print(np.shape(q))		
	snew = np.tile(s,np.shape(q)[0])
	tnew = np.repeat(t,np.shape(q)[1])
	#print(snew)
	qnew = q.flatten()
	#print('converted to scatter')
	return(snew,tnew,qnew)

def rotate(uu,vv,latc,lonc,distance):

	if np.shape(uu) == ():
		sys.exit('cant rotate single values!')
	#do a dimension check
	if np.shape(uu) == np.shape(vv):
		print('uu and vv have the same dimensions')
		if len(latc) == len(lonc) and len(latc) == len(distance):
			if np.shape(uu)[-1]==len(latc):
				print('all dimensions are compatible!')
			else:
				print('dimesions of uu/vv and latc/lonc/distance dont match')
		else:
			sys.exit('dimensions of latc, lonc and distance dont match')
	else:
		sys.exit('dimensions of uu and vv dont match')
	
	#check im time:
	if len(np.shape(uu)) == 3:
	
		dlon = lonc[1:]-lonc[:-1] #compute the lon/lat differences
		dlat = latc[1:]-latc[:-1]
	
		
		u_parallel = np.zeros(shape=(np.shape(uu)[0],np.shape(uu)[1],np.shape(uu)[2]-1)) #define the new variables
		u_perpendicular = np.zeros(shape=np.shape(u_parallel))
		alpha = np.zeros(shape=np.shape(u_parallel)) #angle, the velocities have to be rotated
	
		latlonbetrag = np.sqrt(dlon*dlon+dlat*dlat) #
		#print(np.shape(latlonbetrag))
	
		uu = uu[:,:,:-1] #we lose the last value since we need to compute a difference
		vv = vv[:,:,:-1]
		distance = distance[1:]-distance[:-1]
		
		tt = 0
		while tt < np.shape(uu)[0]:
			uvbetrag=np.sqrt(uu[tt,:,:]*uu[tt,:,:]+vv[tt,:,:]*vv[tt,:,:])
			i=0 
			while i < np.shape(u_parallel)[1]:
				j=0
				while j < np.shape(u_parallel)[2]:
					if vv[tt,i,j] == 0:
						vv[tt,i,j] = 1.e-10
					if uu[tt,i,j] == 0:
						uu[tt,i,j] = 1.e-10
					beta = np.arccos((uu[tt,i,j]*dlon[j]+vv[tt,i,j]*dlat[j])/(latlonbetrag[j]*uvbetrag[i][j]))
					if np.sign(dlat[j]) == np.sign(vv[tt,i,j]): #check if r and u are both above or below 0
						#check if 1,2 or 3,4
						if np.sign(dlat[j]) >=0: #1,2 quadrant
							if np.sign(dlon[j]) != np.sign(uu[tt,i,j]): # check if r und u are in same quadrant
								if dlon[j] >= 0:
									alpha[tt,i,j] = beta
								else:
									alpha[tt,i,j] = -beta
							if np.sign(dlon[j]) == np.sign(uu[tt,i,j]): #if in the same quadrant, check angle to uu-axis aka x-axis
								phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
								phi_uu=np.arccos(uu[tt,i,j]/uvbetrag[i][j])
								if phi_lon<phi_uu:
									alpha[tt,i,j] = beta
								else:
									alpha[tt,i,j] = -beta
	
						else: #3,4 quadrant
							if np.sign(dlon[j]) != np.sign(uu[tt,i,j]): # check if r und u are in same quadrant
								if dlon[j] >= 0:
									alpha[tt,i,j] = -beta
								else:
									alpha[tt,i,j] = beta
							if np.sign(dlon[j]) == np.sign(uu[tt,i,j]): #if in the same quadrant, check angle to uu-axis aka x-axis
								phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
								phi_uu=np.arccos(uu[tt,i,j]/uvbetrag[i][j])
								if phi_lon<phi_uu:
									alpha[tt,i,j] = -beta
								else:
									alpha[tt,i,j] = beta
					else:#u and r are above and below 0 line
						#check if 1 and 4 or 2 and 3:
						if np.sign(dlon[j]) == np.sign(uu[tt,i,j]):
							if np.sign(dlon[j]) >= 0: #1 and 4
								if np.sign(dlat[j]) <=0 and np.sign(vv[tt,i,j]): #lat below, vv above
									alpha[tt,i,j] = beta
								else:	#vv below lat above
									alpha[tt,i,j] = -beta
							else:
								if np.sign(dlat[j]) <=0 and np.sign(vv[tt,i,j]): ##lat below, vv above
									alpha[tt,i,j] = -beta
								else: #vv below lat above
									alpha[tt,i,j] = beta
						else: #1 and 3 or 2 and 4:
							if np.sign(dlon[j]) >=0 and np.sign(uu[tt,i,j]) <0 and np.sign(dlat[j]) > 0 and np.sign(vv[tt,i,j]) <0: #1 and 3 with r in 1 and uu in 3 
								phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
								phi_uu=np.arccos(uu[tt,i,j]/uvbetrag[i][j])
								if phi_uu+phi_lon <= np.pi:
									alpha[tt,i,j] = -beta
								else:
									alpha[tt,i,j] = beta
							elif np.sign(dlon[j]) <=0 and np.sign(uu[tt,i,j]) > 0 and np.sign(dlat[j]) < 0 and np.sign(vv[tt,i,j]) >0: #1 and 3 with r in 3 and uu in 1
								phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
								phi_uu=np.arccos(uu[tt,i,j]/uvbetrag[i][j])
								if phi_uu+phi_lon <= np.pi:
									alpha[tt,i,j] = beta
								else:
									alpha[tt,i,j] = -beta
							elif np.sign(dlon[j]) <=0 and np.sign(uu[tt,i,j]) > 0 and np.sign(dlat[j]) > 0 and np.sign(vv[tt,i,j]) <0: #2,4 r in 2, uu in 4
								phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
								phi_uu=np.arccos(uu[tt,i,j]/uvbetrag[i][j])
								if phi_uu+phi_lon <= np.pi:
									alpha[tt,i,j] = -beta
								else:
									alpha[tt,i,j] = beta
							elif np.sign(dlon[j]) >=0 and np.sign(uu[tt,i,j]) < 0 and np.sign(dlat[j]) < 0 and np.sign(vv[tt,i,j]) >0: #2,4 r in 4, uu in 2
								phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
								phi_uu=np.arccos(uu[tt,i,j]/uvbetrag[i][j])
								if phi_uu+phi_lon <= np.pi:
									alpha[tt,i,j] = beta
								else:
									alpha[tt,i,j] = -beta
	
					j+=1
				i+=1
			#print(alpha[t,:,:])
			u_parallel[tt,:,:]=uvbetrag*np.cos(alpha[tt,:,:])
			u_perpendicular[tt,:,:]=uvbetrag*np.sin(alpha[tt,:,:])
			tt+=1
	elif len(np.shape(uu)) == 2:
		
		dlon = lonc[1:]-lonc[:-1] #compute the lon/lat differences
		dlat = latc[1:]-latc[:-1]

		u_parallel = np.zeros(shape=(np.shape(uu)[0],np.shape(uu)[1]-1)) #define the new variables
		u_perpendicular = np.zeros(shape=np.shape(u_parallel))
		alpha = np.zeros(shape=np.shape(u_parallel)) #angle velocities have to be rotated

		latlonbetrag = np.sqrt(dlon*dlon+dlat*dlat) #
		#print(np.shape(latlonbetrag))

		uu = uu[:,:-1] #we lose the last value since we need to compute a difference
		vv = vv[:,:-1]
		distance = distance[1:]-distance[:-1]


		uvbetrag=np.sqrt(uu[:,:]*uu[:,:]+vv[:,:]*vv[:,:])
		i=0 
		while i < np.shape(u_parallel)[0]:
			j=0
			while j < np.shape(u_parallel)[1]:
				if vv[i,j] == 0:
					vv[i,j] = 1.e-10
				if uu[i,j] == 0:
					uu[i,j] = 1.e-10
				beta = np.arccos((uu[i,j]*dlon[j]+vv[i,j]*dlat[j])/(latlonbetrag[j]*uvbetrag[i][j]))
				if np.sign(dlat[j]) == np.sign(vv[i,j]): #check if r and u are both above or below 0
					#check if 1,2 or 3,4
					if np.sign(dlat[j]) >=0: #1,2 quadrant
						if np.sign(dlon[j]) != np.sign(uu[i,j]): # check if r und u are in same quadrant
							if dlon[j] >= 0:
								alpha[i,j] = beta
							else:
								alpha[i,j] = -beta
						if np.sign(dlon[j]) == np.sign(uu[i,j]): #if in the same quadrant, check angle to uu-axis aka x-axis
							phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
							phi_uu=np.arccos(uu[i,j]/uvbetrag[i][j])
							if phi_lon<phi_uu:
								alpha[i,j] = beta
							else:
								alpha[i,j] = -beta

					else: #3,4 quadrant
						if np.sign(dlon[j]) != np.sign(uu[i,j]): # check if r und u are in same quadrant
							if dlon[j] >= 0:
								alpha[i,j] = -beta
							else:
								alpha[i,j] = beta
						if np.sign(dlon[j]) == np.sign(uu[i,j]): #if in the same quadrant, check angle to uu-axis aka x-axis
							phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
							phi_uu=np.arccos(uu[i,j]/uvbetrag[i][j])
							if phi_lon<phi_uu:
								alpha[i,j] = -beta
							else:
								alpha[i,j] = beta
				else:#u and r are above and below 0 line
					#check if 1 and 4 or 2 and 3:
					if np.sign(dlon[j]) == np.sign(uu[i,j]):
						if np.sign(dlon[j]) >= 0: #1 and 4
							if np.sign(dlat[j]) <=0 and np.sign(vv[i,j]): #lat below, vv above
								alpha[i,j] = beta
							else:	#vv below lat above
								alpha[i,j] = -beta
						else:
							if np.sign(dlat[j]) <=0 and np.sign(vv[i,j]): ##lat below, vv above
								alpha[i,j] = -beta
							else: #vv below lat above
								alpha[i,j] = beta
					else: #1 and 3 or 2 and 4:
						if np.sign(dlon[j]) >=0 and np.sign(uu[i,j]) <0 and np.sign(dlat[j]) > 0 and np.sign(vv[i,j]) <0: #1 and 3 with r in 1 and uu in 3 
							phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
							phi_uu=np.arccos(uu[i,j]/uvbetrag[i][j])
							if phi_uu+phi_lon <= np.pi:
								alpha[i,j] = -beta
							else:
								alpha[i,j] = beta
						elif np.sign(dlon[j]) <=0 and np.sign(uu[i,j]) > 0 and np.sign(dlat[j]) < 0 and np.sign(vv[i,j]) >0: #1 and 3 with r in 3 and uu in 1
							phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
							phi_uu=np.arccos(uu[i,j]/uvbetrag[i][j])
							if phi_uu+phi_lon <= np.pi:
								alpha[i,j] = beta
							else:
								alpha[i,j] = -beta
						elif np.sign(dlon[j]) <=0 and np.sign(uu[i,j]) > 0 and np.sign(dlat[j]) > 0 and np.sign(vv[i,j]) <0: #2,4 r in 2, uu in 4
							phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
							phi_uu=np.arccos(uu[i,j]/uvbetrag[i][j])
							if phi_uu+phi_lon <= np.pi:
								alpha[i,j] = -beta
							else:
								alpha[i,j] = beta
						elif np.sign(dlon[j]) >=0 and np.sign(uu[i,j]) < 0 and np.sign(dlat[j]) < 0 and np.sign(vv[i,j]) >0: #2,4 r in 4, uu in 2
							phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
							phi_uu=np.arccos(uu[i,j]/uvbetrag[i][j])
							if phi_uu+phi_lon <= np.pi:
								alpha[i,j] = beta
							else:
								alpha[i,j] = -beta

				j+=1
			i+=1
		#print(alpha[t,:,:])
		u_parallel[:,:]=uvbetrag*np.cos(alpha[:,:])
		u_perpendicular[:,:]=uvbetrag*np.sin(alpha[:,:])
	
	elif len(np.shape(uu)) == 1:
		
		dlon = lonc[1:]-lonc[:-1] #compute the lon/lat differences
		dlat = latc[1:]-latc[:-1]

		u_parallel = np.zeros(shape=(np.shape(uu)[0]-1)) #define the new variables
		u_perpendicular = np.zeros(shape=np.shape(u_parallel))
		alpha = np.zeros(shape=np.shape(u_parallel)) #angle velocities have to be rotated

		latlonbetrag = np.sqrt(dlon*dlon+dlat*dlat) #
		#print(np.shape(latlonbetrag))

		uu = uu[:-1] #we lose the last value since we need to compute a difference
		vv = vv[:-1]
		distance = distance[1:]-distance[:-1]


		uvbetrag=np.sqrt(uu*uu+vv*vv)		
		j=0
		while j < np.shape(u_parallel)[0]:
			if vv[j] == 0:
				vv[j] = 1.e-10
			if uu[j] == 0:
				uu[j] = 1.e-10
			beta = np.arccos((uu[j]*dlon[j]+vv[j]*dlat[j])/(latlonbetrag[j]*uvbetrag[j]))
			if np.sign(dlat[j]) == np.sign(vv[j]): #check if r and u are both above or below 0
				#check if 1,2 or 3,4
				if np.sign(dlat[j]) >=0: #1,2 quadrant
					if np.sign(dlon[j]) != np.sign(uu[j]): # check if r und u are in same quadrant
						if dlon[j] >= 0:
							alpha[j] = beta
						else:
							alpha[j] = -beta
					if np.sign(dlon[j]) == np.sign(uu[j]): #if in the same quadrant, check angle to uu-axis aka x-axis
						phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
						phi_uu=np.arccos(uu[j]/uvbetrag[j])
						if phi_lon<phi_uu:
							alpha[j] = beta
						else:
							alpha[j] = -beta

				else: #3,4 quadrant
					if np.sign(dlon[j]) != np.sign(uu[j]): # check if r und u are in same quadrant
						if dlon[j] >= 0:
							alpha[j] = -beta
						else:
							alpha[j] = beta
					if np.sign(dlon[j]) == np.sign(uu[j]): #if in the same quadrant, check angle to uu-axis aka x-axis
						phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
						phi_uu=np.arccos(uu[j]/uvbetrag[j])
						if phi_lon<phi_uu:
							alpha[j] = -beta
						else:
							alpha[j] = beta
			else:#u and r are above and below 0 line
				#check if 1 and 4 or 2 and 3:
				if np.sign(dlon[j]) == np.sign(uu[j]):
					if np.sign(dlon[j]) >= 0: #1 and 4
						if np.sign(dlat[j]) <=0 and np.sign(vv[j]): #lat below, vv above
							alpha[j] = beta
						else:	#vv below lat above
							alpha[j] = -beta
					else:
						if np.sign(dlat[j]) <=0 and np.sign(vv[j]): ##lat below, vv above
							alpha[j] = -beta
						else: #vv below lat above
							alpha[j] = beta
				else: #1 and 3 or 2 and 4:
					if np.sign(dlon[j]) >=0 and np.sign(uu[j]) <0 and np.sign(dlat[j]) > 0 and np.sign(vv[j]) <0: #1 and 3 with r in 1 and uu in 3 
						phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
						phi_uu=np.arccos(uu[j]/uvbetrag[j])
						if phi_uu+phi_lon <= np.pi:
							alpha[j] = -beta
						else:
							alpha[j] = beta
					elif np.sign(dlon[j]) <=0 and np.sign(uu[j]) > 0 and np.sign(dlat[j]) < 0 and np.sign(vv[j]) >0: #1 and 3 with r in 3 and uu in 1
						phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
						phi_uu=np.arccos(uu[j]/uvbetrag[j])
						if phi_uu+phi_lon <= np.pi:
							alpha[j] = beta
						else:
							alpha[j] = -beta
					elif np.sign(dlon[j]) <=0 and np.sign(uu[j]) > 0 and np.sign(dlat[j]) > 0 and np.sign(vv[j]) <0: #2,4 r in 2, uu in 4
						phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
						phi_uu=np.arccos(uu[j]/uvbetrag[j])
						if phi_uu+phi_lon <= np.pi:
							alpha[j] = -beta
						else:
							alpha[j] = beta
					elif np.sign(dlon[j]) >=0 and np.sign(uu[j]) < 0 and np.sign(dlat[j]) < 0 and np.sign(vv[j]) >0: #2,4 r in 4, uu in 2
						phi_lon=np.arccos(dlon[j]/latlonbetrag[j])
						phi_uu=np.arccos(uu[j]/uvbetrag[j])
						if phi_uu+phi_lon <= np.pi:
							alpha[j] = beta
						else:
							alpha[j] = -beta

			j+=1
		u_parallel=uvbetrag*np.cos(alpha)
		u_perpendicular=uvbetrag*np.sin(alpha)

	return(u_parallel, u_perpendicular)