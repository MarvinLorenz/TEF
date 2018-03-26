import numpy as np
import netCDF4
import sys
from math import radians, cos, sin, asin, sqrt


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers is 6371
    m = 6371000* c #meters
    return m

def TEF_salt(ncfile,t,lat_or_lon,lat_or_lon_index,i_start,i_end,s,outputtype='normal',exchange_out=True,time_mean='no',land_depth=-10,variable_list=[],heat_exchange=False,rotation=True,slice_k = []): #args: time_mean
	
	#definitions of constants
	c_p = 4182.0 # J/(kg K)specific heat content of water

	#declare variables which change their names in the getm output if you use mean
	if len(variable_list)==0:
		if outputtype == 'mean':
			if lat_or_lon=='lon':
				varlist = ['saltmean','hmean','uumean','dyc','tempmean','sigma_t']
			elif lat_or_lon == 'lat':
				varlist = ['saltmean','hmean','vvmean','dxc','tempmean','sigma_t']
			elif lat_or_lon == 'transect':
				varlist = ['saltmean','hmean','uumean','vvmean','tempmean','sigma_t']
			else:
				sys.exit('lat_or_lon var is not known')
		elif outputtype == 'normal':
			if lat_or_lon=='lon':
				varlist = ['salt','h','uu','dyc','temp','sigma_t']
			elif lat_or_lon == 'lat':
				varlist = ['salt','h','vv','dxc','temp','sigma_t']
			elif lat_or_lon == 'transect':
				varlist = ['salt','h','uu','vv','temp','sigma_t']
			else:
				sys.exit('lat_or_lon var is not known')
				
		else:
			print('outputtype not known')
	else:
		print('using custom variable list')
		varlist = variable_list

	#check if t is an array or and index
	if type(t) is int:
		print('t is an index, nice')
	elif type(t) is list:
		if np.shape(t)[0] >1:
			print('t is an array, that means work')
		else:
			sys.exit('what kind of list is this?')

	elif type(t) is str:
		if t == 'all':
			print('okay, lets do all times')
		else:
			sys.exit('string not known') 

			
	########load necessary data depending on the input:
	print('loading data...')
	if type(t) is int:
		if lat_or_lon == 'lon':
			salt = ncfile.variables[varlist[0]][t,:,i_start:i_end+1,lat_or_lon_index]
			h = ncfile.variables[varlist[1]][t,:,i_start:i_end+1,lat_or_lon_index]
			vel = ncfile.variables[varlist[2]][t,:,i_start:i_end+1,lat_or_lon_index]
			ds = ncfile.variables[varlist[3]][i_start:i_end+1] if varlist[3] == 'latc' else ncfile.variables[varlist[3]][i_start:i_end+1,lat_or_lon_index] 
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][t,:,i_start:i_end+1,lat_or_lon_index]
				dens = ncfile.variables[varlist[5]][t,:,i_start:i_end+1,lat_or_lon_index]+1000.0
		elif lat_or_lon == 'lat':
			salt = ncfile.variables[varlist[0]][t,:,lat_or_lon_index,i_start:i_end+1]
			h = ncfile.variables[varlist[1]][t,:,lat_or_lon_index,i_start:i_end+1]
			vel = ncfile.variables[varlist[2]][t,:,lat_or_lon_index,i_start:i_end+1]
			ds = ncfile.variables[varlist[3]][i_start:i_end+1] if varlist[3] == 'lonc' else ncfile.variables[varlist[3]][lat_or_lon_index,i_start:i_end+1]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][t,:,lat_or_lon_index,i_start:i_end+1]
				dens = ncfile.variables[varlist[5]][t,:,lat_or_lon_index,i_start:i_end+1]+1000.0
		elif lat_or_lon == 'transect':
			ds = ncfile.variables['distance'][:]
			if rotation:
				lonc = ncfile.variables['lonc'][:]
				latc = ncfile.variables['latc'][:]
			salt = ncfile.variables[varlist[0]][t,1:,i_start:i_end+1] #bei nctransect from getm one level is added, therefore 1:
			h = ncfile.variables[varlist[1]][t,1:,i_start:i_end+1]
			uu = ncfile.variables[varlist[2]][t,1:,i_start:i_end+1]
			vv = ncfile.variables[varlist[3]][t,1:,i_start:i_end+1]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][t,1:,i_start:i_end+1]
				dens = ncfile.variables[varlist[5]][t,1:,i_start:i_end+1]+1000.0
			#sys.exit('to be done')
		else:
			print('lat_or_lon var is not known')
	elif type(t) is list:
		if lat_or_lon == 'lon':
			salt = ncfile.variables[varlist[0]][t[0]:t[1]+1,:,i_start:i_end+1,lat_or_lon_index]
			h = ncfile.variables[varlist[1]][t[0]:t[1]+1,:,i_start:i_end+1,lat_or_lon_index]
			vel = ncfile.variables[varlist[2]][t[0]:t[1]+1,:,i_start:i_end+1,lat_or_lon_index]
			ds = ncfile.variables[varlist[3]][i_start:i_end+1] if varlist[3] == 'latc' else ncfile.variables[varlist[3]][i_start:i_end+1,lat_or_lon_index]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][t[0]:t[1]+1,:,i_start:i_end+1,lat_or_lon_index]
				dens = ncfile.variables[varlist[5]][t[0]:t[1]+1,:,i_start:i_end+1,lat_or_lon_index]+1000.0
		elif lat_or_lon == 'lat':
			salt = ncfile.variables[varlist[0]][t[0]:t[1]+1,:,lat_or_lon_index,i_start:i_end+1]
			h = ncfile.variables[varlist[1]][t[0]:t[1]+1,:,lat_or_lon_index,i_start:i_end+1]
			vel = ncfile.variables[varlist[2]][t[0]:t[1]+1,:,lat_or_lon_index,i_start:i_end+1]
			ds = ncfile.variables[varlist[3]][i_start:i_end+1] if varlist[3] == 'lonc' else ncfile.variables[varlist[3]][lat_or_lon_index,i_start:i_end+1]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][t[0]:t[1]+1,:,lat_or_lon_index,i_start:i_end+1]
				dens = ncfile.variables[varlist[5]][t[0]:t[1]+1,:,lat_or_lon_index,i_start:i_end+1]+1000.0
		elif lat_or_lon == 'transect':
			ds = ncfile.variables['distance'][:]
			if rotation:
				lonc = ncfile.variables['lonc'][:]
				latc = ncfile.variables['latc'][:]
			salt = ncfile.variables[varlist[0]][t[0]:t[1]+1,1:,i_start:i_end+1]
			h = ncfile.variables[varlist[1]][t[0]:t[1]+1,1:,i_start:i_end+1]
			uu = ncfile.variables[varlist[2]][t[0]:t[1]+1,1:,i_start:i_end+1]
			vv = ncfile.variables[varlist[3]][t[0]:t[1]+1,1:,i_start:i_end+1]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][t[0]:t[1]+1,1:,i_start:i_end+1]
				dens = ncfile.variables[varlist[5]][t[0]:t[1]+1,1:,i_start:i_end+1]+1000.0
			#sys.exit('to be done')
		else:
			print('lat_or_lon var is not known')
	elif t == 'all':
		if lat_or_lon == 'lon':
			salt = ncfile.variables[varlist[0]][:,:,i_start:i_end+1,lat_or_lon_index]
			h = ncfile.variables[varlist[1]][:,:,i_start:i_end+1,lat_or_lon_index]
			vel = ncfile.variables[varlist[2]][:,:,i_start:i_end+1,lat_or_lon_index]
			ds = ncfile.variables[varlist[3]][i_start:i_end+1] if varlist[3] == 'latc' else ncfile.variables[varlist[3]][i_start:i_end+1,lat_or_lon_index]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][:,:,i_start:i_end+1,lat_or_lon_index]
				dens = ncfile.variables[varlist[5]][:,:,i_start:i_end+1,lat_or_lon_index]+1000.0
		elif lat_or_lon == 'lat':
			salt = ncfile.variables[varlist[0]][:,:,lat_or_lon_index,i_start:i_end+1]
			h = ncfile.variables[varlist[1]][:,:,lat_or_lon_index,i_start:i_end+1]
			vel = ncfile.variables[varlist[2]][:,:,lat_or_lon_index,i_start:i_end+1]
			ds = ncfile.variables[varlist[3]][i_start:i_end+1] if varlist[3] == 'lonc' else ncfile.variables[varlist[3]][lat_or_lon_index,i_start:i_end+1]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][:,:,lat_or_lon_index,i_start:i_end+1]
				dens = ncfile.variables[varlist[5]][:,:,lat_or_lon_index,i_start:i_end+1]+1000.0
		elif lat_or_lon == 'transect':
			ds = ncfile.variables['distance'][:]
			if rotation:
				lonc = ncfile.variables['lonc'][:]
				latc = ncfile.variables['latc'][:]
			salt = ncfile.variables[varlist[0]][:,1:,i_start:i_end+1]
			h = ncfile.variables[varlist[1]][:,1:,i_start:i_end+1]
			uu = ncfile.variables[varlist[2]][:,1:,i_start:i_end+1]
			vv = ncfile.variables[varlist[3]][:,1:,i_start:i_end+1]
			if heat_exchange:
				temp = ncfile.variables[varlist[4]][:,1:,i_start:i_end+1]
				dens = ncfile.variables[varlist[5]][:,1:,i_start:i_end+1]+1000.0
			#sys.exit('to be done')
		else:
			print('lat_or_lon var is not known')
	print('data is loaded')
	
	if len(slice_k) != 0:
		salt = salt[:,slice_k[0]:slice_k[1],:]
		vel = vel[:,slice_k[0]:slice_k[1],:]
		h = h[:,slice_k[0]:slice_k[1],:]
		if heat_exchange:
			temp = temp[:,slice_k[0]:slice_k[1],:]
			dens = dens[:,slice_k[0]:slice_k[1],:]
	#print(ds)

	#check if ds is lonc or latc:
	if variable_list != []:
		if variable_list[3] == 'lonc':
			print('doing haversine function')
			latc = ncfile.variables['latc'][lat_or_lon_index]
			dist_betw_2_p = haversine(ds[1],latc,ds[0],latc)
			ds[:] = dist_betw_2_p
			print(ds)
		elif variable_list[3] == 'latc':
			print('doing haversine function')
			lonc = ncfile.variables['lonc'][lat_or_lon_index]
			dist_betw_2_p = haversine(lonc,ds[1],lonc,ds[0])
			ds[:] = dist_betw_2_p
			print(ds)

	#do a land check:
	if lat_or_lon == 'lon':
		bathy = np.ma.getdata(ncfile.variables['bathymetry'][i_start:i_end+1,lat_or_lon_index])
		do_check=True
	elif lat_or_lon == 'lat':
		bathy = np.ma.getdata(ncfile.variables['bathymetry'][lat_or_lon_index,i_start:i_end+1])
		do_check=True
	elif lat_or_lon == 'transect':
		bathy = np.ma.getdata(ncfile.variables['bathymetry'][i_start:i_end+1])
		do_check=True
	else:
		print('WARNING: no bathymetry for land check found')
		do_check=False

	if do_check:
		print(bathy)
		i=0
		while i < len(bathy):
			if bathy[i] == land_depth:
				sys.exit('land point found, ignoring land not implemented yet')
			else:
				i+=1
		print('no land point found')

	#####calculations

	##### if transect: rotate velocities

	#rotate velocities:
	if lat_or_lon == 'transect':
		if not type(t) is int and rotation:

			print('doing rotation')
			dlon = lonc[1:]-lonc[:-1] #compute the lon/lat differences
			dlat = latc[1:]-latc[:-1]

			u_new = np.zeros(shape=(np.shape(uu)[0],np.shape(uu)[1],np.shape(uu)[2]-1)) #define the new variables
			vel = np.zeros(shape=np.shape(u_new))
			alpha = np.zeros(shape=np.shape(u_new)) #angle velocities have to be rotated

			latlonbetrag = np.sqrt(dlon*dlon+dlat*dlat) #
			#print(np.shape(latlonbetrag))

			uu = uu[:,:,:-1] #we lose the last value since we need to compute a difference
			vv = vv[:,:,:-1]
			ds = ds[1:]-ds[:-1]

			#####checking if rotation is right
			#uu[0,0,0] = 1
			#vv[0,0,0] = -1
			#print(uu[0,0,0],vv[0,0,0])
			#print(dlon[0],dlat[0])
			#####

			tt = 0
			while tt < np.shape(uu)[0]:
				uvbetrag=np.sqrt(uu[tt,:,:]*uu[tt,:,:]+vv[tt,:,:]*vv[tt,:,:])
				i=0 
				while i < np.shape(u_new)[1]:
					j=0
					while j < np.shape(u_new)[2]:
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
				u_new[tt,:,:]=uvbetrag*np.cos(alpha[tt,:,:])
				vel[tt,:,:]=uvbetrag*np.sin(alpha[tt,:,:])
				tt+=1
			salt=salt[:,:,:-1] 
			if heat_exchange:
				temp=temp[:,:,:-1] 
				dens=dens[:,:,:-1] 
			h=h[:,:,:-1] #slice salt to mathc other parameters
			#print(u_new[0,0,0],vel[0,0,0])
			print('rotation done')

		elif rotation:
			print('doing rotation')
			dlon = lonc[1:]-lonc[:-1] #compute the lon/lat differences
			dlat = latc[1:]-latc[:-1]

			u_new = np.zeros(shape=(np.shape(uu)[0],np.shape(uu)[1]-1)) #define the new variables
			vel = np.zeros(shape=np.shape(u_new))
			alpha = np.zeros(shape=np.shape(u_new)) #angle velocities have to be rotated

			latlonbetrag = np.sqrt(dlon*dlon+dlat*dlat) #
			#print(np.shape(latlonbetrag))

			uu = uu[:,:-1] #we lose the last value since we need to compute a difference
			vv = vv[:,:-1]
			ds = ds[1:]-ds[:-1]


			uvbetrag=np.sqrt(uu[:,:]*uu[:,:]+vv[:,:]*vv[:,:])
			i=0 
			while i < np.shape(u_new)[0]:
				j=0
				while j < np.shape(u_new)[1]:
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
			u_new[:,:]=uvbetrag*np.cos(alpha[:,:])
			vel[:,:]=uvbetrag*np.sin(alpha[:,:])
			salt=salt[:,:-1] #slice salt to fit other parameters
			if heat_exchange:
				temp=temp[:,:-1]
				dens=dens[:,:-1]
			h=h[:,:,:-1]
			print('rotation done')

		else:
			print('else rot')
			vel=uu

	##### time mean?!

	if time_mean != 'no' and type(time_mean) is int:

		print('doing time average over', time_mean, 'timesteps')
		print('time dimension before:', np.shape(salt)[0])
		print('time dimension after:', int(np.shape(salt)[0]/time_mean))
		#create new variables with smaller time dimension:
		salt_new = np.zeros(shape=(int(np.shape(salt)[0]/time_mean),np.shape(salt)[1],np.shape(salt)[2])) #using int to cut the last timesteps if they are not enough a compute a full mean
		vel_new = np.zeros(shape=np.shape(salt_new))
		if heat_exchange:
			temp_new = np.zeros(shape=np.shape(salt_new))
			dens_new = np.zeros(shape=np.shape(salt_new))
		h_new = np.zeros(shape=np.shape(salt_new))
		time_indices = []
		i=0
		while i < int(np.shape(salt)[0]/time_mean):
			time_indices.append((i+1)*time_mean-1) #time index of the last timestep 
			#temporary salt and vel variable
			salt_temp = salt[i:(i+1)*time_mean,:,:]*h[i:(i+1)*time_mean,:,:]
			if heat_exchange:
				temp_temp = temp[i:(i+1)*time_mean,:,:]*h[i:(i+1)*time_mean,:,:]
				dens_temp = dens[i:(i+1)*time_mean,:,:]*h[i:(i+1)*time_mean,:,:]
			vel_temp = vel[i:(i+1)*time_mean,:,:]*h[i:(i+1)*time_mean,:,:]
			h_temp = h[i:(i+1)*time_mean,:,:]
			#now store the correct values into new variables
			salt_new[i,:,:] = np.sum(salt_temp, axis=0)/np.sum(h_temp, axis=0) #with h weighted mean 
			if heat_exchange:
				temp_new[i,:,:] = np.sum(temp_temp, axis=0)/np.sum(h_temp, axis=0)
				dens_new[i,:,:]= np.sum(dens_temp, axis=0)/np.sum(h_temp, axis=0)
			vel_new[i,:,:] = np.sum(vel_temp, axis=0)/np.sum(h_temp, axis=0)
			h_new[i,:,:] = np.mean(h_temp, axis=0)
			i+=1 
		#rename the variables
		if heat_exchange:
			temp = temp_new
			dens = dens_new
		salt = salt_new
		vel = vel_new

	if time_mean == 'all':

		print('doing time average over', time_mean, 'timesteps')
		#create new variables with smaller time dimension:
		salt_new = np.zeros(shape=(np.shape(salt)[1],np.shape(salt)[2])) #using int to cut the last timesteps if they are not enough a compute a full mean
		if heat_exchange:
			temp_new = np.zeros(shape=np.shape(salt_new))
			dens_new = np.zeros(shape=np.shape(salt_new))
		vel_new = np.zeros(shape=np.shape(salt_new))
		h_new = np.zeros(shape=np.shape(salt_new))
		time_indices = []
		
		#now store the correct values into new variables
		salt_new[:,:] = np.sum(salt*h, axis=0)/np.sum(h, axis=0)
		if heat_exchange:
			temp_new[:,:] = np.sum(temp*h, axis=0)/np.sum(h, axis=0) #with h weighted mean 
			dens_new[:,:] = np.sum(dens*h, axis=0)/np.sum(h, axis=0)
		vel_new[:,:] = np.sum(vel*h, axis=0)/np.sum(h, axis=0)
		h_new[:,:] = np.mean(h, axis=0)

		#change t so it is a 'timestep again' so it goes into the right routine!
		t = 0
	#rename the variables
		if heat_exchange:
			temp = temp_new
			dens = dens_new
		salt = salt_new
		vel = vel_new
		h = h_new 
	#print(np.shape(salt)[0])


	##### do the TEF
	#define the variables which should be calculated

	if not type(t) is int: #if t is not an index, we have a higher dimension

		Q_iso = np.zeros(shape=(np.shape(salt)[0],len(s)))
		dqds = np.zeros(shape=(np.shape(salt)[0],len(s)-1)) #-1 because we lose 1 value for s if we do the derivative
		if heat_exchange:
			heat_cont = np.zeros(shape=np.shape(Q_iso)) #unit: J m-3
		#print(np.shape(salt))
		#print(np.shape(Q_iso))

		#calculate Q_iso
		tt=0
		while tt < np.shape(salt)[0]: #time
			ss = 0
			while ss < len(s)-1: #salinity
				Q_temp = []
				if heat_exchange:
					Heat_temp = []
					layer_height_temp = []
				i=0
				while i < np.shape(salt)[2]: #lon/lat
					k=0
					while k < np.shape(salt)[1]: # depth
						if heat_exchange:
							if salt[tt,k,i]>=s[ss]-0.5*(np.abs(s[ss+1])-np.abs(s[ss])) and salt[tt,k,i]<s[ss]+0.5*(np.abs(s[ss+1])-np.abs(s[ss])): #only add heat if it is within in the intervall s-0.5 ds <= s < s+0.5ds
								Heat_temp.append(h[tt,k,i]*temp[tt,k,i]*dens[tt,k,i]*c_p)
								layer_height_temp.append(h[tt,k,i])

						if salt[tt,k,i] >s[ss]: #add transport if salinity > than s[ss] which is the new axis
							Q_temp.append(h[tt,k,i]*vel[tt,k,i]*ds[i])
							k+=1
						else:
							k+=1
					i+=1
				Q_iso[tt,ss] = np.sum(Q_temp) #definition of Q_iso
				if heat_exchange:
					if np.sum(layer_height_temp) != 0:
						heat_cont[tt,ss] = np.sum(Heat_temp)/np.sum(layer_height_temp)
					else:
						heat_cont[tt,ss] = np.sum(Heat_temp)
				ss+=1
			tt+=1

		# compute the derivative
		# define a by one shorter salinity list (we lose one point when we compute the derivative)
		s_new = s[0:(len(s)-1)]
		#print(len(s_new))
		tt = 0
		while tt < np.shape(salt)[0]:
			dqds[tt,:] = -(Q_iso[tt,1:]-Q_iso[tt,0:(len(s)-1)])/(s[1:]-s[0:(len(s)-1)]) #computation of q = -dQ/ds
			tt+=1

		if not exchange_out:
			if time_mean == 'no':
				return(Q_iso,dqds,s_new)
			else:
				return(Q_iso,dqds,s_new,time_indices)

		if exchange_out:

			#calculate exchange flow properties

			#this works at the moment only if s has constant step sizes
			Q_in=np.zeros(shape=(np.shape(dqds)[0],))
			Q_out=np.zeros(shape=(np.shape(dqds)[0],))
			S_in=np.zeros(shape=(np.shape(dqds)[0],))
			S_out=np.zeros(shape=(np.shape(dqds)[0],))
			s_in=np.zeros(shape=(np.shape(dqds)[0],))
			s_out=np.zeros(shape=(np.shape(dqds)[0],))

			tt=0
			while tt < np.shape(dqds)[0]:

				Q_in[tt] = np.sum((dqds[tt,:]*(s_new[1]-s_new[0])).clip(min=0))
				Q_out[tt] = np.sum((dqds[tt,:]*(s_new[1]-s_new[0])).clip(max=0))
				S_in[tt] = np.sum((dqds[tt,:]*(s_new[1]-s_new[0])*s_new).clip(min=0)) # Salinity flux
				S_out[tt] = np.sum((dqds[tt,:]*(s_new[1]-s_new[0])*s_new).clip(max=0)) # Salinity flux
				s_in[tt] = S_in[tt]/Q_in[tt]
				s_out[tt] = S_out[tt]/Q_out[tt]

				tt+=1

			ex_out = [Q_in, Q_out, S_in, S_out, s_in, s_out]
			if not heat_exchange:
				if time_mean == 'no':
					return(Q_iso,dqds,s_new, ex_out)
				else:
					return(Q_iso,dqds,s_new, ex_out, time_indices)
			if heat_exchange:
				#calculate stuff:
				Heat_in=np.zeros(shape=(np.shape(dqds)[0],))
				Heat_out=np.zeros(shape=(np.shape(dqds)[0],))
				heat_new = heat_cont[:,:-1]
				tt=0
				while tt < np.shape(dqds)[0]:
					Heat_in[tt] = np.sum((dqds[tt,:]*(s_new[1]-s_new[0])*heat_new).clip(min=0)) # J/s
					Heat_out[tt] = np.sum((dqds[tt,:]*(s_new[1]-s_new[0])*heat_new).clip(max=0)) # J/s
					tt+=1
				ex_out.append(Heat_in)
				ex_out.append(Heat_out)

				if time_mean == 'no':
					return(Q_iso,dqds,s_new, ex_out)
				else:
					return(Q_iso,dqds,s_new, ex_out, time_indices)


	else:
		Q_iso = np.zeros(shape=(len(s),))
		dqds = np.zeros(shape=(len(s)-1,)) #-1 because we lose 1 value for s if we do the derivative
		if heat_exchange:
			heat_cont = np.zeros(shape=np.shape(Q_iso)) #unit: J m-3
		#print(np.shape(salt))
		#print(np.shape(Q_iso))

		#computation:
		ss = 0
		while ss < len(s)-1: #salinity -1 is needed for heat exchange. This should not change anything since Q_iso is zero anyways for great salinities
			Q_temp = []
			if heat_exchange:
				Heat_temp = []
				layer_height_temp = []
			i=0
			while i < np.shape(salt)[1]: #lon/lat
				k=0
				while k < np.shape(salt)[0]: # depth
					if heat_exchange:
						if salt[k,i]>=s[ss]-0.5*(np.abs(s[ss+1])-np.abs(s[ss])) and salt[k,i]<s[ss]+0.5*(np.abs(s[ss+1])-np.abs(s[ss])): #only add heat if it is within in the intervall s-0.5 ds <= s < s+0.5ds
							Heat_temp.append(h[k,i]*temp[k,i]*dens[k,i]*c_p)
							layer_height_temp.append(h[k,i])
					if salt[k,i] >s[ss]: #add transport if salinity > than s[ss] which is the new axis
						Q_temp.append(h[k,i]*vel[k,i]*ds[i])
						k+=1
					else:
						k+=1
				i+=1
			Q_iso[ss] = np.sum(Q_temp) #definition of Q_iso
			if heat_exchange:
				if np.sum(layer_height_temp) != 0:
					heat_cont[ss] = np.sum(Heat_temp)/np.sum(layer_height_temp)
				else:
					heat_cont[ss] = np.sum(Heat_temp)
			ss+=1
		#print(Q_iso)

		# compute the derivative
		# define a by one shorter salinity list (we lose one point when we compute the derivative)
		s_new = s[0:(len(s)-1)]
		#print(len(s_new))
		dqds[:] = -(Q_iso[1:]-Q_iso[0:(len(s)-1)])/(s[1:]-s[0:(len(s)-1)]) #computation of q = -dQ/ds
		#print(dqds)

		if not exchange_out:
			if time_mean == 'no':
				return(Q_iso,dqds,s_new)
			else:
				return(Q_iso,dqds,s_new,time_indices)

		if exchange_out:

			#calculate exchange flow properties

			#this works at the moment only if s has constant step sizes

			Q_in = np.sum((dqds*(s_new[1]-s_new[0])).clip(min=0))
			Q_out = np.sum((dqds*(s_new[1]-s_new[0])).clip(max=0))
			S_in = np.sum((dqds*(s_new[1]-s_new[0])*s_new).clip(min=0)) # Salinity flux
			S_out = np.sum((dqds*(s_new[1]-s_new[0])*s_new).clip(max=0)) # Salinity flux
			s_in = S_in/Q_in
			s_out = S_out/Q_out

			ex_out = [Q_in, Q_out, S_in, S_out, s_in, s_out]
			if not heat_exchange:
				if time_mean == 'no':
					return(Q_iso,dqds,s_new, ex_out)
				else:
					return(Q_iso,dqds,s_new, ex_out,time_indices)
			if heat_exchange:
				heat_new = heat_cont[:-1]
				Heat_in = np.sum((dqds*(s_new[1]-s_new[0])*heat_new).clip(min=0))
				Heat_out = np.sum((dqds*(s_new[1]-s_new[0])*heat_new).clip(max=0))
				ex_out.append(Heat_in)
				ex_out.append(Heat_out)
				if time_mean == 'no':
					return(Q_iso,dqds,s_new, ex_out)
				else:
					return(Q_iso,dqds,s_new, ex_out,time_indices)
