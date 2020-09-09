#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:37:58 2020

@author: noahbrauer
"""

import h5py 
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
import pyart

file = '2A.GPM.DPR.V8-20180723.20150620-S043651-E060924.007430.V06A.HDF5'


#Read in file, extract lat, lons, PSD

DPR = h5py.File(file, 'r')

#nscan is the nummber of scans in each granule 
#nray is the number of angle bins in each  scan; think of these as footprint scans (5 km in diameter for each footprint)
#nbin is the number of range bins in each ray (angle bins)
#nDSD: Parameters are N0 (number concentration) and D0 (mean drop diameter)


lat = DPR['NS']['Latitude'][:,:]    
lon = DPR['NS']['Longitude'][:,:]
z = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:] #nscan x nray (7934,49)
dsd = DPR['NS']['SLV']['paramDSD'][:] #nscan x nray x nbin x DSDmoment (7934,49,176,2)
precip = DPR['NS']['SLV']['precipRateNearSurface'][:]
freezing = DPR['NS']['VER']['heightZeroDeg'][:]
precip_type = DPR['NS']['CSF']['typePrecip'][:]/10**7

###Define cross-section length

ind1 = np.where((lon[:,0]>=-90.5)) #Where are the DPR lons >= -100
ind2 = np.where((lon[:,0])<=-85.5) #Where are lons <= -85
ind3 = np.intersect1d(ind1,ind2) #Both conditions need to be satisfied here


#Change the ray to change cross-section location (i.e. the "27" in this case)
#%%

###Setup to 2D grid for plotting

x = 2.* 17 #48 degrees (from -17 to 17)
re = 6378. #radius of the earth
theta = -1*(x/2.) + (x/48.)*np.arange(0,49) #Split into equal degrees (from -17 to 17)
theta2  = np.ones(theta.shape[0]+1)*np.nan #Define an empty array (NaNs) with same shape as ray dimension
theta = theta - 0.70833333/2. #Shift degrees for plotting pruposes
theta2[:-1] = theta #remove last array dimension in theta so python won't get angry about shape
theta2[-1] = theta[-1] + 0.70833333
theta = theta2*(np.pi/180.) #Convert from degrees to radians

prh = np.ones((177,50))*np.nan #Define empty grid

for i in range(prh.shape[0]): #Loop over each range gate
    for j in range(prh.shape[1]): #Loop over each scan
            a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #Orbit height of 407 km 
            
            prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a)
            
h2 = prh #where prh is the (range bins,ray)-space
h3 =np.ones((177,50))*np.nan

for i in range(h3.shape[1]):
    h3[:,i] = h2[::-1,i] #This reverses the vertical dimension so as indices increase, height increases

#%%
    
ku = DPR['NS']['SLV']['zFactorCorrected'][ind3,:,:] #Read in ku-band reflectivity; nscan x nray (554,49,176)
n0 = dsd[ind3,:,:,0] #Read in the number concentration
d0 = dsd[ind3,:,:,1] #Read in the mean drop diameter  #Both have dimensions nscan x nray x nbin (554,49,176)
zeroDeg = freezing[ind3,:]
#Cut all parameters so they are at same ray as above
ray = 34

ku = ku[:,ray,:]
n0 = n0[:,ray,:]
d0 = d0[:,ray,:]
zero_deg_isotherm = zeroDeg[:,ray]/1000 #Convert from meters to km

#Take lats and lons along same ray
lons = DPR['NS']['Longitude'][ind3,ray]
lats = DPR['NS']['Latitude'][ind3,ray]




#Choose a starting point, then calculate distance
lat0 = lats[0]
lon0 = lons[0]


p = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0) #Define a projection and set starting lat an lon to same point as above

#Define a 2D array for plotting purposes

lat_3d = np.ones(ku.shape)*np.nan
lon_3d = np.ones(ku.shape)*np.nan

for i in range(ku.shape[0]):
    lat_3d[i,:] = lats[i] 
    lon_3d[i,:] = lons[i]  
        

x,y = p(lon_3d,lat_3d) #Now convert degrees to distance (in km)
R_gpm = np.sqrt(x**2 + y**2)*np.sign(x) #Keeps sign of number; converts to radial distance 

#Reverse range gate order for all parameters

ku = ku[:,::-1]
n0 = n0[:,::-1]
d0 = d0[:,::-1]



ku = np.ma.masked_where(ku<=12, ku) #Mask all the bad points in ku data
y = np.ones([ku.shape[0], ku.shape[1]]) #Define an empty array

#Define the appropriate range bins
h4 = h3[:,ray] #This number MUST match the same ray being used
for i in range(y.shape[1]):
    y[:,i] = h4[i]


#Remove the values less than or equal to zero

n0_nan = np.ones((168,176))*np.nan
d0_nan = np.ones((168,176))*np.nan

for i in range(n0.shape[0]):
    for j in range(n0.shape[1]):
        
        if n0[i,j] <= 0:
            n0_nan[i,j] = np.nan
        elif d0[i,j] <= 0:
            d0[i,j] = np.nan
        else:
            n0_nan[i,j] = n0[i,j]
            d0_nan[i,j] = d0[i,j]
            



precip_nan = np.ones((7934,49))*np.nan

for i in range(precip_nan.shape[0]):
    for j in range(precip_nan.shape[1]):
        
        if precip[i,j] == 0:
            precip_nan[i,j] = np.nan
            
        else:
            precip_nan[i,j] = precip[i,j] 








    
    
#%% 
    
#Now we plot! #N-S (along track first)  


#Plot mean drop size first


plt.figure(figsize=(10,10))

vmax = 3
vmin = 0

R_min = R_gpm.min()
R_max = R_gpm.max()

label_size = 30
title_size = 32


cmin =0.; cmax = 3.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, d0_nan, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k', label = r'$0^{o}$C isotherm')
plt.xlabel('Along Track Distance (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)
plt.title(r'GPM Overpass 6/20 0436 UTC $D_{M}$ ', size = title_size)
plt.xlim(0,250)
plt.ylim(0,15)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[mm]',size = label_size)
plt.clim(0,3)
plt.xticks(fontsize = label_size)
plt.yticks(fontsize = label_size)

plt.show()    
    


####Number concentration (liquid water content)

plt.figure(figsize=(10,10))

vmax = 70
vmin = 0

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =0.; cmax = 70.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, n0_nan, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k', label = r'$0^{o}$C isotherm')
plt.xlabel('Along Track Distance (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)
plt.title(r'GPM Overpass 6/20 0436 UTC $N_{W}$ ', size = title_size)
plt.xlim(0,250)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[mm $m^{-3}$]',size = label_size)
plt.clim(0,80)
plt.xticks(fontsize = label_size)
plt.yticks(fontsize = label_size)

plt.show()    
    

###And lastly Ku-band


plt.figure(figsize=(10,10))

vmax = 60
vmin =12

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, ku, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k', label = r'$0^{o}$C isotherm')
plt.xlabel('Along Track Distance (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)
plt.title(r'GPM Overpass 6/20 0436 UTC KuPR ', size = title_size)
#plt.xlim(300,450)
plt.xlim(0,250)
plt.ylim(0,15)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.clim(12,60)
plt.xticks(fontsize = label_size)
plt.yticks(fontsize = label_size)

plt.show()   






#Determine median value of ind3 (intersection between longitudes)
cross_track_index = int(len(ind3)/2)+ind3[0]



#Let's plot a map with the center point location for our cross-sections

#Mask near surface reflectivity values

z = np.ma.masked_where(z<=12, z)



cmin = 12.; cmax = 70.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-91,-85]); ylim = np.array([35,40])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,z,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray]+0.65,lat[ind3[0],ray],'*w',markersize = 20, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[dBZ]',name='Calibri',size=label_size)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(label_size)
   
plt.title('GPM Overpass 6/20 0436 UTC', size = title_size)



#%%

#Surface rainfall rate

plt.figure(figsize=(10,10))

cmin = 0.; cmax = 50.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
xlim = np.array([-91,-85]); ylim = np.array([35,40])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,precip_nan,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray]+0.65,lat[ind3[0],ray],'*w',markersize = 20, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[mm/hour]',name='Calibri',size=label_size)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(label_size)
plt.title('GPM Overpass 6/20 0436 UTC Surface Rainfall Rate', size = title_size)

#%%

#Set categories for precip type:

precip_type[precip_type<1] = np.nan





colormap = ['white','green', 'yellow', 'red']
 
from matplotlib.colors import ListedColormap
cmap_precip = ListedColormap(colormap)     

print(np.nanmax(precip_type))

#%%

plt.figure(figsize=(10,10))

cmin = 1.; cmax = 4.; cint = 1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_precip,lut=nlevs)
xlim = np.array([-91,-85]); ylim = np.array([35,40])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,precip_type,clevs,cmap=cmap_precip,extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray]+0.65,lat[ind3[0],ray],'*w',markersize = 20, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_yticklabels(['None','Stratiform', 'Convective', 'Other'], size = label_size)

plt.title('GPM Overpass 6/20 0436 UTC Precipitation Category', size = title_size)
   
         


#%%

#Now lets plot CFADs

def create_histogram(input_values, num_bins, min_value, max_value):
    """Creates a histogram with uniform bin-spacing.
    N = number of input values
    K = number of bins
    :param input_values: length-N numpy array of input values (to be binned).
    :param num_bins: Number of bins.
    :param min_value: Minimum value to include in histogram.  Any input value <
        `min_value` will be assigned to the first bin.
    :param max_value: Maximum value to include in histogram.  Any input value >
        `max_value` will be assigned to the last bin.
    :return: input_to_bin_indices: length-N numpy array of bin indices.  If
        input_values[i] = j, the [i]th input value belongs in the [j]th bin.
    :return: num_examples_by_bin: length-K numpy array, where the [j]th value is
        the number of inputs assigned to the [j]th bin.
    """

    bin_cutoffs = np.linspace(min_value, max_value, num=num_bins + 1)
    input_to_bin_indices = np.digitize(
        input_values, bin_cutoffs, right=False) - 1
    input_to_bin_indices[input_to_bin_indices < 0] = 0
    input_to_bin_indices[input_to_bin_indices > num_bins - 1] = num_bins - 1

    num_examples_by_bin = np.full(num_bins, -1, dtype=int)
    for j in range(num_bins):
        num_examples_by_bin[j] = np.sum(input_to_bin_indices == j)

    return input_to_bin_indices, num_examples_by_bin    

#Convert the NaNs to zeros 
    

ku_nonan = np.nan_to_num(ku)
print(np.nanmax(ku_nonan))




bins = 28

ku_hist = np.zeros((ku_nonan.shape[1],bins))


for i in range(ku.shape[0]):
    ku_hist[i,:] = create_histogram(ku_nonan[:,i],28,18,46 )[1]
    
    
print(y.shape)    
print(R_gpm.shape)    
print(ku_hist.T.shape)
print(ku.shape)


#%%
#Now convert 0 dBZ back to NaNs
ku_hist[ku_hist<=12] =np.nan

num_bins = ku_hist.shape[1]
num_altitudes = ku_hist.shape[1]

#Now we can plot

fig,ax = plt.subplots(figsize = [10,5])
pm = plt.pcolormesh(ku_hist, cmap = 'pyart_NWSRef')

x_tick_values = np.linspace(0, num_bins - 1, num = num_bins, dtype = float)
x_tick_strings = [None] * num_bins


#num_altitudes = np.arange(0,44.25,step = 0.25)

min_z = 18

for k in range(num_bins):
    if k == 0:
        x_tick_strings[k] = '< 18'
    elif k == num_bins - 1:
        x_tick_strings[k] = '>= 46'
    else:
        x_tick_strings[k] = '[{0:d}, {1:d})'.format(k + min_z, (k + min_z) +1 )

plt.xticks(x_tick_values[::2], x_tick_strings[::2], rotation = 90, size = 16) 

plt.ylim(0,100)
ax.set_yticks([i for i in range(0,101,4)])
ax.set_yticklabels([i/4 for i in range(0,101,4)])
  


plt.ylabel('Altitude (km)', size = 18)

plt.xlim(0,25)

plt.clim(0,5)
plt.xlabel('dBZ', size = 18)
plt.title(r'GPM $Z_{M}$(Ku) 6/20 0436 UTC', size = 20)
plt.colorbar().set_label(label = 'Frequency', size = 26)
#plt.ylim(0,20)
plt.show()








    









####Now plot the cross-track cross-section

#%%



#now lets do the ~W-E cross-section 
ku = DPR['NS']['SLV']['zFactorCorrected'][cross_track_index,:,:]
lons = DPR['NS']['Longitude'][cross_track_index,:]
lats = DPR['NS']['Latitude'][cross_track_index,:]
# pick the starting point to be the reference point to calc dist. 
lat0 = lats[0]
lon0 = lons[0]

p = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0)
#make it 2d for the pcolormesh 
lat_3d = np.zeros(ku.shape)
lon_3d = np.zeros(ku.shape)
for i in np.arange(0,ku.shape[0]):
    lat_3d[i,:] = lats[i]
    lon_3d[i,:] = lons[i]
    
#convert to distances 
x,y = p(lon_3d,lat_3d)
#get the radial distance
R_gpm = np.sqrt(x**2  + y**2)*np.sign(x)
#flip the order of the range gates
ku = ku[:,::-1]

#mask bad data
ku = np.ma.masked_where(ku <= 12,ku)


y = h3.T

plt.figure(figsize=(10,10))
R_min = R_gpm.min()
R_max = R_gpm.max()


cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., h3[1:,1:].T, ku, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
plt.xlabel('Cross-Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 6/20 0436 UTC Ku-Band ', size = 20)
#plt.xlim(300,450)
plt.xlim(0,400)
plt.ylim(0,15)
#plt.colorbar(pm, label = 'dBZ')
plt.colorbar().set_label(label='dBZ',size=18)
plt.clim(12,60)


plt.show()   

