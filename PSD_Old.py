#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:03:51 2020

Source code: Randy Chase (@DopplerChase); Additions and Modifications by noahbrauer (@NOAABrauer)
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

file = '2A.GPM.DPR.V8-20180723.20150617-S053856-E071129.007384.V06A.HDF5'

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

###Define cross-section length

ind1 = np.where((lon[:,0]>=-100.)) #Where are the DPR lons >= -100
ind2 = np.where((lon[:,0])<=-85.) #Where are lons <= -85
ind3 = np.intersect1d(ind1,ind2) #Both conditions need to be satisfied here


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


    
ku = DPR['NS']['SLV']['zFactorCorrected'][ind3,:,:] #Read in ku-band reflectivity; nscan x nray (554,49,176)
n0 = dsd[ind3,:,:,0] #Read in the number concentration
d0 = dsd[ind3,:,:,1] #Read in the mean drop diameter  #Both have dimensions nscan x nray x nbin (554,49,176)

#Cut all parameters so they are at same ray as above
ku = ku[:,27,:]
n0 = n0[:,27,:]
d0 = d0[:,27,:]

#Take lats and lons along same ray
lons = DPR['NS']['Longitude'][ind3,27]
lats = DPR['NS']['Latitude'][ind3,27]


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
h4 = h3[:,27] #This number MUST match the same ray being used
for i in range(y.shape[1]):
    y[:,i] = h4[i]
    
#%% 
    
#Now we plot!     

#Plot the along track cross-sections
    
plt.figure(figsize=(10,10))

vmax = 60
vmin = 12

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, ku, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
plt.xlabel('N-S Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title('GPM Overpass 6/17 0538 UTC N-S Cross-Section', size = 20)
plt.xlim(570,650)
plt.ylim(0,15)
plt.colorbar(pm, label = 'dBZ')
plt.clim(12,60)


plt.show()



plt.figure(figsize=(10,10))

vmax = 4
vmin = 0

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, ku, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
plt.xlabel('N-S Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 6/17 0538 UTC N-S Cross-Section $D_{0}$', size = 20)
plt.xlim(570,650)
plt.ylim(0,15)
plt.colorbar(pm, label = 'mm')
plt.clim(12,60)


plt.show()



#%%

#Now for the cross-track cross-sections

ku = DPR['NS']['SLV']['zFactorCorrected'][2884,:,:] #ind3 provides list of indices; select median index and hardcode in here for center point
lons = DPR['NS']['Longitude'][2884,:]
lats = DPR['NS']['Latitude'][2884,:]
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
    


    




    
    
    











              









