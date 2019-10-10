#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:03:49 2019

@author: noahbrauer
"""

import h5py 
import numpy as np
import matplotlib.pyplot as plt

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
import pyart


file = '2A.GPM.DPR.V8-20180723.20150617-S053856-E071129.007384.V06A.HDF5'


DPR = h5py.File(file, 'r')

lat2 = DPR['NS']['Latitude'][:,:]    
lon2 = DPR['NS']['Longitude'][:,:]
z = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:]  #nray, nscan
precip_rate = DPR['NS']['SLV']['precipRateNearSurface'][:]
alt = DPR['NS']['PRE']['elevation'][:] #in meters

#%%
#Choose our lat-lon of interest

ind1 = np.where((lon2[:,0]>=-100))
ind2 = np.where((lon2[:,0])<=-85)
ind3 = np.intersect1d(ind1,ind2)


#Make coordinates into 2 dimensions


x2 = 2.*17 #48 degrees
re = 6378 #radius of earth in km

theta = -1 *(x2/2.) + (x2/48.)*np.arange(0,49) #Get equal degree increments

theta2 = np.zeros(theta.shape[0]+1)
theta = theta - 0.70833333/2.
theta2[:-1] = theta
theta2[-1] = theta[-1] + 0.70833333

theta = theta2 * (np.pi/180.) #convert to radians

prh = np.zeros([177,50]) #set up matrix 
for i in np.arange(0,177): #loop over num range gates
    for j in np.arange(0,50): #loop over scans 
        a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #407 km is the orbit height 
        prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a) #more geometry 
        
h2 = prh

h3 = np.zeros([h2.shape[0],h2.shape[1]])
for i in np.arange(0,h3.shape[1]):
    h3[:,i] = h2[::-1,i] #reverse order so as index go up the height goes up
    
    
from pyproj import Proj

#Extract the vertical profile of reflectivity

ku = DPR['NS']['SLV']['zFactorCorrected'][ind3,:,:]  #176 levels

ku = ku[:,27,:]
lons = DPR['NS']['Longitude'][ind3,27]
lats = DPR['NS']['Latitude'][ind3,27]

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


y = np.zeros([ku.shape[0],ku.shape[1]])

#make sure this number matches the one above 
h4 = h3[:,27]
    
for i in np.arange(y.shape[1]):
    y[:,i] = h4[i]


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
plt.title('GPM Overpass 6/17 1352 UTC N-S Cross-Section', size = 20)
plt.xlim(570,650)
plt.ylim(0,15)
plt.colorbar(pm, label = 'dBZ')
plt.clim(12,60)


plt.show()




#%%

