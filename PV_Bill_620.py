#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:47:46 2020

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter

import pyart

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

file = 'pv.nc'
era = {}

nc = Dataset(file, 'r', unpack = True)

lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:] - 360

xlim = [-93,-85]; ylim = [36,41]
ilat = np.where((lat>ylim[0])&(lat<ylim[1]))[0]
ilon = np.where((lon>xlim[0])&(lon<xlim[1]))[0]

latitude = lat[ilat]
longitude = lon[ilon]


lat2,lon2 = np.meshgrid(latitude,longitude)


#%%
level = nc.variables['level'][:]*100
level_mb = level/100

time = nc.variables['time'][:]
timeUnits = nc.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
era['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
era['day'] = np.asarray([d.day for d in era['date']])
era['month'] = np.asarray([d.month for d in era['date']])
era['year'] = np.asarray([d.year for d in era['date']])


pv = nc.variables['pv'][:,:,ilat,ilon]
temp = nc.variables['t'][:,:,ilat,ilon]




#%%
#Calculate potential temperature from temperature

def potential_temp(temp,pres):
    theta = temp*(100000/pres)**(287/1004)
    return theta


T,Z,I,J = temp.shape
tempsquish = temp.reshape(T,Z,I*J, order = 'F')


theta_squished = np.ones((3,37,589))*np.nan

for i in range(theta_squished.shape[0]):
    for j in range(theta_squished.shape[1]):
        for k in range(theta_squished.shape[2]):
            theta_squished[i,j,k] = potential_temp(tempsquish[i,j,k], level[j])
 

#Now reshape back into oringal form
        
theta = theta_squished.reshape(T,Z,I,J)        
#%%

#Define a constant latitude and take cross section along this; varies by longitude (-99.5 -93.5 )
#constant latitude is 30degN

theta_bill = theta[:,:,11,:]
pv_bill = pv[:,:,11,:]*10**6






sigma = 1.25
theta_smooth = gaussian_filter(theta_bill, sigma)



fig = plt.figure(figsize = [10,10])
ax = plt.axes()
#plt.yscale('log')
plt.ylim(level_mb[36],level_mb[11])
clevs = np.arange(0,10,0.5)
plt.contour(longitude,level_mb,pv_bill[2,:,:],clevs,colors = 'black')
cp = plt.contourf(longitude, level_mb,pv_bill[2,:,:],clevs, cmap = 'pyart_NWSRef')
clevs2 = np.arange(210,370,10)
plt.contour(longitude,level_mb,theta_smooth[2,:,:],colors = 'red',linewidths = 2)
cs = plt.contour(longitude,level_mb,theta_smooth[2,:,:], clevs2,colors = 'red',linewidths = 2)


plt.clabel(cs,inline = 1, fontsize = 10, fmt='%4.0f')
cbar = plt.colorbar(cp,ticks = clevs[::2], orientation = 'horizontal')
cbar.set_label('PVU', size = 16)
plt.xlabel('Longitude', size = 18)
plt.ylabel('Pressure (mb)', size = 18)
plt.title('Potential Vorticity and $\Theta$ at Latitude = 38$^o$N 6/20 06 UTC', size = 20)


plt.show()
