# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 20:38:56 2020

@author: noahb
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter
import pyart

file = 'omega.nc'
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

level = nc.variables['level'][:]*100
level_mb = level/100

time = nc.variables['time'][:]
timeUnits = nc.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
era['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
era['day'] = np.asarray([d.day for d in era['date']])
era['month'] = np.asarray([d.month for d in era['date']])
era['year'] = np.asarray([d.year for d in era['date']])


omega = nc.variables['w'][:,:,ilat,ilon]
temp = nc.variables['t'][:,:,ilat,ilon]
#Calculate potential temperature from temperature

def potential_temp(temp,pres):
    theta = temp*(100000/pres)**(287/1004)
    return theta


T,Z,I,J = temp.shape
tempsquish = temp.reshape(T,Z,I*J, order = 'F')


theta_squished = np.ones((tempsquish.shape[0],tempsquish.shape[1],tempsquish.shape[2]))*np.nan

for i in range(theta_squished.shape[0]):
    for j in range(theta_squished.shape[1]):
        for k in range(theta_squished.shape[2]):
            theta_squished[i,j,k] = potential_temp(tempsquish[i,j,k], level[j])


#%%
theta = theta_squished.reshape(T,Z,I,J)        

#Define a constant latitude and take cross section along this; varies by longitude (-99.5 -93.5 )
#constant latitude is 30degN

longitude_index = 11

theta_bill = theta[:,:,longitude_index,:]
omega_bill = omega[:,:,longitude_index,:]


#%%
#Now we can plot les data! 
sigma = 1.25
theta_smooth = gaussian_filter(theta_bill, sigma)

time = 5

fig = plt.figure(figsize = [10,10])
ax = plt.axes()
#plt.yscale('log')
plt.ylim(level_mb[24],level_mb[0])
clevs = np.arange(-6,6,.25)
plt.contour(longitude,level_mb,omega_bill[time,:,:],clevs,colors = 'black')
cp = plt.contourf(longitude, level_mb,omega_bill[time,:,:],clevs, cmap = 'bwr')
clevs2 = np.arange(210,370,10)
plt.contour(longitude,level_mb,theta_smooth[time,:,:],colors = 'red',linewidths = 2)
cs = plt.contour(longitude,level_mb,theta_smooth[time,:,:], clevs2,colors = 'red',linewidths = 2)


plt.clabel(cs,inline = 1, fontsize = 10, fmt='%4.0f')
cbar = plt.colorbar(cp,ticks = clevs[::2], orientation = 'horizontal')
cbar.set_label(r'$\mu bs^{-1}$', size = 16)
plt.xlabel('Longitude', size = 18)
plt.ylabel('Pressure (mb)', size = 18)
plt.title('Vertical Velocity and $\Theta$ at Latitude = 38$^o$N 6/20 0009 UTC', size = 20)


plt.show()











