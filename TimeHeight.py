#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 16:25:38 2019

@author: noahbrauer
"""

from Import_Radar import read_file, diff_reflect,reflect_ncdc, cmap_zdr, cmap_kdp, zdrcolor

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob
from matplotlib.colors import LinearSegmentedColormap
import pyart

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


sample_file = 'nexrad_3d_v4_1_20150616T222000Z.nc'
nc_sample = Dataset(sample_file, 'r')
lon = nc_sample.variables['Longitude'][:]-360
lat = nc_sample.variables['Latitude'][:]
altitude = nc_sample.variables['Altitude'][:]

data2 = read_file(sample_file)
z_sample = data2['Z_H']['values'][:]


#Now mesh the grids

lat2,lon2 = np.meshgrid(lat,lon)

files = glob.glob('*.nc')

#loop through each time; select coordinate for Port O'Connor, TX on 6/16; -96.43W (218), 28.49 N (167); index


#El Campo, TX: -96.233, 29.177; lon(180), lat (200)



###Try 194(lat) 237 (lon)




#altitude, lat,lon

z = np.ones((274,29,8,8))*np.nan
zdr = np.ones((274,29,8,8))*np.nan
kdp = np.ones((274,29,8,8))*np.nan
cc = np.ones((274,29,8,8))*np.nan
adp = np.ones((274,29,8,8))*np.nan

#Alright, loop now 

for n,j in enumerate(files):
    data = read_file(j)
    print(j)
    

    
        
    z[n,:,:,:] = data['Z_H']['values'][:,190:198,233:241]
    zdr[n,:,:,:] = data['zdr']['values'][:,190:198,233:241]
    kdp[n,:,:,:] = data['kdp']['values'][:,190:198,233:241]
    cc[n,:,:,:] = data['cc']['values'][:,190:198,233:241]
    adp[n,:,:,:] = data['adp']['values'][:,190:198,233:241]
    
print(np.nanmax(z))    
#%%
    
#Reshape arrays such that dimensions are: (time, altitude, x*y)
    
T,Z,I,J  = z.shape   

z_reshape = z.reshape(T,Z,I*J, order = 'F')
zdr_reshape = zdr.reshape(T,Z,I*J, order = 'F')
kdp_reshape = kdp.reshape(T,Z,I*J, order = 'F')
cc_reshape = cc.reshape(T,Z,I*J, order = 'F')
adp_reshape = adp.reshape(T,Z,I*J, order = 'F')
   

print(np.nanmax(z_reshape))
#Loop through time for each lat-lon coordinate and calculate 5-point radial mean around each elevation scan    


z_mean = np.ones((274,29))*np.nan
zdr_mean = np.ones((274,29))*np.nan
kdp_mean = np.ones((274,29))*np.nan
cc_mean = np.ones((274,29))*np.nan
adp_mean = np.ones((274,29))*np.nan


for i in range(z_mean.shape[0]):
    z_mean[i,:] = np.nanmean(z_reshape[i,:,:], axis = 1)
    zdr_mean[i,:] = np.nanmean(zdr_reshape[i,:,:], axis = 1)
    kdp_mean[i,:] = np.nanmean(kdp_reshape[i,:,:], axis = 1)
    cc_mean[i,:] = np.nanmean(cc_reshape[i,:,:], axis = 1)
    adp_mean[i,:] = np.nanmean(adp_reshape[i,:,:], axis = 1)



#%%
    
#Define new altitude array with evenly spaced levels

altlist = np.array([1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.])    
altitude_interp = np.arange(1,22.5,0.5)    


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 26)
#plt.xlabel('Time', size = 14)
plt.title(r'El Campo $Z_{H}$  6/16', size = 32)  
k=plt.imshow(z_mean.T, cmap='pyart_NWSRef', interpolation = 'hanning', origin = 'lower', extent = [0,274,0,22], vmin = 0, vmax = 75, aspect = 'auto' )
plt.ylim(altlist[0], altlist[27])
plt.yticks(altlist[::6],altlist[::6], fontsize = 16)
plt.colorbar().set_label(label='dBZ',size=26)

#plt.xticks(np.arange(0,288, step = 72),time_cut2930[::72], fontsize = 16)
plt.tight_layout()



fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 26)
#plt.xlabel('Time', size = 14)
plt.title(r'El Campo $Z_{DR}$  6/16', size = 32)  
k=plt.imshow(zdr_mean.T, cmap=zdrcolor, interpolation = 'hanning', origin = 'lower', extent = [0,274,0,22], vmin = -3, vmax = 4, aspect = 'auto' )
plt.ylim(altlist[0], altlist[27])
plt.yticks(altlist[::6],altlist[::6], fontsize = 16)
plt.colorbar().set_label(label='dB',size=26)

#plt.xticks(np.arange(0,288, step = 72),time_cut2930[::72], fontsize = 16)
plt.tight_layout()


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 26)
#plt.xlabel('Time', size = 14)
plt.title(r'El Campo $K_{DP}$  6/16', size = 32)  
k=plt.imshow(kdp_mean.T, cmap=cmap_kdp, interpolation = 'hanning', origin = 'lower', extent = [0,274,0,22], vmin = 0, vmax = 4, aspect = 'auto' )
plt.ylim(altlist[0], altlist[27])
plt.yticks(altlist[::6],altlist[::6], fontsize = 16)
plt.colorbar().set_label(label='deg/km',size=26)

#plt.xticks(np.arange(0,288, step = 72),time_cut2930[::72], fontsize = 16)
plt.tight_layout()




fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 26)
#plt.xlabel('Time', size = 14)
plt.title(r'El Campo $\rho_{hv}$  6/16', size = 32)  
k=plt.imshow(cc_mean.T, cmap=cmap_kdp, interpolation = 'hanning', origin = 'lower', extent = [0,274,0,22], vmin = 0.9, vmax = 1, aspect = 'auto' )
plt.ylim(altlist[0], altlist[27])
plt.yticks(altlist[::6],altlist[::6], fontsize = 16)
plt.colorbar().set_label(label = r'$\rho_{hv}$',size=26)

#plt.xticks(np.arange(0,288, step = 72),time_cut2930[::72], fontsize = 16)
plt.tight_layout()


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 26)
#plt.xlabel('Time', size = 14)
plt.title(r'El Campo Specific Attenuation  6/16', size = 32)  
k=plt.imshow(adp_mean.T, cmap=cmap_kdp, interpolation = 'hanning', origin = 'lower', extent = [0,274,0,22], vmin = 0, vmax = 0.003, aspect = 'auto' )
plt.ylim(altlist[0], altlist[27])
plt.yticks(altlist[::6],altlist[::6], fontsize = 16)
plt.colorbar().set_label(label = r'dB $km^{-1}$', size=26)

#plt.xticks(np.arange(0,288, step = 72),time_cut2930[::72], fontsize = 16)
plt.tight_layout()


   
    
#%%


cmin = 0.; cmax = 75.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  

#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-97,-95.5]); ylim = np.array([28.5,30])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m.contourf(lon2,lat2,z_sample[3,:,:].T,clevs,cmap='pyart_NWSRef',extend='both') 


#m.drawcounties()

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dBZ',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
   
plt.title('GPM Precipitation Rate 9/15/2019 1604 UTC', size = 20)

x2star,y2star = m(-96.05,29.05)
m.plot(x2star,y2star,'ro',markersize=7, color = 'k')


plt.show(block=False)





    