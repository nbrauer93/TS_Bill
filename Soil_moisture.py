#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 06:34:42 2019

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import netCDF4 as nc
import numpy as np
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob



import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


files = nc.MFDataset('*.nc')

#Import l'attributes from le dataset

soil_moisture = files.variables['soill'][:]
lat = files.variables['lat'][:]
lon = files.variables['lon'][:]
level = files.variables['level'][:]


#Define a dictionary to extract temporal periods

narr = {}

time = files.variables['time'][:]
timeUnits = files.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])

#This is the indexing for TS Bill June 9-16, 2015
bill_index = np.where((narr['day']<17)&(narr['day']>8)& (narr['year']==2015))[0]
#And the indexing for the 1979-2015 climatology for June 9-16
climo_index = np.where((narr['day']<17)&(narr['day']>8))[0]

#Soil moisture for the Bill period (100 cm depth)

soil_bill = soil_moisture[bill_index,3,:,:]

#Soil moisture for climo period (100 cm depth)

soil_climo = soil_moisture[climo_index,3,:,:]


###Okay, now compute mean and standard deviation over entire time period...will do percentiles later for comparison

soil_mean = np.nanmean(soil_climo, axis = 0)
soil_std = np.nanstd(soil_climo, axis = 0)


#Now average mean soil moisture during Bill period (Average June 9-16, 2015 soil moisture)

bill_mean = np.nanmean(soil_bill, axis = 0)

#Now compute standardized anomalies during TS Bill period

bill_ano = (bill_mean - soil_mean)/soil_std

#%%
#Now we plot les data

cmin = -4.1; cmax = 4.1; cint = 0.2; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(10,10))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-104,-92]); ylim = np.array([24,38])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon,lat,bill_ano,clevs,cmap='RdBu',extend='both') 
m.drawcounties()
#m.pcolor(lon, lat, significant30, hatch='.', alpha = 0.)
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel(r'[$\sigma$]',size=32)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(26)

plt.title('100 cm Soil Moisture Anomaly (June 9-16, 2015)',name='Calibri',size=34)
plt.show(block=False) 












