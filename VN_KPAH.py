#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 18:37:21 2020

@author: noahbrauer
"""

import netCDF4
import gzip
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


with gzip.open('GRtoDPR.KPAH.150620.7430.V05A.DPR.MS.1_21.nc.gz') as gz:
    with netCDF4.Dataset('file', mode = 'r', memory=gz.read()) as nc:
        print(nc.variables)
        latitude = nc.variables['latitude'][:]
        longitude = nc.variables['longitude'][:]
        z = nc.variables['ZFactorCorrected'][:]




lat2,lon2 = np.meshgrid(latitude,longitude)  



#Setup plotting 


cmin = 0; cmax = 60; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(10,10))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-91,-88]); ylim = np.array([36,39])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(longitude,latitude,z,clevs,cmap='pyart_NWSRef',extend='both') 
m.drawcounties()
#m.pcolor(lon, lat, significant30, hatch='.', alpha = 0.)
'''
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dBZ',size=32)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(26)
'''
plt.title(r'KPAH $Z_{H}$ 6/20 0743 UTC ',name='Calibri',size=28)
plt.show(block=False) 
  
