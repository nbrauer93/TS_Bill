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
from mpl_toolkits.basemap import shiftgrid
import pyart
import seaborn as sns


with gzip.open('GRtoDPR.KPAH.150620.7430.V05A.DPR.MS.1_21.nc.gz') as gz:
    with netCDF4.Dataset('file', mode = 'r', memory=gz.read()) as nc:
        print(nc.variables)
        latitude = nc.variables['latitude'][:].T
        longitude = nc.variables['longitude'][:].T
        z = nc.variables['ZFactorCorrected'][:]
        #zdr = nc.variables['']
        elev_angle = nc.variables['timeSweepStart'][:]



elev_angle = 0
z_lowest = z[elev_angle,:]
#z_lowest[z_lowest<=0] = np.nan





#%%

#from matplotlib.colors import ListedColormap

colormap = ['white','dodgerblue', 'deepskyblue', 'lawngreen', 'lightgreen', 'green', 'gold', 'darkorange', 'red', 'firebrick']

z_color = np.empty(667,dtype = 'str')
z_color = []

#for i in range(len(z_color)):
for i in range(len(z_lowest)):

    if z_lowest[i]<=5:
        #z_color[i] = colormap[0]
        z_color.append(colormap[0])
    
    elif z_lowest[i]>5 and z_lowest[i]<=10:
        #z_color[i] = colormap[1]
        z_color.append(colormap[1])
    
    elif z_lowest[i]>10 and z_lowest[i]<=15:
        #z_color[i] = colormap[2]
        z_color.append(colormap[2])
    
    elif z_lowest[i]>15 and z_lowest[i]<=20:
        #z_color[i] = colormap[3]
        z_color.append(colormap[3])
    
    elif z_lowest[i]>20 and z_lowest[i]<=25:
        #z_color[i] = colormap[4]
        z_color.append(colormap[4])
    
    elif z_lowest[i]>25 and z_lowest[i]<=30:
        #z_color[i] = colormap[5]
        z_color.append(colormap[5])
        
    elif z_lowest[i]>30 and z_lowest[i]<=35:
        #z_color[i] = colormap[6]
        z_color.append(colormap[6])
    
    elif z_lowest[i]>35 and z_lowest[i]<=40:
        #z_color[i] = colormap[7]
        z_color.append(colormap[7])
        
    elif z_lowest[i]>40 and z_lowest[i]<=45:
        #z_color[i] = colormap[8]
        z_color.append(colormap[8])
        
    elif z_lowest[i] == -100:
        #z_color[i] = colormap[0]
        z_color.append(colormap[0])
        
        
from matplotlib.colors import ListedColormap
cmap_z = ListedColormap(colormap)        
    



#%%
import matplotlib

elev_angle = 0
#Setup plotting 
cmin = 0; cmax = 50; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_z,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(z_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = z_color, vmin = 0, vmax = 50, cmap = cmap, edgecolors = 'k')
c2 = plt.scatter(longitude[:,elev_angle], latitude[:,elev_angle], facecolors = 'none' )
plt.xlabel('Longitude', size = 24)
plt.ylabel('Latitude', size = 24)
plt.title(r'KPAH $Z_{H}$ 6/20 0743 UTC ',name='Calibri',size=26)

colour_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

colour_bar_object.set_label('dBZ', size = 24)

plt.show()
plt.close(figure_object)


