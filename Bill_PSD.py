# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:11:57 2019

@author: noahb
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


file = '2A.GPM.DPR.V8-20180723.20150617-S053856-E071129.007384.V06A.HDF5'

DPR = h5py.File(file, 'r')

#7934 bins, 49 scanning levels, 176 rays

lat2 = DPR['NS']['Latitude'][:,:]
lon2 = DPR['NS']['Longitude'][:,:]
z_u = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:] #load in ku data
z_a = DPR['MS']['SLV']['zFactorCorrectedNearSurface'][:] #load in ka data; nray x nscan
dsd = DPR['NS']['SLV']['paramDSD'][:] #nbin x nray x nscan x DSDmoment



#Assign longitudinal domain to plot data
'''
ind1 = np.where(lon2[:,0] >= -100)
ind2 = np.where(lon2[:,0] <= -90)
ind3 = np.intersect1d(ind1,ind2)
'''
#Define an area to average DSD moments over; calculate mean over this domain

#lat: 30-31 deg N
#lon: -96 to -95 W

lon_ind1 = np.where(lon2[:,0]>-96)
lon_ind2 = np.where(lon2[:,0]>-95)
lon_ind3 = np.intersect1d(lon_ind1, lon_ind2)

lat_ind1 = np.where(lat2[:,0]>30)
lat_ind2 = np.where(lat2[:,0]>31)
lat_ind3 = np.intersect1d(lat_ind1, lat_ind2)



#Assign domain to DSD data

dsd_domain = dsd[lon_ind3,:,:,:]

#extract number concentration (0th moment), and mean drop diameter (1st order moment)
#Now extract Nw and D_o and assign to latitudinal domain

psd_0 = dsd_domain[lat_ind3,:,:,0]
psd_1 = dsd_domain[lat_ind3,:,:,1]

#Now psd at 1 and 2 km

nw_1km = psd_0[:,4,:]
nw_2km = psd_0[:,8,:]

d0_1km = psd_1[:,4,:]
d0_2km = psd_1[:,8,:]

###Now extract bad data

nw_1km[nw_1km == -9999.9] = np.nan
nw_2km[nw_2km == -9999.9] = np.nan
d0_1km[d0_1km == -9999.9] = np.nan
d0_2km[nw_2km == -9999.9] = np.nan


#Now compute spatial mean over domain

nw_1km_mean = np.nanmean(nw_1km, axis = 1)
nw_2km_mean = np.nanmean(nw_2km, axis = 1)








        
        
        




















