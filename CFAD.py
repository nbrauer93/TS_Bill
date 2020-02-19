#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:58:48 2020

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


sample_file = 'nexrad_3d_v4_1_20150617T222000Z.nc'
nc_sample = Dataset(sample_file, 'r')
lon = nc_sample.variables['Longitude'][:]-360
lat = nc_sample.variables['Latitude'][:]
altitude = nc_sample.variables['Altitude'][:]

data2 = read_file(sample_file)
z_sample = data2['Z_H']['values'][:]


#Now mesh the grids

lat2,lon2 = np.meshgrid(lat,lon)

files = glob.glob('*617*.nc')

#%%

#Grady, OK: 433 (lat), 160 (lon)




#altitude, lat,lon

z = np.ones((288,29,8,8))*np.nan
zdr = np.ones((288,29,8,8))*np.nan
kdp = np.ones((288,29,8,8))*np.nan
cc = np.ones((288,29,8,8))*np.nan
adp = np.ones((288,29,8,8))*np.nan

#Alright, loop now 

for n,j in enumerate(files):
    data = read_file(j)
    print(j)
    

    
        
    z[n,:,:,:] = data['Z_H']['values'][:,429:437,156:164]
    zdr[n,:,:,:] = data['zdr']['values'][:,429:437,156:164]
    kdp[n,:,:,:] = data['kdp']['values'][:,429:437,156:164]
    cc[n,:,:,:] = data['cc']['values'][:,429:437,156:164]
    adp[n,:,:,:] = data['a']['values'][:,429:437,156:164]
    
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


z_mean = np.ones((288,29))*np.nan
zdr_mean = np.ones((288,29))*np.nan
kdp_mean = np.ones((288,29))*np.nan
cc_mean = np.ones((288,29))*np.nan
adp_mean = np.ones((288,29))*np.nan


for i in range(z_mean.shape[0]):
    z_mean[i,:] = np.nanmean(z_reshape[i,:,:], axis = 1)
    zdr_mean[i,:] = np.nanmean(zdr_reshape[i,:,:], axis = 1)
    kdp_mean[i,:] = np.nanmean(kdp_reshape[i,:,:], axis = 1)
    cc_mean[i,:] = np.nanmean(cc_reshape[i,:,:], axis = 1)
    adp_mean[i,:] = np.nanmean(adp_reshape[i,:,:], axis = 1)



time = []    


for j in (files):
    data2 = read_file(j)
    timejunk = data2['Analysis_time']
    time.append(timejunk)

#%%
    
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
    
    
    
#%%    
#Plot the CFAD
z_nonan = np.nan_to_num(z_mean)
zdrnonan = np.nan_to_num(zdr_mean)

print(np.max(z_nonan))
print(z_nonan.shape)

z_hist = np.zeros((10,29))

for i in range(z_hist.shape[1]):
    z_hist[:, i] = create_histogram(z_nonan[:, i], 10, 0, 50)[1]






#%%


#bins_z = np.arange(0,70, step = 5)
#bins_zdr = np.arange(-1,3, step = 0.25)

'''

z_hist = np.zeros((29,11))
z_bins = np.zeros((29,10))
zdr_hist = np.zeros((29,11))
zdr_bins = np.zeros((29,10))

  

for i in range(z_hist.shape[0]):
    z_hist[i,:] = np.histogram(z_nonan[:,i])[1]
    zdr_hist[i,:] = np.histogram(zdrnonan[:,i])[1]
    z_bins[i,:] = np.histogram(z_nonan[:,i])[0]
    zdr_bins[i,:] = np.histogram(zdrnonan[:,i])[0]
    
    
    
z_percentiles = np.zeros((29,10))
zdr_percentiles = np.zeros((29,10))

for i in range(z_percentiles.shape[0]):
    for j in range(z_percentiles.shape[1]):
        
        z_percentiles[i,j] = (z_bins[i,j]/288)*100
        zdr_percentiles[i,j] = (zdr_bins[i,j]/288)*100
'''            




    
    
z_hist[z_hist <= 0] = np.nan
#zdr_hist[zdr_hist <= 0] = np.nan

#z_percentiles[z_percentiles <= 0] = np.nan
#zdr_percentiles[zdr_percentiles <= 0] = np.nan



num_bins = z_hist.shape[0]
num_altitudes = z_hist.shape[1]

fig = plt.figure(figsize = [10,5]) 
plt.xlim(0, num_bins)
plt.pcolormesh(z_hist.T, cmap = 'pyart_NWSRef')
plt.clim(0,100)

x_tick_values = np.linspace(0, num_bins - 1, num=num_bins, dtype=float) + 0.5
x_tick_strings = [None] * num_bins

for k in range(num_bins):
    if k == 0:
        x_tick_strings[k] = '< 5'
    elif k == num_bins - 1:
        x_tick_strings[k] = '>= 45'
    else:
        x_tick_strings[k] = '[{0:d}, {1:d})'.format(k * 5, (k + 1) * 5)

plt.xticks(x_tick_values, x_tick_strings, rotation=90, size=16)

y_tick_values = np.linspace(0, num_altitudes - 1, num=num_altitudes, dtype=float) + 0.5
plt.yticks(y_tick_values[::2], altitude[::2], size=16)

plt.ylabel('Altitude (km)', size = 16)
plt.xlabel('Bin range (dBZ)', size = 16)
plt.title('Grady, OK 6/17 Contoured Frequency by Altitude', size = 18)
plt.colorbar().set_label(label = 'Frequency',size=26)
plt.show()

'''
fig = plt.figure(figsize = [10,5]) 
plt.xlim(0,5)
plt.pcolormesh(zdr_hist.T, cmap = 'pyart_NWSRef')
plt.clim(0,10)
plt.ylabel('Altitude (km)', size = 16)
plt.xlabel('dB', size = 16)
plt.title('Grady, OK 6/17 Contoured Frequency by Altitude', size = 18)
plt.colorbar().set_label(label = 'Frequency',size=26)
plt.show()        
'''   