#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 12:13:12 2019

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

#Import file

file = '2A.GPM.DPR.V8-20180723.20150617-S053856-E071129.007384.V06A.HDF5'


DPR = h5py.File(file, 'r')

lat2 = DPR['NS']['Latitude'][:,:]    
lon2 = DPR['NS']['Longitude'][:,:]
z = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:]  #nray, nscan
precip_rate = DPR['NS']['SLV']['precipRateNearSurface'][:]
alt = DPR['NS']['PRE']['elevation'][:] #in meters

#Extract DSD

dsd = DPR['NS']['SLV']['paramDSD'][:]  #nDSD x nbin x nray x nscan

#%%
#Choose our lat-lon of interest

ind1 = np.where((lon2[:,0]>=-100))
ind2 = np.where((lon2[:,0])<=-85)
ind3 = np.intersect1d(ind1,ind2)


#Make coordinates into 2D


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
        prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a) 
        
h2 = prh

h3 = np.zeros([h2.shape[0],h2.shape[1]])
for i in np.arange(0,h3.shape[1]):
    h3[:,i] = h2[::-1,i] #Make vertical indexing more intuitive...ground up
    
    
    