#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import matplotlib as mpl
import matplotlib.tri as mplt
import matplotlib.pyplot as plt

import sys, os
import numpy as np
import scipy as sp
import scipy.io as sio
from scipy.io import netcdf
import glob
from multiprocessing import Pool

LOC='gp'
DATE='aug01_13'
GRID='dngridCSR'
TYPE='3d'
CMIN=0
CMAX=5
STARTTIME=0
ENDTIME=100

PATH2SIM='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/GP/2013_Aug_01_3D/output/subdomain_GP1_0001.nc'
PATH2OBS='/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_F_20130801_78_2_001_sE15.mat'
SAVEPATH='~/karsten/data/'+REGION+'/timeseries/speed/'

"""
The program is translated from Mitchell O'Flaherty's speed-plotting matlab code
and adapted from some of his python workspace functions.

It creates equally spaced vectors onto a spatially-varying speed map. The region
is defined, as are the drifter and FVCOM files.

USAGE: python modelSpeedQuiver.py
"""



