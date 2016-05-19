from __future__ import division,print_function
import matplotlib as mpl
import scipy as sp
from datatools import *
from gridtools import *
from misctools import *
from plottools import *
from projtools import *
import interptools as ipt
import matplotlib.tri as mplt
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
import os as os
import sys
np.set_printoptions(precision=8,suppress=True,threshold=np.nan)
from scipy.interpolate import interp1d
from matplotlib import patches as pp
from osgeo import osr, gdal
from matplotlib.colors import LinearSegmentedColormap
import collections

# Define names and types of data
name='vh_high_3d_profile'
grid='vh_high'
datatype='2d'

### load the .nc file #####
data = loadnc('runs/'+grid+'/'+name+'/output/',singlename=grid + '_0001.nc')
print('done load')
data = ncdatasort(data)
print('done sort')


idx={}
for node in range(data['node']):
    print(node)
    idx[str(node)]=np.argwhere(data['nv']==node)[:,0]
    

mu=np.empty((len(data['time']),data['node']))
mv=np.empty((len(data['time']),data['node']))
for node in range(data['node']):
    print(node)
    mu[:,node]=np.mean(data['ua'][:,idx[str(node)]],axis=1)
    mv[:,node]=np.mean(data['va'][:,idx[str(node)]],axis=1)
