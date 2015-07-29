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

PATH2SIM='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/GP/BFRIC_0.015/2013_Aug_01_3D/output/subdomain_GP1_0001.nc'
PATH2OBS='/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_F_20130801_78_2_001_sE15.mat'
SAVEPATH='~/karsten/data/'+REGION+'/timeseries/speed/'

"""
The program is translated from Mitchell O'Flaherty's speed-plotting matlab code
and adapted from some of his python workspace functions.

It creates equally spaced vectors onto a spatially-varying speed map. The region
is defined, as are the drifter and FVCOM files.

USAGE: python modelSpeedQuiver.py
"""

def loadnc(datadir):
    """
    Loads a .nc data file.
    """

    # initialize dictionary for the data; store filepath
    data = {}
    data['filepath'] = datadir

    # load data
    try:
        ncid = netcdf.netcdf_file(datadir, 'r', mmap=True)
    except IOError:
        sys.exit('FVCOM file not found!')

    for i in ncid.variables.keys():
        data[i]=ncid.variables[i].data

    if data.has_key('nv'):
        data['nv']=data['nv'].astype(int).T-1
        data['trigrid'] = mplt.Triangulation(data['lon'], data['lat'], data['nv'])
        data['trigridxy'] = mplt.Triangulation(data['x'], data['y'], data['nv'])
    if data.has_key('nbe'):
        data['nbe']=data['nbe'].astype(int).T-1
    if ncid.dimensions.has_key('nele'):
        data['nele'] = ncid.dimensions['nele']
    if ncid.dimensions.has_key('node'):
        data['node'] = ncid.dimensions['node']

    ncid.close()

    # Gets the lon / lat data
    # worth trying to find other lon/lat data within and round directory?
    if (data.has_key('lon') and data.has_key('x')):
        if ((data['lon'].sum()==0).all() or (data['x']==data['lon']),all()):
            data['lon'] = data['x']
            data['lat'] = data['y']

    if data.has_key('nv'):
        data['trigrid'] = mplt.Triangulation(data['lon'], data['lat'], data['nv'])
        data['trigridxy'] = mplt.Triangulation(data['x'], data['y'], data['nv'])

    return data


def ncsortdata(data, uvhset=False, trifinder=False):
    """
    From the nc data provided, common variables are produced.
    """
    nodexy = np.zeros((len(data['x']),2))
    x = data['x']
    y = data['y']
    nv = data['nv']
    lon = data['lon']
    lat = data['lat']

    #make uvnodes by averaging the values of ua/va over the nodes of an element
    nodexy = np.empty((len(lon),2))
    nodexy[:,0] = x
    nodexy[:,1] = y
    uvnode = np.empty((len(nv[:,0]),2))
    uvnode[:,0] = (x[nv[:,0]] + x[nv[:,1]] + x[nv[:,2]]) / 3.0
    uvnode[:,1] = (y[nv[:,0]] + y[nv[:,1]] + y[nv[:,2]]) / 3.0

    nodell = np.empty((len(lon),2))
    nodell[:,0] = lon
    nodell[:,1] = lat
    uvnodell = np.empty((len(nv[:,0]),2))
    uvnodell[:,0] = (lon[nv[:,0]] + lon[nv[:,1]] + lon[nv[:,2]]) / 3.0
    uvnodell[:,1] = (lat[nv[:,0]] + lat[nv[:,1]] + lat[nv[:,2]]) / 3.0

    if uvhset:
        uvh= np.empty((len(nv[:,0]),1))
        uvh[:,0] = (data['h'][nv[:,0]] + data['h'][nv[:,1]] \
                 + data['h'][nv[:,2]]) / 3.0
        data['uvh']=uvh

    data['uvnode'] = uvnode
    data['uvnodell'] = uvnodell
    data['nodell'] = nodell
    data['nodexy'] = nodexy

    if data.has_key('time'):
        data['time']=data['time']+678576

    if data.has_key('trigrid')==False:
        if (data.has_key('nv') and data.has_key('lat') and data.has_key('lon')):
            data['trigrid'] = mplt.Triangulation(data['lon'], data['lat'], \
                data['nv'])

    if data.has_key('trigridxy')==False:
        if (data.has_key('nv') and data.has_key('x') and data.has_key('y')):
            data['trigridxy'] = mplt.Triangulation(data['x'], data['y'], \
                data['nv'])

    if trifinder:
        data['trigrid_finder']=data['trigrid'].get_trifinder()
        data['trigridxy_finder']=data['trigridxy'].get_trifinder()

    return data

def speed_plot(i):
    """
    Plots the time steps for the FVCOM file.
    """

    print i
    f = plt.figure()
    ax = plt.axes([.125,.1,.775,.8])
    triax=ax.tripcolor(data['trigrid'],np.sqrt(data['ua'][i,:]**2 \
            + data['va'][i,:]**2), vmin=CMIN, vmax=CMAX)

    # pretty_plotll

    f.savefig(savepath + LOC + DATE + '_' + GRID + ("%04d" %i) + '.png', dpi=150)
    plt.close(f)



if __name__ == '__main__':

    # identify the fvcom .nc file to load
    filepath = PATH2SIM
    name = LOC + '_d_' + DATE

    data = loadnc(filepath)
    data = ncsortdata(data)

    regions = {'pp' : [-66.225, -66.195, 44.37, 44.41],
               'mp' : [-64.52, -64.3, 45.3, 45.4],
               'gp' : [-66.38, -66.29, 44.213, 44.32],
               'dg' : [-65.79, -65.73, 44.65, 44.7] }

    region = regions(LOC)

    savepath = SAVEPATH
    if not os.path.exists(savepath): os.makedirs(savepath)

    pool = Pool(8)
    pool.map(speed_plot, range(STARTTIME, ENDTIME))
