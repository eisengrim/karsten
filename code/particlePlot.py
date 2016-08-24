#! /usr/env/python2.7

"""
usage:
python particlePlot.py path/2/fvcom/output.nc path/2/pyticle/tracker/output.nc
"""

# library imports
import sys, os
import os.path as osp
import scipy as sp
import numpy as np
import netCDF4 as nc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import matplotlib.tri as Tri
from datetime import datetime, timedelta
from pyseidon import *

# local imports
from drifterUtils import *
from createColorMap import createColorMap

path2pytrkr = '/EcoII/acadia_uni/projects/drifters/pyticle_tracker/' + \
        'HH_20111022-29/HH_2011_Oct_22-29_AF_3D_output.nc'
path2fvcom = '/EcoII/acadia_uni/projects/acadia_force_numerical_model_r20/' + \
        '2011-10-22_2011-10-29/output/acadia_force_3d_0001.nc'
loc = 'HH'
sim = '2011_Oct_22-29_AF_3D'
bfric = '0.015'

if __name__ == '__main__':

    if not osp.exists(path2pytrkr):
        sys.exit('nc file not found / invalid.')

    pytkl = nc.Dataset(path2pytrkr, 'r', format='NETCDF4_CLASSIC')

    sns.set(font="serif")

    if loc == 'GP':
        centre = [-66.33906, 44.26898]
    elif loc == 'DG':
        centre = [-65.76000, 44.67751]
    elif loc == 'PP':
        centre = [-65.206924, 44.389368]
    elif loc == 'HH':
        centre = [-64.625361, 45.206072]
    else:
        sys.exit('location tag not recognized.')

    if bfric not in ['0.015', '0.012', '0.009']:
        sys.exit('bottom friction tag not valid.')

    if not osp.exists(path2fvcom) or not osp.isfile(path2fvcom):
        sys.exit('simulation path not found / valid.')

    print 'opening fvcom file...'
    ncfile = FVCOM(path2fvcom, debug=False)
    # print 'calculating ebb/flood split at centre of location...'
    # fI,eI,_,_ = ncfile.Util2D.ebb_flood_split_at_point(centre[0],centre[1])
    # print 'calculating model velocity norm...'
    # ncfile.Util3D.velo_norm()

    # print 'creating time window...'
    ## this is in julian time
    # tModel = ncfile.Variables.julianTime
    # tPytkl = pytkl.variables['time'][:]
    # win1 = (np.abs(tModel - tPytkl.min())).argmin()
    # win2 = (np.abs(tModel - tPytkl.max())).argmin()

    # tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
    bathym = ncfile.Grid.h.clip(max=150)

    print 'preparing to create color map...'
    fig = createColorMap(ncfile, bathym, mesh=False, \
            title = 'Particle Track for ('+loc+','+sim+')',
            label='Depth (m)')
    ax = fig.gca()
    lon = pytkl.variables['lon'][:]
    lat = pytkl.variables['lat'][:]# + 0.055
    print 'adding scatter plot...'
    plt.scatter(lon, lat, c='k')
    ax.set_xlim(pytkl.variables['lon'][:].min()-0.01, \
            pytkl.variables['lon'][:].max()+0.01)
    ax.set_ylim(pytkl.variables['lat'][:].min()-0.05, \
            pytkl.variables['lat'][:].max()+0.05)
    plt.show()


#    savename = raw_input('enter a save name: ')
#    if not savename:

#        sys.exit()
#
#    print 'creating save directory...'
#    savedir = PATH2OUT + loc + '_' + simdate + '/'
#    if not osp.exists(savedir):
#        os.makedirs(savedir)
#
#    plt.savefig(savedir + savename)
#
