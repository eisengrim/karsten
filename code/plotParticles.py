#! /usr/env/python2.7

"""
usage:
python plotParticles.py [GP/PP/DG] YYYY_Mmm_DD_3D [0.015/0.012/0.009] path/2/pyticle/tracker/output.nc
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

PATH2SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
PATH2OBS="/EcoII/acadia_uni/workspace/observed/"
PATH2OUT="/EcoII/acadia_uni/projects/drifters/plots/pytrkr/"

# local imports
from drifterUtils import *
from createColorMap import createColorMap


if __name__ == '__main__':

    if len(sys.argv) < 5:
        sys.exit('not enough command line args. see usage.')

    try:
        filename = str(sys.argv[4])
        bfric = str(sys.argv[3])
        simdate = str(sys.argv[2])
        loc = str(sys.argv[1])
    except:
        sys.exit('error parsing command line argument.')

    if not osp.exists(filename):
        sys.exit('nc file not found / invalid.')

    pytkl = nc.Dataset(filename, 'r', format='NETCDF4_CLASSIC')

    sns.set(font="serif")

    if loc == 'GP':
        centre = [-66.33906, 44.26898]
    elif loc == 'DG':
        centre = [-65.76000, 44.67751]
    elif loc == 'PP':
        centre = [-65.206924, 44.389368]
    else:
        sys.exit('location tag not recognized.')

    if bfric not in ['0.015', '0.012', '0.009']:
        sys.exit('bottom friction tag not valid.')

    simpath = PATH2SIM + 'BFRIC_' + bfric + '/' + loc + '/' + simdate + \
            '/output/subdomain_' + loc + '1_0001.nc'
    obspath = PATH2OBS + loc + '/Drifter/'

    if not osp.exists(simpath) or not osp.isfile(simpath):
        sys.exit('simulation path not found / valid.')

    print 'opening fvcom file...'
    ncfile = FVCOM(simpath, debug=False)
    print 'calculating ebb/flood split at centre of location...'
    fI, eI, _, _ = ncfile.Util2D.ebb_flood_split_at_point(centre[0], centre[1])
    print 'calculating model velocity norm...'
    ncfile.Util3D.velo_norm()

    # for now, drifter = false. add this later

    print 'creating time window...'
    # this is in julian time
    tModel = ncfile.Variables.julianTime
    tPytkl = pytkl.variables['time'][:]

    win1 = (np.abs(tModel - tPytkl.min())).argmin()
    win2 = (np.abs(tModel - tPytkl.max())).argmin()

    tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)

    print 'preparing to create color map...'

    fig = createColorMap(ncfile, tideNorm[0,:], mesh=False, \
            title = 'Pyticle Track for ' +loc+'_'+ simdate)

    lon = pytkl.variables['lon'][:]
    lat = pytkl.variables['lat'][:]
    print 'adding scatter plot...'
    plt.scatter(lon, lat, c='k')

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
