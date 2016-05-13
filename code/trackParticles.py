#! /usr/env/python2.7

"""
usage:
python trackParticles.py [GP/PP/DG] YYYY_Mmm_DD_3D [0.015/0.012/0.009] path/2/pyticle/tracker/output.nc
"""

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
PATH2OUT="/array/home/119865c/karsten/plots/pytrkr/"


def dn2dt(datenum):
    """
    Convert matlab datenum to python datetime.
    input:
        - matlab datenum
    returns:
        - python datetime
    """
    print 'converting matlab datenum to python datetime...'

    return datetime.fromordinal(int(datenum)) + \
           timedelta(days=datenum % 1) - timedelta(days=366)


def createColorMap(model, var, title='', mesh=True, bounds=[], debug=True):
    """
    2D colormap plot of a given variable and mesh. This function is adapted from
    PySeidon's colormap_var, except it is customized to add the plot to an
    existing figure. Holds the plot.

    input:
        - var = gridded variable, 1D numpy array (nele or nnode)
        - title = plot title, string
        - mesh = boolean, True with mesh, False without mesh
        - bounds = list, constricted region subdomain in form of
            [lon.min, lon.max, lat.min, lat.max]
    returns:
        - figure for future plotting
    """

    if debug:
        print '\tplotting grid...'
    # figure if var has nele or nnode dimensions
    if var.shape[0] == model.Grid.nele:
        dim = model.Grid.nele
    elif var.shape[0] == model.Grid.nnode:
        dim = model.Grid.nnode
    else:
        sys.exit('variable has the wrong dimension, shape not equal to grid ' \
                + 'nummber of elements or nodes')

    # bounding box nodes, elements and variables
    lon = model.Grid.lon[:]
    lat = model.Grid.lat[:]
    if debug:
        print '\tcomputing bounding box...'
    if bounds:
        bb = bounds
    else:
        bb = [lon.min(), lon.max(), lat.min(), lat.max()]

    if not hasattr(model.Grid, 'triangleLL'):
        # mesh triangle
        if debug:
            print '\tcomputing triangulation...'
        trinodes = model.Grid.trinodes[:]
        tri = Tri.Triangulation(lon, lat, triangles=trinodes)
    else:
        tri = model.Grid.triangleLL

    # setting limits and levels of colormap
    if debug:
        print '\tcomputing cmin...'
    cmin = var[:].min()
    if debug:
        print '\tcomputing cmax...'
    cmax = var[:].max()
    step = (cmax-cmin) / 50.0

    # depth contours to plot
    levels = np.arange(cmin, (cmax+step), step)

    # define figure window
    if debug:
        print '\tcreating subplot...'

    fig = plt.figure(figsize=(18,10))
    plt.rc('font', size='22')
    ax = fig.add_subplot(111, aspect=(1.0/np.cos(np.mean(lat) * np.pi/180.0)))

    if debug:
        print '\tcomputing colormap...'
    cmap = plt.cm.jet
    f = ax.tripcolor(tri, var[:], vmax=cmax, vmin=cmin, cmap=cmap)

    if mesh:
        plt.triplot(tri, color='white', linewidth=0.5)

    # label and axis parameters
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.patch.set_facecolor('0.5')
    cbar = fig.colorbar(f, ax=ax)
    cbar.set_label('Mean Velocity Norm (m/s)', rotation=-90, labelpad=30)
    scale = 1

    # ticker for coordinate degree axis
    if debug:
        print '\tconfiguring axis...'
    ticks = tic.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
    ax.xaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_formatter(ticks)
    ax.set_xlim([bb[0], bb[1]])
    ax.set_ylim([bb[2], bb[3]])
    ax.grid()
    plt.title(title)

    plt.hold('on')

    if debug:
        print '...colormap passed.'

    return fig


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
