#! /usr/env/python2.7

"""
Creates plots and computes stats for pyticle tracker output.

usage: python run_all_drifters.py [-p -a -s -w] -r N -n N

cmd line args:
    -p :: plot
    -a :: compute stats
    -s :: save plots
    -w :: overwrite current plots
    -b :: select a bottom friction
    -l :: select a location
    -d :: select a simulation date
    -n :: number of particles (int)
    -r :: radius of initialisation (float)
"""


# lib imports
from pyticle_tracker import pyticle
import sys, os
import os.path as osp
import scipy as sp
import numpy as np
import time
import scipy.io as sio
import netCDF4 as nc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import matplotlib.tri as Tri
from matplotlib import rc
from datetime import datetime, timedelta
from pyseidon import *
import argparse as arp

LOC = ['DG', 'GP', 'PP']
BFRIC = '0.015'
SIM = ['2014_Aug_12_3D', '2013_Aug_08_3D', '2013_Aug_01_3D', '2012_Jun_05_3D', \
       '2012_Jul_05_3D', '2013_Oct_10_3D', '2013_Nov_05_3D', '2013_Nov_06_3D']

PATH2OUT = '/EcoII/acadia_uni/projects/drifters/plots/pytrkr/'
PATH2SIM = '/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/' + \
                        'drifter_runs/BFRIC_'
OUTPATH = '/EcoII/acadia_uni/projects/drifters/pyticle_tracker/'
PATH2INIT = '/array/home/119865c/karsten/drifters/start_info/'
PATH2OBS = '/EcoII/acadia_uni/workspace/observed/'


def parseArgs():
    """
    Parses command line args.
    """

    parser = arp.ArgumentParser(prog='run_all_drifters.py')

    parser.add_argument('-p', action='store_true', help='generate plots.')
    parser.add_argument('-s', action='store_true', help='save plots.')
    parser.add_argument('-a', action='store_true', help='run analysis.')
    parser.add_argument('-w', action='store_true', help='overwrite plots.')
    parser.add_argument('-b', nargs=1, choices=('0.009','0.012','0.015'), \
            help='select a bottom friction.', default='0.015', type=str)
    parser.add_argument('-l', nargs='*', choices=LOC, \
            help='select a location tag.')
    parser.add_argument('-d', nargs='*', choices=SIM,\
            metavar='YYYY_Mmm_DD_3D', help='select a simulation date.')
    multiple = parser.add_argument_group('multiple pyticle options')
    multiple.add_argument('-r', nargs=1, type=float, metavar='#', \
            help='define an initial radius in degrees.')
    multiple.add_argument('-n', type=int, metavar='#', nargs=1, \
            help='define a number of particles.')
    # parser.add_argument('-sd', nargs=1, metavar='#', type=float,\
    #         help='adds random noise as std. dev. to pyticle velocities.')
    parser._optionals.title = 'optional flag arguments'

    args = parser.parse_args()

    return args


def dn2dt(datenum):
    """
    Convert matlab datenum to python datetime.
    input:
        - matlab datenum
    return:
        - python datetime
    """

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

    # set cmd line options
    print 'parsing command line options...'
    args = parseArgs()

    if args.l:
        LOC = args.l
    if args.d:
        SIM = args.d

    if args.b:
        bfric = args.b
    else:
        bfric = BFRIC

    for loc in LOC:
        for sim in SIM:

            filename = PATH2SIM + bfric + '/' + loc + '/' + sim + \
                        '/output/subdomain_' + loc + '1_0001.nc'
            start_info = PATH2INIT + 'pyticle_info_' + loc + '_' + sim + '.txt'
            path2drift = PATH2OBS + loc + '/Drifter/'
            outpath = OUTPATH

            print 'looking for ncfile...'
            if not osp.exists(filename):
                continue

            # define output path
            outpath = outpath + loc + '_' + sim

            if args.n:
                outpath = outpath + '_n{}'.format(args.n[0])
            if args.r:
                outpath = outpath + '_r{}'.format(args.r[0])

            outpath += '/'

            if not osp.exists(outpath):
                os.makedirs(outpath)

            # get starting locations and timesteps of drifters
            indata=np.genfromtxt(start_info, dtype=None)
            print str(indata.shape[0]) + ' drifters...'

            # set centre of location
            if loc == 'GP':
                centre = [-66.33906, 44.26898]
            elif loc == 'DG':
                centre = [-65.76000, 44.67751]
            elif loc == 'PP':
                centre = [-65.206924, 44.389368]
            else:
                sys.exit('location tag not recognized.')

            if not osp.exists(filename) or not osp.isfile(filename):
                sys.exit('simulation path not found / valid.')

            # open the FVCOM file
            print 'opening fvcom file...'
            ncfile = FVCOM(filename, debug=False)
            print 'calculating ebb/flood split at centre of location...'
            fI, eI, _, _ = ncfile.Util2D.ebb_flood_split_at_point(centre[0], \
                            centre[1])
            print 'calculating model velocity norm...'
            ncfile.Util3D.velo_norm()

            # set seaborn sns font
            sns.set(font="serif")
            # activate latex text rendering
            # rc('text', usetex=True)

            if bfric not in ['0.015', '0.012', '0.009']:
                sys.exit('bottom friction tag not valid.')

            for row in indata:
                drifter = path2drift + row[0]
                inloc = [row[1], row[2]]
                inlocs = inloc
                savedir = outpath + row[0][:-4] + '_output.nc'

                # added radius and number of particles
                if args.n:
                    num = args.n[0]
                    inlocs = np.tile(inloc, (num, 1))
                    print 'starting with {} particles...'.format(num)
                else:
                    num = 1

                if args.r:
                    # 1 deg lat = 110575 m
                    # 1 deg lon = 111303 m
                    inlocs = np.vstack((np.random.uniform( \
                            inloc[0]-args.r[0]/111303.0, \
                            inloc[0]+args.r[0]/111303.0, num), \
                            np.random.uniform(inloc[1]-args.r[0]/110575.0, \
                                              inloc[1]+args.r[0]/110575.0, \
                                              num))).T
                    print 'randomizing starting locations...'

                # if the run exists, skip it
                if not osp.exists(savedir):
                    # set options of drifters
                    # note: interpolation ratio is how many timesteps per
                    # model timestep to linearly interpolate nc data
                    # output ratio is how often to output particle potision
                    options={}
                    options['starttime']=row[3]
                    options['endtime']=row[4]
                    options['interpolationratio']=60
                    options['outputratio']=2
                    options['ncformat']='NETCDF4_CLASSIC'
                    options['useLL']=True
                    options['layer']=0
                    options['gridDim']='2D'
                    options['projstr']='lcc +lon_0=-64.55880 +lat_0=41.78504 '+ \
                                   '+lat_1=39.69152 +lat_2=43.87856'

                    # run mitchell's particle tracker
                    start = time.clock()
                    mypy=pyticle(filename, inlocs, savedir, options=options)
                    mypy.run()
                    print('run in: %f' % (time.clock() - start))

                # open pytkl structure
                pytkl = nc.Dataset(savedir, 'r', format='NETCDF4_CLASSIC')

                print 'creating time window...'
                # this is in julian time
                tModel = ncfile.Variables.julianTime
                tPytkl = pytkl.variables['time'][:]

                win1 = (np.abs(tModel - tPytkl.min())).argmin()
                win2 = (np.abs(tModel - tPytkl.max())).argmin()

                tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)

                print 'opening drifter file...'
                drift = Drifter(path2drift+row[0], debug=False)

                # do things based on command line args
                if args.p:
                    # check whether or not to overwrite
                    savename = row[0][:-4] + '_plot.png'
                    saveplot = PATH2OUT + loc + '_' + sim + '/'
                    if args.s:
                        save = True
                        if args.w:
                            ow = True
                        else:
                            ow = False
                    else:
                        save = False
                        ow = False

                    if save and not ow:
                        if osp.exists(saveplot+savename):
                            continue

                    print 'preparing to create color map...'
                    fig = createColorMap(ncfile, tideNorm[0,:], mesh=False, \
                            title = 'Pyticle Track for ' + row[0])

                    # add scatter plots of data
                    lon = pytkl.variables['lon'][:]
                    lat = pytkl.variables['lat'][:]
                    print 'adding scatter plot...'
                    pytkl_path = plt.scatter(lon, lat, c='m', lw=0, alpha=0.6, \
                        marker="^", s=25)

                    lonD = drift.Variables.lon
                    latD = drift.Variables.lat
                    drift_path = plt.scatter(lonD, latD, c='k', lw=0, s=25)

                    plt.legend([drift_path, pytkl_path],
                           ['Drifter', 'Model'])

                    if save:
                        print 'creating save directory...'
                        if not osp.exists(saveplot):
                            os.makedirs(saveplot)

                        plt.savefig(saveplot + savename)

                    else:
                        plt.show()

                    plt.close(fig)

                if args.a:
                    pass
