#! /usr/env/python2.7

"""
Creates plots and computes stats for pyticle tracker output.

usage: python trackDrifters.py [-p -a -s -w] -r N -n N

cmd line args:
    -p :: plot
    -a :: compute stats
    -s :: save plots
    -w :: overwrite current model and plots
    -b :: select a bottom friction
    -l :: select a location
    -d :: select a simulation date
    -n :: number of particles (int)
    -r :: radius of initialisation (float)
    -g :: additive white Gaussian noise
    -t :: change the starting time step
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
from sklearn.preprocessing import normalize

# local imports
from createColorMap import createColorMap
from drifterUtils import *
from drifterPlotUtils import plotTimeSeries

LOC = ['DG', 'GP', 'PP']
BFRIC = '0.015'
SIM = ['2014_Aug_12_3D', '2013_Aug_08_3D', '2013_Aug_01_3D', '2012_Jun_05_3D', \
       '2012_Jul_05_3D', '2013_Oct_10_3D', '2013_Nov_05_3D', '2013_Nov_06_3D']

PATH2PLOT = '/EcoII/acadia_uni/projects/drifters/plots/pytrkr/'
PATH2SIM = '/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/' + \
                        'drifter_runs/BFRIC_'
PATH2OUT = '/EcoII/acadia_uni/projects/drifters/pyticle_tracker/'
PATH2INIT = '/array/home/119865c/karsten/drifters/start_info/'
PATH2OBS = '/EcoII/acadia_uni/workspace/observed/'


def parseArgs():
    """
    Parses command line args.
    """

    parser = arp.ArgumentParser(prog='trackDrifters.py')

    parser.add_argument('-p', nargs='*', metavar='#', choices=[1,2,3], \
            default=0, type=int, help='generate plots.\n\t1 :: basic trajectory' \
            + '\n\t2 :: lon/lat discrepancy' + '\n\t3 :: probability density')
    parser.add_argument('-s', action='store_true', help='save plots.')
    parser.add_argument('-a', action='store_true', help='run analysis.')
    parser.add_argument('-w', action='store_true', help='overwrite data.')
    parser.add_argument('-b', nargs=1, choices=('0.009','0.012','0.015'), \
            help='select a bottom friction.', default='0.015', type=str)
    parser.add_argument('-l', nargs='*', choices=LOC, \
            help='select a location tag.')
    parser.add_argument('-d', nargs='*', choices=SIM,\
            metavar='YYYY_Mmm_DD_3D', help='select a simulation date.')
    multiple = parser.add_argument_group('pyticle options')
    multiple.add_argument('-g', action="store_true", \
            help='additive white Gaussian noise.')
    multiple.add_argument('-t', action='store_true', \
            help='alter the initial timestep.')
    multiple.add_argument('-r', nargs=1, type=float, metavar='#', \
            help='define an initial radius in degrees.')
    multiple.add_argument('-n', type=int, metavar='#', nargs=1, \
            help='define a number of particles.')
    # parser.add_argument('-sd', nargs=1, metavar='#', type=float,\
    #         help='adds random noise as std. dev. to pyticle velocities.')
    parser._optionals.title = 'optional flag arguments'

    args = parser.parse_args()

    return args


def plotTracks(row, pytkl, drift, ncfile, tideNorm, sim, loc, args, saveplot):
    """
    Plots the tracks of the two objects.
    """
    # check whether or not to overwrite
    savename = row[0][:-4] + '_plot.png'

    if args.s and not args.w:
        if osp.exists(saveplot + savename):
             return

    print 'preparing to create color map...'
    fig = createColorMap(ncfile, tideNorm[0,:], mesh=False, \
                        title = 'Pyticle Track for ' + row[0])

    # add scatter plots of data
    lon = pytkl.variables['lon'][:]
    lat = pytkl.variables['lat'][:]
    print 'adding scatter plot...'
    pytkl_path = plt.scatter(lon, lat, c='m', lw=0, alpha=0.6, marker="^", s=10)

    lonD = drift.Variables.lon
    latD = drift.Variables.lat
    drift_path = plt.scatter(lonD, latD, c='k', lw=0, s=25)

    plt.legend([drift_path, pytkl_path],['Drifter', 'Model'])

    if args.s:
        print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def analyseTracks(pytkl, drift, ncfile, sim, loc, row, args):
    """
    Runs a bunch of statistics on the two tracks.
    """
    pass


def geographicError(pytkl, drift, ncfile, sim, loc, row, args, saveplot):
    """
    Calculates errors in lon/lat coordinates graphically.
    """
    # check whether or not to overwrite
    savename = row[0][:-4] + '_lonlatErr.png'

    if args.s and not args.w:
        if osp.exists(saveplot + savename):
             return

    print 'calculating lon/lat error...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat

    print str(lon.shape[1]) + ' pyticles...'
    print 'no. pts pytkl: ' + str(len(lon))
    print 'no. pts drift: ' + str(len(lonD))

    uD = drift.Variables.u
    vD = drift.Variables.v
    u = np.array(pytkl.variables['u'][:])
    v = np.array(pytkl.variables['v'][:])

    dtime = drift.Variables.matlabTime
    ptime = mjd2num(pytkl.variables['time'][:])

    dtime = [dn2dt(x) for x in dtime]
    ptime = [dn2dt(x) for x in ptime]

    fig = plt.figure()
    fig.suptitle('Positional Timeseries :: Model vs. Drifter\n' + \
                 'BFRIC={} | Filename={}'.format(str(bfric), row[0][:-4]),
                 fontsize=14)

    print 'creating figures...'
    result = plotTimeSeries(fig, (ptime, dtime), (lon, lonD), loc, \
            label=['Model Drifter','Observed Drifter'], where=121, \
            title='Longitudal Comparison', \
            axis_label='Longitude $^\circ$', styles=['#1DA742','#900C3F'], \
            debug=True, legend=False)

    if not result:
        sys.exit('error plotting longitudal data.')

    result = plotTimeSeries(fig, (ptime, dtime), (lat, latD), loc, \
            label=['Model Drifter','Observed Drifter'], where=122, \
            title='Latitudal Comparison', \
            axis_label='Latitude $^\circ$', styles=['#FFC300','#581845'], \
            debug=True, legend=False)

    if not result:
        sys.exit('error plotting latitudal data.')

    if args.s:
        print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def spatialProbability(pytkl, drift, ncfile, sim, loc, row, args, saveplot):
    """
    Creates a color-varying spatial plot of the model drifters.
    """
    # check whether or not to overwrite
    savename = row[0][:-4] + '_prob.png'

    if args.s and not args.w:
        if osp.exists(saveplot + savename):
             return

    print 'calculating domain frequencies...'
    elems = np.array(pytkl.variables['indomain'][:])

    # calculate spatial probability density
    # count the number of occurrences of elements
    freq = np.zeros(ncfile.Grid.nele)
    idx, cnt = np.unique(elems, return_counts=True)
    freq[idx] = np.log(cnt)

    # freq = freq.normalize(freq[:,np.newaxis], axis=0).ravel()
    fig = createColorMap(ncfile, freq/np.amax(freq), mesh=False, \
            label = "Probability of Model Drifter Transit",
            title = "Spatially-Varying Probability Plot of Model Drifter\n" + \
            'BFRIC={} | Filename={}'.format(str(bfric), row[0][:-4]))

    plt.scatter(drift.Variables.lon, drift.Variables.lat, c='k', lw=0)

    if args.s:
        print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


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
            outpath = PATH2OUT

            print 'looking for ncfile...'
            if not osp.exists(filename):
                continue

            # define output path
            if args.n:
                outpath += 'n{}'.format(args.n[0])
            if args.r:
                outpath += '_r{}'.format(args.r[0])
            if args.t:
                outpath += '_t'
            if args.g:
                outpath += '_g'

            if outpath != '/':
                outpath += '/'

            outpath += loc + '_' + sim
            if args.n:
                outpath += '_n{}'.format(args.n[0])
            if args.r:
                outpath += '_r{}'.format(args.r[0])
            if args.t:
                outpath += '_t'
            if args.g:
                outpath += '_g'
            outpath += '/'

            if not osp.exists(outpath):
                os.makedirs(outpath)

            # get starting locations and timesteps of drifters
            indata = np.genfromtxt(start_info, dtype=None)
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

                if args.t:
                    print 'randomizing starting times...'
                    intime = row[3] + np.random.choice([-1,1])
                else:
                    intime = row[3]

                # if the run exists skip it
                if not osp.exists(savedir):
                    # set options of drifters
                    # note: interpolation ratio is how many timesteps per
                    # model timestep to linearly interpolate nc data
                    # output ratio is how often to output particle potision
                    options={}
                    options['starttime']=intime
                    options['endtime']=row[4]
                    options['interpolationratio']=60
                    options['outputratio']=2
                    options['ncformat']='NETCDF4_CLASSIC'
                    options['useLL']=True
                    options['layer']=0
                    options['gridDim']='2D'
                    if args.g:
                        options['awgn']=True
                    options['projstr']='lcc +lon_0=-64.55880 +lat_0=41.78504 '+ \
                                   '+lat_1=39.69152 +lat_2=43.87856'

                    # run mitchell's particle tracker
                    print 'tracking pyticles...'
                    start = time.clock()
                    mypy=pyticle(filename, inlocs, savedir, options=options)
                    mypy.run()
                    print('run in: %f' % (time.clock() - start))

                # open pytkl structure
                print 'opening pytkl file...'
                pytkl = nc.Dataset(savedir, 'r', format='NETCDF4_CLASSIC')

                print 'creating time window...'
                # this is in julian time
                tModel = ncfile.Variables.julianTime
                tPytkl = pytkl.variables['time'][:]

                win1 = (np.abs(tModel - tPytkl.min())).argmin()
                win2 = (np.abs(tModel - tPytkl.max())).argmin()

                # calculate time norm
                tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
                print 'opening drifter file...'
                drift = Drifter(path2drift + row[0], debug=False)

                # do things based on command line args
                if 0 not in args.p:
                    save = PATH2PLOT
                    if args.n:
                        save += 'n{}'.format(args.n[0])
                    if args.r:
                        save += '_r{}'.format(args.r[0])
                    if args.t:
                        save += '_t'
                    if args.g:
                        save += '_g'

                    if save[-1] != '/':
                        save += '/'

                    save = save + loc + '_' + sim
                    if args.n:
                        save += '_n{}'.format(args.n[0])
                    if args.r:
                        save += '_r{}'.format(args.r[0])
                    if args.t:
                        save += '_t'
                    if args.g:
                        save += '_g'
                    save += '/'

                    if 1 in args.p:
                        plotTracks(row, pytkl, drift, ncfile, tideNorm, sim, \
                                    loc, args, save)
                    if 2 in args.p:
                        geographicError(pytkl, drift, ncfile, sim, loc, row, \
                                    args, save)
                    if 3 in args.p:
                        spatialProbability(pytkl, drift, ncfile, sim, loc, \
                                    row, args, save)

                if args.a:
                    analyseTracks(pytkl,drift,ncfile,sim,loc,row,args)

            print '...all done!'
