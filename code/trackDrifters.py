#! /usr/env/python2.7

"""
Creates plots and computes stats for pyticle tracker output.

usage: python trackDrifters.py [-p -a -s -w] -r N -n N

cmd line args:
    -p :: plot
    -a :: compute stats
    -s :: save plots
    -v :: add verbosity
    -o :: overwrite current model and plots
    -b :: select a bottom friction
    -l :: select a location
    -d :: select a simulation file date
    -n :: number of particles (int)
    -r :: radius of initialisation (float)
    -w :: add diffusion in calculations from wiener process
    -f :: diffusion fudge factor
    -t :: change the starting time step
"""


# lib imports
from pyticle_tracker import pyticle
from __future__ import division
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
import pandas as pd
import pylab

# local imports
from createColorMap import createColorMap
from drifterUtils import *
from drifterPlotUtils import plotTimeSeries

LOC = ['DG', 'GP', 'PP']
BFRIC = '0.015'
SIM = ['2014_Aug_12_3D', '2013_Aug_08_3D', '2013_Aug_01_3D', '2012_Jun_05_3D', \
       '2012_Jul_05_3D', '2013_Oct_10_3D', '2013_Nov_05_3D', '2013_Nov_06_3D',
       '2013_Aug_02_3D']

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
            + '\n\t2 :: mean lon/lat error' + '\n\t3 :: log-probability density' \
            + '\n\t4 :: pyticle dispersion' + '\n\t5 :: individual error boxplots')
    parser.add_argument('-s', action='store_true', help='save plots.')
    parser.add_argument('-a', action='store_true', help='run analysis.')
    parser.add_argument('-o', action='store_true', help='overwrite data.')
    parser.add_argument('-v', action='store_true', help='verbose output.')
    parser.add_argument('-b', nargs=1, choices=('0.009','0.012','0.015'), \
            help='select a bottom friction.', default='0.015', type=str)
    parser.add_argument('-l', nargs='*', choices=LOC, \
            help='select a location tag.')
    parser.add_argument('-d', nargs='*', choices=SIM,\
            metavar='YYYY_Mmm_DD_3D', help='select a simulation date.')
    multiple = parser.add_argument_group('pyticle options')
    multiple.add_argument('-w', nargs=1, metavar='#', type=float, \
            help='add diffusion from wiener process. requires a seed')
    multiple.add_argument('-f', nargs=1, help='diffusion fudge factor', \
            metavar='ff', type=float, default=1)
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

    if args.s and not args.o:
        if osp.exists(saveplot + savename):
             return

    if args.v:
        print 'preparing to create color map...'
    fig = createColorMap(ncfile, tideNorm[0,:], mesh=False, \
            title = 'Pyticle Track for ' + row[0][:-4],
                        label = 'Mean Velocity Norm (m/s)')

    # add scatter plots of data
    lon = pytkl.variables['lon'][:]
    lat = pytkl.variables['lat'][:]
    if args.v:
        print 'adding scatter plot...'
    pytkl_path = plt.scatter(lon, lat, c='m', lw=0, alpha=0.6, marker="^", s=10)

    lonD = drift.Variables.lon
    latD = drift.Variables.lat
    drift_path = plt.scatter(lonD, latD, c='k', lw=0, s=25)

    plt.legend([drift_path, pytkl_path],['Drifter', 'Model'])

    if args.s:
        if args.v:
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


def geographicError(pytkl, drift, ncfile, sim, loc, row, args, \
                    saveplot, tn=None):
    """
    Calculates errors in lon/lat coordinates and velocities graphically.

    pytkl :: pyticle object
    ncfile :: FVCOM mode
    drift :: drifter data
    sim :: sim date
    loc :: sim location
    row :: drifter file name
    args :: command line args
    saveplot :: name of save plot dir
    tn :: optional tideNorm, prints side by side
    """
    # check whether or not to overwrite
    savename = row[0][:-4] + '_lonlatErr.png'

    if args.s and not args.o:
        if osp.exists(saveplot + savename):
             return

    if args.v:
        print 'calculating lon/lat error...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat
    uD = drift.Variables.u
    vD = drift.Variables.v
    u = np.array(pytkl.variables['u'][:])
    v = np.array(pytkl.variables['v'][:])

    dtime = drift.Variables.matlabTime
    ptime = mjd2num(pytkl.variables['time'][:])

    if args.v:
        print str(lon.shape[1]) + ' pyticles...'
        print 'no. pts pytkl: ' + str(len(lon))
        print 'no. pts drift: ' + str(len(lonD))
        print 'broadcasting...'

    # truncate longer list to broadcast
    if len(lon) > len(lonD):
        idx = broadcast(ptime, dtime)
        lon = lon[idx]
        lat = lat[idx]
        u = u[idx]
        v = v[idx]
        ptime = ptime[idx]
    elif len(lonD) > len(lon):
        idx = broadcast(dtime, ptime)
        dtime = dtime[idx]
        lonD = lonD[idx]
        latD = latD[idx]
        uD = uD[idx]
        vD = vD[idx]

    dtime = [dn2dt(x) for x in dtime]
    ptime = [dn2dt(x) for x in ptime]

    # if multiple particles, take the mean
    if len(lon.shape) > 1:
        lonErr = np.squeeze(np.array([np.abs(row-lonD) for row in lon.T])).T
        latErr = np.squeeze(np.array([np.abs(row-latD) for row in lat.T])).T
        uErr = np.squeeze(np.array([np.abs(row-uD) for row in u.T])).T
        vErr = np.squeeze(np.array([np.abs(row-vD) for row in v.T])).T

        lonDiff = np.mean(lonErr, axis=1)
        latDiff = np.mean(latErr, axis=1)
        uDiff = np.mean(uErr, axis=1)
        vDiff = np.mean(vErr, axis=1)
        lab1 = ['Mean Lon. Diff.', 'Mean Lat. Diff.']
        lab2 = ['Mean u Diff.', 'Mean v Diff.']
    else:
        lonDiff = np.abs(lonD - lon)
        latDiff = np.abs(latD - lat)
        uDiff = np.abs(uD - u)
        vDiff = np.abs(vD - v)
        lab1 = ['Lon. Diff.', 'Lat. Diff.']
        lab2 = ['u Diff.', 'v Diff.']

    if tn != None:
        ## **Needs work done to ensure proper scaling of subplots on figure!
        where1 = 221
        where2 = 223
        fig = createColorMap(ncfile, tn[0,:], mesh=False, \
                title='Trajectory for {}'.format(row[0][:-4]), \
                label='Mean Velocity Norm (m/s)')
        fig.suptitle('Positional Timeseries :: ' + \
                 'BFRIC={} | Filename={}'.format(str(bfric), row[0][:-4]),
                 fontsize=10)
        plt.scatter(lon, lat)
        plt.scatter(lonD, latD)
    else:
        fig = plt.figure()
        where1 = 211
        where2 = 212

    if args.v:
        print 'creating figures...'
    result, axis = plotTimeSeries(fig, (ptime, dtime), (lonDiff, latDiff), loc, \
            label=lab1, where=where1, \
            title='Longitudal/Latitudal Timeseries Error', \
            axis_label='Lon/Lat Error ($^\circ$)', styles=['#1DA742','#900C3F'], \
            debug=True, legend=True)

    if not result:
        sys.exit('error plotting longitudal/latitudal data.')

    result, axis = plotTimeSeries(fig, (ptime, dtime), (uDiff, vDiff), loc, \
            label=lab2, where=where2, \
            title='u/v-Velocity Timeseries Error', \
            axis_label='u/v-Velocity Error (m/s)', styles=['#FFC300','#581845'], \
            axx=axis, debug=True, legend=True)

    if not result:
        sys.exit('error plotting velocity data.')

##    result, axis = plotTimeSeries(fig, (ptime, dtime), (lon, lonD), loc, \
##            label=lab, where=where1, \
##            title='Longitudal Timeseries Comparison', \
##            axis_label='Longitude ($^\circ$)', styles=['#1DA742','#900C3F'], \
##            debug=True, legend=leg)
##
##    if not result:
##        sys.exit('error plotting longitudal data.')
##
##    result, axis = plotTimeSeries(fig, (ptime, dtime), (lat, latD), loc, \
##            label=['Model Drifter','Observed Drifter'], where=where2, \
##            title='Latitudal Timeseries Comparison', \
##            axis_label='Latitude ($^\circ$)', styles=['#FFC300','#581845'], \
##            axx=axis, debug=True, legend=leg)
##
##    if not result:
##        sys.exit('error plotting latitudal data.')

    if args.s:
        if args.v:
            print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def boxplotError(pytkl, drift, ncfile, sim, loc, row, args, saveplot):
    """
    Plots a box plot for each model-drifter time series error.
    """
    # check whether or not to overwrite
    savename = row[0][:-4] + '_lonlatErr.png'

    if args.s and not args.o:
        if osp.exists(saveplot + savename):
             return

    if args.v:
        print 'calculating lon/lat error...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat
    uD = drift.Variables.u
    vD = drift.Variables.v
    u = np.array(pytkl.variables['u'][:])
    v = np.array(pytkl.variables['v'][:])

    dtime = drift.Variables.matlabTime
    ptime = mjd2num(pytkl.variables['time'][:])

    if args.v:
        print str(lon.shape[1]) + ' pyticles...'
        print 'no. pts pytkl: ' + str(len(lon))
        print 'no. pts drift: ' + str(len(lonD))
        print 'broadcasting...'

    # truncate longer list to broadcast
    if len(lon) > len(lonD):
        idx = broadcast(ptime, dtime)
        lon = lon[idx]
        lat = lat[idx]
        u = u[idx]
        v = v[idx]
        ptime = ptime[idx]
    elif len(lonD) > len(lon):
        idx = broadcast(dtime, ptime)
        dtime = dtime[idx]
        lonD = lonD[idx]
        latD = latD[idx]
        uD = uD[idx]
        vD = vD[idx]

    dtime = [dn2dt(x) for x in dtime]
    ptime = [dn2dt(x) for x in ptime]

    if len(lon.shape) > 1:
        lonErr = np.squeeze(np.array([np.abs(row-lonD) for row in lon.T])).T
        latErr = np.squeeze(np.array([np.abs(row-latD) for row in lat.T])).T
        uErr = np.squeeze(np.array([np.abs(row-uD) for row in u.T])).T
        vErr = np.squeeze(np.array([np.abs(row-vD) for row in v.T])).T

    if args.v:
        print 'creating figure...'
    # PANDAS? PYLAB?
    # Is the error normal??
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    bp = ax1.boxplot(lonErr, showmeans=True, patch_artist=True)
    for box in bp['boxes']:
        box.set(color='#1C2833',linewidth=2)
        box.set(facecolor='lightblue', alpha=1)
    for whisker in bp['whiskers']:
        whisker.set(color='#1C2833', linewidth=2)
    for cap in bp['caps']:
        cap.set(color='#1C2833', linewidth=2)
    for median in bp['medians']:
        median.set(color='#1C2833')
    for flier in bp['fliers']:
        flier.set(marker='+')
        flier.set(color='k', alpha=0.5)
    for mean in bp['means']:
        mean.set(marker='o')
        mean.set(color='k')
    ax1.set_ylabel('Longitudal Error')
    ax1.set_title('Particle-Drifter Error')

    ax2 = fig.add_subplot(212)
    bp = ax2.boxplot(lonErr, showmeans=True, patch_artist=True)
    for box in bp['boxes']:
        box.set(color='#1C2833',linewidth=2)
        box.set(facecolor='lightblue', alpha=1)
    for whisker in bp['whiskers']:
        whisker.set(color='#1C2833', linewidth=2)
    for cap in bp['caps']:
        cap.set(color='#1C2833', linewidth=2)
    for median in bp['medians']:
        median.set(color='#1C2833')
    for flier in bp['fliers']:
        flier.set(marker='+')
        flier.set(color='k', alpha=0.5)
    for mean in bp['means']:
        mean.set(marker='o')
        mean.set(color='k')
    ax2.set_xlabel('Particle Label')
    ax2.set_ylabel('Latitudal Error')

    if args.s:
        if args.v:
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

    if args.s and not args.o:
        if osp.exists(saveplot + savename):
             return

    if args.v:
        print 'calculating domain frequencies...'
    elems = np.array(pytkl.variables['indomain'][:])

    # calculate spatial probability density
    # count the number of occurrences of elements
    freq = np.zeros(ncfile.Grid.nele)
    idx, cnt = np.unique(elems, return_counts=True)
    freq[idx] = np.log(cnt)

    # freq = freq.normalize(freq[:,np.newaxis], axis=0).ravel()
    fig = createColorMap(ncfile, freq/np.amax(freq), mesh=False, \
            label = "Log Probability of Model Drifter Transit",
            title = "Spatially-Varying Probability Plot of Model Drifter\n" + \
            'BFRIC={} | Filename={}'.format(str(bfric), row[0][:-4]))

    plt.scatter(drift.Variables.lon, drift.Variables.lat, c='k', lw=0)

    if args.s:
        if args.v:
            print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def dispersionPlots(pytkl, drift, ncfile, sim, loc, row, args, saveplot):
    """
    Plots mean long./lat. and relative displacement of the pyticles
    against time, with a 95% confidence interval.
    """
    # check whether or not to overwrite
    savename = row[0][:-4] + '_disp.png'

    if args.s and not args.o:
        if osp.exists(saveplot + savename):
             return

    if args.v:
        print 'collecting lon/lats...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat

    dtime = drift.Variables.matlabTime
    ptime = mjd2num(pytkl.variables['time'][:])
    dtime = [dn2dt(x) for x in dtime]
    ptime = [dn2dt(x) for x in ptime]

    if args.v:
        print str(lon.shape[1]) + ' pyticles...'
        print 'no. pts pytkl: ' + str(len(lon))
        print 'no. pts drift: ' + str(len(lonD))

    if args.v:
        print 'calculating statistics...'
    lonM = np.mean(lon, axis=1)
    latM = np.mean(lat, axis=1)
    lonSD = np.std(lon, axis=1)
    latSD = np.std(lat, axis=1)

    dlonM = lonM - lonM[0]
    dlatM = latM - latM[0]

    disp = np.sqrt(np.square(dlonM*79868) + np.square(dlatM*111117))
    dispSD = np.sqrt(np.square(lonSD*79868) + np.square(latSD*111117))

    t = sp.stats.t.ppf(0.95, len(lonM)-1)
    CI = t*dispSD

    # plotting
    if args.v:
        print 'creating plot...'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ptime, disp,'o', color='#34495E', markersize=5, markeredgewidth=1,\
            markeredgecolor='#2E4053', markerfacecolor='#34495E')
    ax.fill_between(ptime, disp+CI, disp-CI, color='#F39C12', edgecolor='')
    # minor hack for labelling CI
    ax.plot(ptime, disp+CI, color='#F39C12', label='95% Confidence Interval')
    ax.set_ylabel('Displacement (m)')
    ax.set_xlabel('Time (HH:MM:SS)')
    ax.set_title('Dispersion of Pyticles')
    plt.gcf().autofmt_xdate()
    plt.grid(True)
    plt.tight_layout()
    plt.ylim(0, np.max(disp)+np.max(dispSD))

    if args.s:
        if args.v:
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

            if args.v:
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
            if args.w:
                outpath += '_w{}'.format(args.w[0])
                if args.f[0] != 1:
                    outpath += '_ff{}'.format(args.f[0])

            if outpath != '/':
                outpath += '/'

            outpath += loc + '_' + sim
            if args.n:
                outpath += '_n{}'.format(args.n[0])
            if args.r:
                outpath += '_r{}'.format(args.r[0])
            if args.t:
                outpath += '_t'
            if args.w:
                outpath += '_w{}'.format(args.w[0])
                if args.f[0] != 1:
                    outpath += '_f{}'.format(args.f[0])

           outpath += '/'

            if not osp.exists(outpath):
                os.makedirs(outpath)

            # get starting locations and timesteps of drifters
            indata = np.genfromtxt(start_info, dtype=None)
            if args.v:
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
            if args.v:
                print 'opening fvcom file...'
            ncfile = FVCOM(filename, debug=False)
            if args.v:
                print 'calculating ebb/flood split at centre of location...'
            fI, eI, _, _ = ncfile.Util2D.ebb_flood_split_at_point(centre[0], \
                            centre[1])
            if args.v:
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
                    if args.v:
                        print 'starting with {} particles...'.format(num)
                else:
                    num = 1

                if args.r:
                    # 1 deg lat = 111117 m
                    # 1 deg lon =  79868 m
                    inlocs = np.vstack((np.random.uniform( \
                            inloc[0]-args.r[0]/79868.0, \
                            inloc[0]+args.r[0]/79868.0, num), \
                            np.random.uniform(inloc[1]-args.r[0]/111117.0, \
                                              inloc[1]+args.r[0]/111117.0, \
                                              num))).T

                    if args.v:
                        print 'randomizing starting locations...'

                if args.t:
                    if args.v:
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
                    if args.w:
                        options['diffusion']=True
                        options['seed']=args.w[0]
                        options['ff']=args.f[0]

                    options['projstr']='lcc +lon_0=-64.55880 +lat_0=41.78504 '+ \
                                   '+lat_1=39.69152 +lat_2=43.87856'

                    # run mitchell's particle tracker
                    if args.v:
                        print 'tracking pyticles...'
                    start = time.clock()
                    mypy=pyticle(filename, inlocs, savedir, options=options)
                    mypy.run()
                    print('run in: %f' % (time.clock() - start))

                # open pytkl structure
                if args.v:
                    print 'opening pytkl file...'
                pytkl = nc.Dataset(savedir, 'r', format='NETCDF4_CLASSIC')

                if args.v:
                    print 'creating time window...'
                # this is in julian time
                tModel = ncfile.Variables.julianTime
                tPytkl = pytkl.variables['time'][:]

                win1 = (np.abs(tModel - tPytkl.min())).argmin()
                win2 = (np.abs(tModel - tPytkl.max())).argmin()

                # calculate time norm
                tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
                if args.v:
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
                    if args.w:
                        save += '_w{}'.format(args.w[0])
                         if args.f[0] != 1:
                            save += '_ff{}'.format(args.f[0])

                    if save[-1] != '/':
                        save += '/'

                    save = save + loc + '_' + sim
                    if args.n:
                        save += '_n{}'.format(args.n[0])
                    if args.r:
                        save += '_r{}'.format(args.r[0])
                    if args.t:
                        save += '_t'
                    if args.w:
                        save += '_w{}'.format(args.w[0])
                        if args.f[0] != 1:
                            save += '_ff{}'.format(args.f[0])

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
                    if 4 in args.p:
                        dispersionPlots(pytkl, drift, ncfile, sim, loc, \
                                    row, args, save)
                    if 5 in args.p:
                        boxplotError(pytkl, drift, ncfile, sim, loc, row, args, \
                                    save)

                if args.a:
                    analyseTracks(pytkl,drift,ncfile,sim,loc,row,args)

            print '...all done!'
