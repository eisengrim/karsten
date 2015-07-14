#!/usr/env/python2.7
#! encoding: utf-8

"""
Filename = drifterBias.py
Author = Kody Crowell
Version = 1.0

Given a directory containing Drifter MATLAB files, and a couple FVCOM ncfiles,
this program calculates the mean and standard deviation of the speed biases.
Plots are created and shown of the bias vs each drifter and model speed.

The drifter directory is defined using the region tag 'GP', 'DG' or 'PP'.
"""

# Library Imports
import sys, os
import argparse as arp
import numpy as np
import scipy as sp
import scipy.special as sps
import scipy.io as sio
import os.path as osp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as tic
from mpl_toolkits.basemap import Basemap
# import seaborn as sns
from pyseidon import *

PATH_TO_SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
PATH_TO_OBS="/EcoII/acadia_uni/workspace/observed/"
GRID='dngridCSR'

def dn2dt(datenum, debug=False):
    """
    Convert matlab dateum to python datetime.

    input:
        - matlab datenum
    returns:
        - python datetime
    """
    if debug:
        print 'converting matlab datenum to python datetime...'

    return datetime.fromordinal(int(datenum)) + \
           timedelta(days=datenum % 1) - \
           timedelta(days=366)

def driftTimes(name, debug=False):
    """
    Identify the timespans for each drifter file. Adapted from Jon Smith's
    timespan.py.

    input:
        - drifter file name
    returns:
        - starting and ending time
    """
    # iterate through the drifter files
    start_time = []
    end_time = []

    try:
        if debug:
            print 'loading file {}...'.format(name)
        drft = sio.loadmat(name)
        times = drft['gps_observation'][0][0][0][0]
    except KeyError:
        times = drft['time'][0][0][0][0]

    # grab the times and convert them to strings
    start, end = times[0], times[-1]
    return start, end


def calculateBias(ncfile, files, loc, debug=False):
    """
    Compiles necessary data, calculates individual biases and
    takes the mean and standard deviation. Also calculates the
    bias for ALL drifter files. Returns the mean biases, the drifter
    numbers, the model and observed speeds the the std deviations.

    input:
        - ncfile : FVCOM object
        - files : list of matlab filenames in directory
        - loc : location tag
        - sim_name : name of simulation directory used
    """

    # sns.set(font="serif")
    if debug:
        print '{} drifters will be analysed...'.format(len(files))

    drifters = {}
    all_bias = []
    all_mean = []
    all_sdev = []
    obs_speed = []
    mod_speed = []
    obs_uspeed = []
    mod_uspeed = []
    o_lon = []
    o_lat = []

    for i, fname in enumerate(files, start=1):
        drifters[i] = fname
        if debug:
            print 'creating drifter object...'
        drift = Drifter(fname, debug=False)

        # create validation structure
        if debug:
            print 'creating validation object...'
        valid = Validation(drift, ncfile, flow='sf', debug=False)

        if debug:
            print 'extracting information...'

        # extract information
        mTimes = valid.Variables.struct['mod_time']
        oU = valid.Variables.struct['obs_timeseries']['u']
        oV = valid.Variables.struct['obs_timeseries']['v']
        mU = valid.Variables.struct['mod_timeseries']['u']
        mV = valid.Variables.struct['mod_timeseries']['v']
        olon = valid.Variables.struct['lon']
        olat = valid.Variables.struct['lat']
        uspeedS = np.asarray(np.sqrt(mU**2 + mV**2))
        uspeedO = np.asarray(np.sqrt(oU**2 + oV**2))

        if debug:
            print 'calculating signed speeds...'

        # datetimes = np.asarray([dn2dt(time) for time in mTimes])

        if not uspeedS.shape == uspeedO.shape:
            if debug:
                print 'drifter {} does not have similar-shaped speeds...'
            continue

        speedO = uspeedO * np.sign(oV)
        speedS = uspeedS * np.sign(mV)
        diffs = np.subtract(uspeedS, uspeedO)

        mean_bias = np.mean(diffs)
        sdev_bias = np.std(diffs)
        all_sdev.append(sdev_bias)
        all_mean.append(mean_bias)
        all_bias.extend(diffs)

        obs_speed.extend(speedO)
        mod_speed.extend(speedS)
        obs_uspeed.extend(uspeedO)
        mod_uspeed.extend(uspeedS)
        o_lon.extend(olon)
        o_lat.extend(olat)

    return drifters, all_mean, all_sdev, obs_speed, mod_speed, obs_uspeed, \
         mod_uspeed, all_bias, o_lat, o_lon


def parseArgs():
    """
    Parses, identifies and resolves command line arguments.
    """

    parser = arp.ArgumentParser(prog='drifterBias.py', description="Opens " \
            + "an FVCOM file and multiple drifter files and calcualtes the "\
            + "mean and standard deviation of the biases.")
    # creates object to be parsed. adds command line optional arguments
    parser.add_argument("--debug", "-v", "--verbose", action="store_true", \
            help="increases output verbosity.")
    parser.add_argument("-V", '--version', action='version', version='v1.0')
    # required options bad form, users expect options to be optional
    require = parser.add_argument_group('required flag arguments')
    require.add_argument("--loc", '-l', help="define the region to operate "\
            + "within.", nargs=1, choices=('GP', 'DG', 'PP'), required=True)
    parser.add_argument("--dir", '-D', nargs=1, help="defines an FVCOM " \
            + "directory. Name is in the form YYYY_Mmm_DD_3D", \
            metavar='dirname', type=str)
    parser._optionals.title = 'optional flag arguments'

    args = parser.parse_args()

    if args.debug:
        # identifies program options
        print '-verbosity turned on.-'

        if args.loc:
            print '\tlocation tag set to {}...'.format(args.loc)
        if args.dir:
            print '\tfvcom directory set to {}...'.format(args.dir)

    if not args.loc:
        sys.exit('a location tag is needed. type --help for more info.')

    return args


def setOptions(args):
    """
    Interprets the options returned by argument parser.
    """
    obs_dir = PATH_TO_OBS + args.loc[0] + '/Drifter/'

    if args.debug:
        print 'looking for drifter directory...'

    if not osp.exists(obs_dir) or not osp.isdir(obs_dir):
        sys.exit('drifter directory not found.')
    elif args.debug:
        print 'drifter directory successfully found.'
        print '\tgathering all files...'

    matfiles = [obs_dir + file for file in os.listdir(obs_dir)]
    if len(matfiles) == 0:
        sys.exit('no files found in drifter directory.')
    elif args.debug:
        print '\tall drifter files found.'

    if args.debug:
        print 'looking for fvcom directory(s)...'

    # locate given fvcom file
    if args.dir:
        sim_path = [PATH_TO_SIM + args.loc[0] + '/' + args.dir[0]]
        if not osp.exists(sim_path) or not osp.isdir(sim_path):
            sys.exit('the directory {} could not be located.'.format(sim_path))
        elif args.debug:
            print '\tfvcom directory found. \n\tloading nc file...'

        sim_path += '/output/subdomain_' + args.loc[0] + '1_0001.nc'
        if not osp.exists(sim_path) or not osp.isfile(sim_path):
            sys.exit('fvcom file not in directory.')
        elif args.debug:
            print '\tfvcom file successfully located.'
    else:
        dirs = os.listdir(PATH_TO_SIM + args.loc[0] + '/')
        sim_path = [PATH_TO_SIM + args.loc[0] + '/' + file + \
            '/output/subdomain_' + args.loc[0] + '1_0001.nc' for file in dirs]

        for path in sim_path:
            if not osp.exists(path) or not osp.isfile(path):
                sys.exit('fvcom file {} is not found.'.format(path))
                sim_path.remove(path)

        if len(sim_path) == 0:
            sys.exit('no ncfiles found in directory.')
        elif args.debug:
            print '\tall ncdirectories found.'

    return args.loc[0], sim_path, obs_dir, matfiles


def createPlots(speedS, speedO, mean, stdev, bias, uspdS, uspdO, lon, lat, \
        model, debug=debug):
    """
    Creates a bunch of plots.
    """
    # estimate cube ratio
    ratio = sps.cbrt(np.mean(np.power(uspdO,3))/np.mean(np.power(uspdS,3)))
    if debug:
        print 'speed ratio is: {}'.format(ratio)

    # create plots
    if debug:
        print 'creating plots...'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(speedS, bias, alpha=0.25)
    ax.set_xlabel('Model Speed (m/s)')
    ax.set_ylabel('Bias')

    # determine line of best fit
    if debug:
        print '\tdetermining line of best fit...'

    par = np.polyfit(speedS, bias, 1)
    m = par[-2]
    b = par [-1]

    variance = np.var(bias)
    residuals = np.var([(m*xx + b - yy)  for xx,yy \
            in zip(speedS, bias)])
    Rsqr = np.round(1-residuals/variance, decimals=5)
    if debug:
        print '\tR^2 value for bias plot is {}...'.format(Rsqr)
    plt.hold('on')
    ax.plot(speedS, m*speedS+b, 'r-')
    plt.grid(True)

    # plot cube speeds
    fi = plt.figure()
    ax1 = fi.add_subplot(111)
    speedS3 = np.power(speedS, 3)
    speedO3 = np.power(speedO, 3)
    ax1.scatter(speedS3, speedO3, alpha=0.25)
    ax1.set_xlabel('Model Speed (m/s)')
    ax1.set_ylabel('Drifter Speed (m/s)')

    coeff = np.polyfit(speedS3, speedO3, 1)
    m = coeff[-2]
    b = coeff[-1]
    if debug:
        print 'coeffs for cube plot are: \n\tm={}\n\tb={}'.format(m,b)
    plt.hold('on')
    ax1.plot(speedS3, m*speedS3+b)
    variance = np.var(speedO3)
    residuals = np.var([(m*xx + b - yy)  for xx,yy \
            in zip(speedS3, speedO3)])
    Rsqr = np.round(1-residuals/variance, decimals=5)
    if debug:
        print 'R^2 for cube plot is {}'.format(Rsqr)

    # plot bias v drifter
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(np.arange(1,num_drift+1), mean, 'go')
    ax2.set_ylabel('Bias')
    ax2.set_xlabel('Drifter')
    plt.hold('on')
    ax2.axhline(y=np.mean(bias), linewidth=2)
    plt.grid(True)

    # spatial plot of ratios...
    if debug:
        print 'starting spatial plot...'
    glon = model.Grid.lon
    glat = model.Grid.lat
    if debug:
        print 'computing bounding boxes...'
    bounds = [np.min(glon), np.max(glon), np.min(glat), np.max(glat)]

    if not hasatr(model.Grid, 'triangleLL'):
        tri = Tri.Triangulation(glon, glat, triangles=model.Grid.trinodes)
    else:
        tri = model.Grid.triangleLL

    if debug:
        print 'creating subplot...'
    # gif = plt.figure()
    # xa = gif.add_subplot(111, aspect=(1.0/np.cos(np.mean(glat) * np.pi/180.0)))

    if debug:
        print 'computing colors...'
        var3 = np.divide(speedO3, speedS3)
        color = [str(ratio/255.0) for ratio in var3]
    if debug:
        print 'creating map...'

    map = Basemap(projection='merc', resolution = 'h', area_thresh = 0.1, \
            llcrnrlon=bounds[0], llcrnrlat=bounds[2], \
            urcrnrlon=bounds[1], urcrnrlat=bounds[3])

    map.drawcoastlines()
    map.drawcountries()
    map.bluemarble()

    f = map.scatter(lon, lat, latlon=True, s=100, c=color)
    cbar = map.colorbar(f, ax=map)
    cbar.set_label('Speed^3 Ratio', rotation=-90, labelpad=30)
    scale = 1

    ticks = tic.FuncFormatter(lambda glon, pos: '{0:g}'.format(lon/scale))
    map.xaxis.set_major_formatter(ticks)
    map.yaxis.set_major_formatter(ticks)
    map.set_xlabel('Longitude')
    map.set_ylabel('Latitiude')


if __name__ == '__main__':

    # parse the command line args and identify parameters
    print '\nparsing command line options...'
    args = parseArgs()
    debug = args.debug

    loc, sim_path, obs_dir, obs_files = setOptions(args)

    if debug:
        print '\n--parameters selected--'
        print 'location: ', loc, '\nsim_path: ', sim_path, \
                '\nobs_path: ', obs_dir

        print '\nloading fvcom object...'

    drifters = {}
    all_mean = []
    all_std = []
    all_speedO = []
    all_speedS = []
    all_uspeedO = []
    all_uspeedS = []
    all_bias = []
    all_lon = []
    all_lat = []
    num_drift = 0

    for dir_name in sim_path:
        ncfile = FVCOM(dir_name, debug=False)
        if debug:
            print 'ncfile for {} loaded.'.format(dir_name)

        # find the relevant time window to work in
        mTimes = ncfile.Variables.matlabTime[:]
        mStart, mEnd= float(mTimes[0]), float(mTimes[-1])

        if debug:
            print 'model time is from {} to {}.'.format(mStart, mEnd)

        # from given drifter files, find files in fvcom runtime window
        if debug:
            print 'gathering all drifter files in model runtime window...'

        files = []
        for matfile in obs_files:
            dStart, dEnd = driftTimes(matfile, debug=debug)
            dStart, dEnd = float(dStart), float(dEnd)
            if dStart > mStart and mEnd > dEnd:
                files.append(matfile)
                if debug:
                    print 'file {} is within runtime window.'.format(matfile)
        if not files:
            sys.exit('drifters given are not within model runtime window.')

        drift, mean, std, speedO, speedS, uspdO, uspdS, bias, lat, lon \
                = calculateBias(ncfile, files, loc, debug=debug)

        if debug:
            print 'adding to cumulative data'

        drifters['dir_name'] = drift
        all_mean.extend(mean)
        all_std.extend(std)
        all_speedS.extend(speedS)
        all_speedO.extend(speedO)
        all_uspeedS.extend(uspdS)
        all_uspeedO.extend(uspdO)
        all_bias.extend(bias)
        all_lat.extend(lat)
        all_lon.extend(lon)
        num_drift = num_drift + len(drift)

    if debug:
        print 'creating as numpy arrays...'

    speedS = np.asarray(all_speedS)
    speedO = np.asarray(all_speedO)
    mean = np.asarray(all_mean)
    stdev = np.asarray(all_std)
    bias = np.asarray(all_bias)
    uspdS = np.asarray(all_uspeedS)
    uspdO = np.asarray(all_uspeedO)
    lat = np.asarray(all_lat)
    lon = np.asarray(all_lon)

    createPlots(speedS, speedO, mean, stdev, bias, uspdS, uspdO, lat, lon, \
            ncfile, debug=debug)

    plt.show()
    if debug:
        print '...all done!'
