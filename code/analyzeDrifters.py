#!/usr/env/python2.7
#! encoding: utf-8

"""
Filename = analyzeDrifters.py
Author = Kody Crowell
Version = 1.4

This code is an updated version of drifterBias.py, tasked with running various
stats that are returned to be used in the command line, and then plotted.
"""

# library imports
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
import seaborn as sns
from pyseidon import *
from interpolation_utils import *
import sklearn.preprocessing as skp

PATH2SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
PATH2OBS="/EcoII/acadia_uni/workspace/observed/"

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
            print 'examining drifter file...'
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

    if debug:
        print '{} drifters will be analysed...'.format(len(files))

    drifters = {}
    all_bias = []
    all_ubias = []
    all_mean = []
    all_sdev = []
    obs_speed = []
    mod_speed = []
    obs_uspeed = []
    mod_uspeed = []
    o_lon = []
    o_lat = []
    lon0 = []
    lat0 = []
    depths = []
    all_erru = []
    all_errv = []
    all_err_mag = []
    all_err_dir = []

    for i, fname in enumerate(files, start=1):
        drifters[i] = fname
        if debug:
            print 'creating drifter object {}...'.format(i)
            print fname
        drift = Drifter(fname, debug=False)
        # create validation structure
        if debug:
            print '\tcreating validation object...'
        try:
            valid = Validation(drift, ncfile, flow='sf', debug=False)
        except IndexError:
            print 'cannot create validation object for drifter %i.' % i
            continue

        if debug:
            print '\textracting information...'

        # extract information
        mTimes = valid.Variables.struct['mod_time']
        oU = valid.Variables.struct['obs_timeseries']['u']
        oV = valid.Variables.struct['obs_timeseries']['v']
        mU = valid.Variables.struct['mod_timeseries']['u']
        mV = valid.Variables.struct['mod_timeseries']['v']
        olon = valid.Variables.struct['lon']
        olat = valid.Variables.struct['lat']
        mlon = ncfile.Grid.lon
        mlat = ncfile.Grid.lat
        mlonc = ncfile.Grid.lonc
        mlatc = ncfile.Grid.latc
        # how to get mlon, mlat?
        uspeedS = np.sqrt(mU**2 + mV**2)
        uspeedO = np.sqrt(oU**2 + oV**2)

        # index finder is only in dvt branch
        # indices = [ncfile.Util2D.index_finder(j,k) for j, k in zip(olon, olat)]
        # c_pts = closest_point(olon, olat, mlon, mlat, mlonc, mlatc, tri??)
        c_pts = closest_points(olon, olat, mlon, mlat)
        #indices = np.where(
        depth = ncfile.Grid.h[indices]

        if debug:
            print '\tcalculating signed speeds...'

        # datetimes = np.asarray([dn2dt(time) for time in mTimes])

        if not uspeedS.shape == uspeedO.shape:
            if debug:
                print 'drifter {} does not have similar-shaped speeds...'
            continue

        if debug:
            print '\tcalculating statistics...'

        speedO = uspeedO * np.sign(oV)
        speedS = uspeedS * np.sign(mV)
        diffs = np.subtract(uspeedS, uspeedO)
        udiffs = np.subtract(speedS, speedO)
        mean_bias = np.mean(diffs)
        sdev_bias = np.std(diffs)
        erru = np.subtract(mU, oU)
        errv = np.subtract(mV, oV)
        err_mag = np.sqrt(erru**2 + errv**2)
        err_dir = np.arctan(np.divide(errv, erru))/np.pi * 180

        all_sdev.append(sdev_bias)
        all_mean.append(mean_bias)
        all_bias.extend(diffs)
        all_ubias.extend(udiffs)
        all_erru.extend(erru)
        all_errv.extend(errv)
        all_err_mag.extend(err_mag)
        all_err_dir.extend(err_dir)
        depths.extend(depth)

        if debug:
            print '\tcompiling data...'
        lon0.append(olon[0])
        lat0.append(olat[0])
        obs_speed.extend(speedO)
        mod_speed.extend(speedS)
        obs_uspeed.extend(uspeedO)
        mod_uspeed.extend(uspeedS)
        o_lon.extend(olon)
        o_lat.extend(olat)

    return drifters, all_mean, all_sdev, obs_speed, mod_speed, obs_uspeed, \
         mod_uspeed, all_bias, o_lat, o_lon, lon0, lat0, all_ubias, depths, \
         erru, errv, all_err_mag, all_err_dir


def parseArgs():
    """
    Parses, identifies and resolves command line arguments. Only required args are
    a region tag and a save path.
    """

    parser = arp.ArgumentParser(prog='drifterBias.py', description="Opens " \
            + "any FVCOM files and drifter files, calculates statistics "\
            + "and generates plots.")
    # adds command line optional arguments
    parser.add_argument("--debug", "-v", "--verbose", action="store_true", \
            help="increases output verbosity.")
    parser.add_argument("-V", '--version', action='version', version='v1.4')
    require = parser.add_argument_group('required')
    require.add_argument("--loc", '-l', help="define the region to operate "\
            + "within.", nargs=1, choices=('GP', 'DG', 'PP'), required=True)
    parser.add_argument("--dir", '-D', nargs='*', help="defines an FVCOM " \
            + "directory. Name is in the form YYYY_Mmm_DD_3D", \
            metavar='dirname', type=str)
    parser.add_argument("--drifter",'-d', nargs='*', help="defines one or " \
            + "more matlab drifter files to be read.", metavar='matfile', \
            type=str, default=False)
    parser.add_argument("--bfric", '-B', help="select a bottom " \
            + "roughness.", nargs=1, choices=('0.009', '0.012', '0.015'), \
            default='0.015', required=True, type=str)
    # save or show data
    parser.add_argument("--noplot", "-x", action="store_true", help='generates ' \
            + 'no plots, just calculates statistics.')
    parser.add_argument("--savepath", '-p', nargs='?', help="defines an  " \
            + "alternate savepath.", metavar='savedir')
    parser._optionals.title = 'optional flag arguments'
    # option to write initial positions
    parser.add_argument("--write", '-w', help='records initial positions of ' \
            + 'drifters to .dat file.', action="store_true")

    args = parser.parse_args()

    if args.debug:
        # identifies program options
        print '-verbosity turned on.-'

        if args.loc:
            print '\tlocation tag set to {}...'.format(args.loc)

    if not args.loc:
        sys.exit('a location tag is needed. type --help for more info.')

    if args.noplot:
        print '\tplots suppressed...'
    if args.write:
        # add write for plots...
        print '\tinitial locations to be recorded...'

    return args


def setOptions(args):
    """
    Interprets the options returned by argument parser.
    """
    obs_dir = PATH2OBS + args.loc[0] + '/Drifter/'

    # look for observed data
    if args.debug:
        print 'looking for drifter directories...'

    if not osp.exists(obs_dir) or not osp.isdir(obs_dir):
        sys.exit('drifter directory not found.')
    elif args.debug:
        print 'drifter directory successfully found.'
        print '\tgathering all files...'

    if args.drifter:
        if args.drifter[0].endswith('.dat') or args.drifter[0].endswith('.txt'):
            matfiles = np.loadtxt(args.drifter[0])
        else:
            matfiles = [obs_dir + file for file in args.drifter]
    else:
        matfiles = [obs_dir + file for file in os.listdir(obs_dir)]

    if len(matfiles) == 0:
        sys.exit('no drifter files found.')
    elif args.debug:
        print '\tall drifter files found.'

    # look for simulated data
    if args.debug:
        print 'looking for fvcom directory(s)...'

    path2sim = PATH2SIM + 'BFRIC_' + args.bfric[0] + '/'

    # locate given fvcom file
    if args.dir:
        dirs = args.dir
    else:
        dirs = os.listdir(path2sim + args.loc[0] + '/')

    sim_path = [path2sim + args.loc[0] + '/' + file + '/output/subdomain_' \
                + args.loc[0] + '1_0001.nc' for file in dirs]

    for path in sim_path:
        if not osp.exists(path) or not osp.isfile(path):
            print '\tfvcom file {} is not found. removing...'.format(path)
            sim_path.remove(path)

    if len(sim_path) == 0:
        sys.exit('no ncfiles found.')
    elif args.debug:
        print '\tnc files found.'

    # look for save directory
    # does the savepath exist?
    if args.noplot:
        plot=False
    else:
        plot=True

    # define savepath
    savepath = args.savepath
    if args.savepath:
        if savepath[-1] != '/':
            savepath = savepath + '/'
        if debug:
            print 'savepath selected: ', savepath
            print 'looking for save directory...'
        if not osp.exists(savepath):
            if debug:
                print 'directory not found.'
                print 'creating directories...'
                print 'directory {} successfully created.'.format(savepath)
            os.makedirs(savepath)
        elif not osp.isdir(savepath):
            sys.exit('{} is not a directory.'.format(savepath))

    return args.loc[0], sim_path, obs_dir, matfiles, plot, savepath


def createPlots(speedS, speedO, mean, stdev, bias, uspdS, uspdO, lon, lat, \
        model, ubias, bfric, erru, errv, depth, drifts, debug=False):
    """
    Creates a bunch of plots.
    """
    # estimate cube ratio
    ratio = sps.cbrt(np.mean(np.power(uspdO,3))/np.mean(np.power(uspdS,3)))
    if debug:
        print 'speed ratio is: {}'.format(ratio)

    sns.set(font="serif")
    # create plots
    if debug:
        print 'creating plots...'

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(speedS, ubias, alpha=0.25)
    ax.set_xlabel('Model Speed (m/s)')
    ax.set_ylabel('Bias')
    ax.set_title('Bias vs. Model Speed for BFRIC={}'.format(bfric))

    # determine line of best fit
    if debug:
        print '\tdetermining line of best fit...'

    par = np.polyfit(speedS, ubias, 1)
    m = par[-2]
    b = par [-1]

    variance = np.var(ubias)
    residuals = np.var([(m*xx + b - yy)  for xx,yy \
            in zip(speedS, ubias)])
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
    ax1.set_title('Model and Drifter Speed Comparison for BFRIC={}'.format(bfric))
    coeff = np.polyfit(speedS3, speedO3, 1)
    m = coeff[-2]
    b = coeff[-1]
    if debug:
        print '\tcoeffs for cube plot are: \n\t\tm={}\n\t\tb={}'.format(m,b)
    plt.hold('on')
    ax1.plot(speedS3, m*speedS3+b)
    variance = np.var(speedO3)
    residuals = np.var([(m*xx + b - yy)  for xx,yy \
            in zip(speedS3, speedO3)])
    Rsqr = np.round(1-residuals/variance, decimals=5)
    if debug:
        print '\tR^2 for cube plot is {}'.format(Rsqr)

    # plot bias v drifter
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    if debug:
        print '\tnum drift is ', np.arange(1,num_drift+1).size
        print '\tlen of mean is ', len(mean)

    # hacky fix for differing lengths
    try:
        ax2.plot(np.arange(1,num_drift+1), mean, 'go')
    except:
        ax2.plot(np.arange(1,num_drift), mean, 'go')
    ax2.set_ylabel('Bias')
    ax2.set_xlabel('Drifter')
    ax2.set_title('Individual Bias vs Drifter for BFRIC={}'.format(bfric))
    plt.hold('on')
    ax2.axhline(y=np.mean(bias), linewidth=2)
    plt.grid(True)

    # spatial plot of ratios...
    if debug:
        print '\tstarting spatial plot...'
    glon = model.Grid.lon
    glat = model.Grid.lat
    if debug:
        print '\tcomputing bounding boxes...'
    bounds = [np.min(glon), np.max(glon), np.min(glat), np.max(glat)]


    uspeedO3 = np.power(uspdO, 3)
    uspeedS3 = np.power(uspdS, 3)
    if debug:
        print '\tcomputing colors...'
    var3 = np.divide(uspeedO3, uspeedS3)
    # color = np.subtract(var3, np.min(var3)) / (np.max(var3) - np.max(var3))

    if debug:
        print '\tcreating map...'
        print '\tcreating scatter plot...'
	f=plt.figure()
	ax = f.add_axes([.125,.1,.775,.8])
	ax.triplot(glon, glat, model.Grid.trinodes, zorder=10, lw=10)
	clim=np.percentile(var3,[5,95])
	cb = ax.scatter(lon, lat, c=var3, s=10, edgecolor='None', \
            vmin=clim[0], vmax=clim[1], zorder=20)
    plt.colorbar(cb)
    if debug:
        print '\tcreating color bar...'


if __name__ == '__main__':

    # parse the command line args and identify parameters
    print '\nparsing command line options...'
    args = parseArgs()
    debug = args.debug

    loc, sim_path, obs_dir, obs_files, plot, savepath = setOptions(args)

    if debug:
        print '\n--parameters selected--'
        print 'location: ', loc, '\nsim_path: ', sim_path, \
                '\nobs_path: ', obs_dir

    # initialize cumulative data arrays
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
    all_lon0 = []
    all_lat0 = []
    all_ubias = []
    all_depth = []
    all_erru = []
    all_errv = []
    all_err_mag = []
    all_err_dir = []
    num_drift = 0

    for dir_name in sim_path:
        print '\nloading fvcom object...'
        ncfile = FVCOM(dir_name, debug=False)
        if debug:
            print 'ncfile for {} loaded.'.format(dir_name)

        # find the relevant time window to work in
        mTimes = ncfile.Variables.matlabTime[:]
        mStart, mEnd = float(mTimes[0]), float(mTimes[-1])

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

        drift, mean, std, speedO, speedS, uspdO, uspdS, bias, lat, lon, \
            lon0, lat0, ubias, depth, erru, errv, err_mag, err_dir \
            = calculateBias(ncfile, files, loc, debug=debug)

        if debug:
            print 'adding to cumulative data...'

        drifters[dir_name] = drift
        all_mean.extend(mean)
        all_std.extend(std)
        all_speedS.extend(speedS)
        all_speedO.extend(speedO)
        all_uspeedS.extend(uspdS)
        all_uspeedO.extend(uspdO)
        all_bias.extend(bias)
        all_ubias.extend(ubias)
        all_lat.extend(lat)
        all_lon.extend(lon)
        all_lon0.extend(lon0)
        all_lat0.extend(lat0)
        all_depth.extend(depth)
        all_erru.extend(erru)
        all_errv.extend(errv)
        all_err_mag.extend(err_mag)
        all_err_dir.extend(err_dir)
        num_drift = num_drift + len(drift)

        if debug:
            print 'continuing loop...'

    if debug:
        print 'creating as numpy arrays...'
    speedS = np.asarray(all_speedS)
    speedO = np.asarray(all_speedO)
    mean = np.asarray(all_mean)
    stdev = np.asarray(all_std)
    bias = np.asarray(all_bias)
    ubias = np.asarray(all_ubias)
    uspdS = np.asarray(all_uspeedS)
    uspdO = np.asarray(all_uspeedO)
    lat = np.asarray(all_lat)
    lon = np.asarray(all_lon)
    depth = np.array(all_depth)
    erru = np.asarray(all_erru)
    errv = np.asarray(all_errv)
    err_mag = np.asarray(all_err_mag)
    err_dir = np.asarray(all_err_dir)

    # create plots and do stats
    # createPlots(speedS, speedO, mean, stdev, bias, uspdS, uspdO, lat, lon, \
    #        ncfile, ubias, args.bfric, erru, errv, depth, drifters, debug=debug)
    # plt.show()

    # write init loc data to text file
    if args.write:
        if debug:
            print 'recording initial positions...'

        with open('init_locs_'+loc+'.dat', 'w') as f:
            for lon, lat in zip(all_lon0, all_lat0):
                f.write(str(lon) + ' ' + str(lat) + '\n')

    if debug:
        print '...all done!'
