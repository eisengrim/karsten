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
import pandas as pd
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
import numpy.linalg as LA
from scipy.stats.stats import pearsonr

# local import
from createColorMap import createColorMap
from drifterUtils import dn2dt, driftTimes
from drifterAnalysisUtils import *

PATH2SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
PATH2OBS="/EcoII/acadia_uni/workspace/observed/"


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
            + "roughness.", nargs=1, choices=('0.009','0.012','0.015','0.020'),\
            default='0.015', required=True, type=str)
    parser._optionals.title = 'optional flag arguments'
    # option to write initial positions
    parser.add_argument("--write", '-w', help='records initial positions of ' \
            + 'drifters to .dat file.', action="store_true")
    parser.add_argument('--tide', '-t', help='adds tidal option', nargs=1, \
            choices=('ebb', 'flood', None), default=None)

    args = parser.parse_args()

    if args.debug:
        # identifies program options
        print '-verbosity turned on.-'

        if args.loc:
            print '\tlocation tag set to {}...'.format(args.loc)

    if not args.loc:
        sys.exit('a location tag is needed. type --help for more info.')

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

    if args.tide:
        tide = args.tide[0]
    else:
        tide = None

    return args.loc[0], sim_path, obs_dir, matfiles, tide


if __name__ == '__main__':

    # use seaborn without Arial font installed
    sns.set(font="serif")
    # parse the command line args and identify parameters
    print '\nparsing command line options...'
    args = parseArgs()
    debug = args.debug

    loc, sim_path, obs_dir, obs_files, tide = setOptions(args)

    if debug:
        print '\n--parameters selected--'
        print 'location: ', loc, '\nsim_path: ', sim_path, \
                '\nobs_path: ', obs_dir, '\ntide: ', str(tide)

    # initialize cumulative data arrays
    drifters = {}
    all_date = []
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
    all_avg_depth = []
    all_erru = []
    all_errv = []
    all_err_mag = []
    all_err_dir = []
    num_drift = 0

    for dir_name in sim_path:
        if debug:
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
            lon0, lat0, ubias, depth, erru, errv, err_mag, err_dir, avg_d \
            = calculateBias(ncfile, files, loc, date=dir_name[82:96], tide_opt=tide, debug=debug)

        if debug:
            print 'adding to cumulative data...'

        # extracts the name of the directory without including the whole path
        drifters[dir_name[82:96]] = drift
        all_date.extend(dir_name[82:96])
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
        all_avg_depth.extend(avg_d)
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
    avg_depth = np.asarray(all_avg_depth)

    print '\n----returned cumulative data----'
    print 'speedS, \nspeedO, \nmean, \nstdev, \nbias, \nubias, \nuspdO, ' \
           + '\nuspdS, \nlat, \nlon, \ndepth, \nerru, \nerrv, \nerr_mag, ' \
           + '\nerr_dir'
    print 'all individul data is encased within the dictionary \'drifters\''
    print 'only drifters for ' + str(tide) + ' are given.'

    # write init loc data to text file
    if args.write:
        if debug:
            print 'recording initial positions...'

        with open('init_locs_'+loc+'.dat', 'w') as f:
            for lon, lat in zip(all_lon0, all_lat0):
                f.write(str(lon) + ' ' + str(lat) + '\n')

    if debug:
        print '...all done!'
