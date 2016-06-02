#!/usr/env/python2.7
#! encoding: utf-8

"""
Filename = classifyDrifters.py
Author = Kody Crowell
Version = 1.0

Given a directory containing Drifter MATLAB files, and a single FVCOM ncfile,
this program classifies drifters based on the depth of their trajectory, the
number of data points, the variability of the data, whether it is during flood
or ebb tide, and where in the tidal cycle the drifter was active.

The drifter directory is defined using the region tag 'GP', 'DG' or 'PP'.
The bottom friction is defined using 0.0XX, where XX is a number.
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
import seaborn as sns
from pyseidon import *
import sklearn.preprocessing as skp

PATH_TO_SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
PATH_TO_OBS="/EcoII/acadia_uni/workspace/observed/"
GRID='dngridCSR'

# local imports
from drifterUtils import dn2dt, driftTimes


def parseArgs():
    """
    Parses, identifies and resolves command line arguments.
    """

    parser = arp.ArgumentParser(prog='classifyDrifters.py', description="Opens " \
            + "an FVCOM file and multiple drifter files and classifies them.")
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
    parser.add_argument("--bfric", '-B', help="select a bottom " \
            + "friction.", nargs=1, choices=('0.009','0.012','0.015','0.020'), \
            default='0.015', type=str)
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

    if args.write:
        # add write for plots...
        print '\tinitial locations to be recorded...'

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

    if args.bfric:
        path2sim = PATH_TO_SIM + 'BFRIC_' + args.bfric + '/'

    # locate given fvcom file
    if args.dir:
        sim_path = [path2sim + args.loc[0] + '/' + args.dir[0]]
        if not osp.exists(sim_path[0]) or not osp.isdir(sim_path[0]):
            sys.exit('the directory {} could not be located.'.format(sim_path))
        elif args.debug:
            print '\tfvcom directory found. \n\tloading nc file...'

        sim_path[0] = sim_path[0]+'/output/subdomain_'+args.loc[0]+'1_0001.nc'
        if not osp.exists(sim_path[0]) or not osp.isfile(sim_path[0]):
            sys.exit('fvcom file not in directory.')
        elif args.debug:
            print '\tfvcom file successfully located.'
    else:
        dirs = os.listdir(path2sim + args.loc[0] + '/')
        sim_path = [path2sim + args.loc[0] + '/' + file + \
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


def classifyDrifts(ncfile, matfiles, loc, data, debug=False):
    """
    Classifies drifters based on several criteria.
    input:
    - ncfile : FVCOM object
    - files : list of matlab filenames in directory
    - loc : location tag
    - data : dictionary of data for drifters
    """
    pass


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
        data = {'start' : [],
                'end' : [],
                'name' : []}

        for matfile in obs_files:
            dStart, dEnd = driftTimes(matfile, debug=debug)
            dStart, dEnd = float(dStart), float(dEnd)
            if dStart > mStart and mEnd > dEnd:
                files.append(matfile)
                if debug:
                    print 'file {} is within runtime window.'.format(matfile)
                data['start'].append(dStart)
                data['end'].append(dEnd)
                data['name'].append(matfile)
        if not files:
            sys.exit('drifters given are not within model runtime window.')

        classifyDrifts(ncfile, files, loc, data, debug=debug)

        if debug:
            print 'adding to cumulative data...'

    if debug:
        print '...all done!'
