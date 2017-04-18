#!/usr/env/python2.7
#! encoding: utf-8

"""
Filename = compareSpeeds.py
Author = Kody Crowell
Version = 1.4

The program creates a series of plots containing a drifter trajectory superposed
onto a time-averaged velocity norm of the flood/ebb tide, and a comparative speed-
time graph for the corresponding observed and modelled data. The program is
called with a specific location tag and FVCOM simulation, in which case the
program will look for all drifter files on the cluster workspace where the drifter
runs are located and plot the ones that fall within the time window. Otherwise, a
specific file or several files can be given. If the tag show is set to True, then
the program will show the plot(s) WITHOUT saving them. An alternate directory can
be defined as the savepath. Debug increase output verbosity.

Usage:
$ python compareSpeeds.py --loc/-l LOCATION_TAG \
         --dir/-D FVCOM_DIRECTORY_NAME_AS_DATE  [--show/-s] \
         [--savepath/-sp SAVEPATH] [--drifter/-d DRIFT_FILE(S)] \
         [--debug/--verbose/-v] [--version/-V] [--ratio]

Examples:
$ python compareSpeeds.py --help
$ python compareSpeeds.py --loc GP --savepath '/array/home/119865c/karsten/data/' --dir '2013_Aug_01_3D'
$ python compareSpeeds.py -s --loc 'GP' --dir '2013_Aug_01_3D' --debug
$ python compareSpeeds.py -l 'GP' -D '2013_Aug_01_3D' -d 'GP_F_20130801_H_009_SE15.mat'

Note that a region tag must be given if a sim_dir is given, and a sim_dir must be
given if a drift_file is given. Both the location and FVCOM directory flags are
required to run the program.
"""

# Library Imports
import sys, os
import argparse as arp
import numpy as np
import scipy.io as sio
import os.path as osp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as tic
import seaborn as sns
from pyseidon import *

# local imports
from createColorMap import createColorMap
from drifterUtils import *
from drifterPlotUtils import *


PATH_TO_SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
PATH_TO_OBS="/EcoII/acadia_uni/workspace/observed/"
SAVEPATH="/EcoII/acadia_uni/projects/drifters/plots/"
GRID='dngridCSR'


def createPlots(ncfile, files, loc, savepath, sim_name, bfric, tight=False, \
                                         ratio=1.0, debug=False, plot=False):
    """
    Compiles necessary data, creates plots and saves / shows them all. The plot is
    a spatially varying color map of speed with a drifter trajectory superposed and
    a comparitive speed-time plot for each drifter in the given directory.

    input:
        - ncfile : FVCOM object
        - files : list of matlab filenames in directory
        - loc : location tag
        - savepath : savepath to be used (if plot=False)
        - sim_name : name of simulation directory used
        - bfric : bottom friction value
        - plot : True if plot is shown and not saved
        - tight : boolean, True if subdomain region is to be constricted
        - ratio : ratio to adjust model data
    """

    sns.set(font="serif")

    # find the location centre for flood/tide split calculation
    # not yet working...
    if loc == 'GP':
        centre = [-66.33906, 44.26898]
        if tight:
            bounds = [-66.355, -66.31, 44.245, 44.2925]
        else:
            bounds = []
    elif loc == 'DG':
        centre = [-65.76000, 44.67751]
        if tight:
            bounds = [-65.775, -65.77, 44.665, 44.69]
        else:
            bounds = []
    elif loc == 'PP':
        centre = [-66.206924, 44.389368]
        # find out the tightness required for PP
        if tight:
            bounds = [-66.225, -66.195, -44.37, -44.41]
        else:
            bounds = []

    if debug:
        print 'calculating ebb/flood split at centre of location...'
        print 'calculating model velocity norm...'

    fI, eI, _, _ = ncfile.Util2D.ebb_flood_split_at_point(centre[0], centre[1])
    ncfile.Util3D.velo_norm()

    if debug:
        print '{} plot(s) will be created...'.format(len(files))

    if not plot:
        savepath = savepath + 'bfric_' + bfric + '/' + loc + '_' + sim_name
        if ratio != 1.0:
            savepath = savepath + '/with_ratio_{}'.format(str(ratio))

        # creates a subdirectory, so as not to overwrite existing files
        if debug:
            print 'creating new subdirectory...'
        now = datetime.now()
        now = now.strftime("%Y%m%d")
        if not osp.exists(savepath):
            os.makedirs(savepath)
        else:
            savepath = savepath + '/_' + now
            os.makedirs(savepath)
        savepath = savepath + '/'

    for i, fname in enumerate(files, start=1):
        if debug:
            print 'creating drifter object...'
        drift = Drifter(fname, debug=False)
        fname = fname[48:-4]

        # creates drifter object window for flow map
        if debug:
            print 'creating drifter object window...'
        tModel = ncfile.Variables.matlabTime
        tDrift = drift.Variables.matlabTime
        win1 = (np.abs(tModel-tDrift.min())).argmin()
        win2 = (np.abs(tModel-tDrift.max())).argmin()

        tide = str(drift.Data['water_level'].tide)
        # averages velocity norm over flood or ebb cycle within drifter window
        if tide == 'flood':
            tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
        elif tide == 'ebb':
            tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)

        # create spatially varying color map of mean velocity norm
        if debug:
            print 'preparing to create colormap...'
        fig = createColorMap(ncfile, tideNorm[0,:], mesh=False, bounds=bounds, \
                title='Trajectory for ' + fname[:-4], debug=debug, \
                label='Mean Velocity Norm during '+tide.capitalize()+' Tide (m/s)')
        # create title
        # fig.suptitle('Data from ' + fname[:-4], fontsize=14)

        x = drift.Variables.lon
        y = drift.Variables.lat
        u = drift.Variables.u
        v = drift.Variables.v

        if debug:
            print 'creating scatter plot...'
        plt.scatter(x,y)

        if debug:
            print 'preparing to plot time series...'

        # create validation structure
        if debug:
            print 'creating validation object...'

        try:
            valid = Validation(drift, ncfile, flow='sf', debug=False)
        except IndexError:
            print 'cannot create validation object for drifter %i.' % i
            continue

        # calculate speed from interpolated and observed date
        mTimes = valid.Variables.struct['mod_time']
        oU = valid.Variables.struct['obs_timeseries']['u']
        oV = valid.Variables.struct['obs_timeseries']['v']
        mU = valid.Variables.struct['mod_timeseries']['u']
        mV = valid.Variables.struct['mod_timeseries']['v']

        if debug:
            print '\tcalculating speeds...'
        speedS = np.asarray(np.sqrt(mU**2 + mV**2))
        speedO = np.asarray(np.sqrt(oU**2 + oV**2))

        # ratio addition
        if debug:
            print '\tadding ratio adjustments...'
        speedS = speedS * ratio
        datetimes = np.asarray([dn2dt(time) for time in mTimes])

        # For now, separate the two plots.
        # Set BFRIC for now
        # WHAT IS WRONG WITH THIS
        # fig2=plt.figure()
        # result, axis = plotTimeSeries(fig2, np.reshape(np.tile(datetimes,2),\
        #         (2, len(datetimes))), np.vstack((speedS, speedO)), \
        #         loc, label=['Simulated','Observed'], where=111, \
        #         title='Path Speeds for ' + fname[:-4] + ' |  BFRIC=0.015', \
        #         axis_label='Speed (m/s)')

        # if not result:
        #     if debug:
        #         print '...error encountered with drifter {}.'.format(i)
        #         print 'continuing...'
        #     plt.close()
        #     continue

        if plot:
            if debug:
                print 'displaying plot...'
            plt.show()
        else:
            if debug:
                print 'saving plot...'
            fig.savefig(savepath + fname + '_traj.png')
            # result.savefig(savepath + fname + '_speed.png')
            if debug:
                print '...plot saved to: ', savepath+fname

        # clear the figure window
        plt.close()


def parseArgs():
    """
    Parses, identifies and resolves command line arguments.
    """

    parser = arp.ArgumentParser(prog='compareSpeeds.py', description="Plots " \
            + "comparative speed-time graphs and trajectories for given drifter " \
            + "data and model data. Can set option to show or save output.")
    # creates object to be parsed. adds command line optional arguments
    parser.add_argument("--debug", "-v", "--verbose", action="store_true", \
            help="increases output verbosity.")
    parser.add_argument("-V", '--version', action='version', version='v1.4')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--show", "--plot", "-s", action="store_true", \
            help='shows the plot(s) WITHOUT saving them.')
    group.add_argument("--savepath", '-sp', nargs='?', help="defines an  " \
            + "alternate savepath.", metavar='savedir', type=str, default=SAVEPATH)
    # required options bad form, users expect options to be optional
    require = parser.add_argument_group('required flag arguments')
    require.add_argument("--loc", '-l', help="defines the region to operate " \
            + "within.", nargs=1, choices=('GP', 'DG', 'PP'), required=True)
    require.add_argument("--dir", '-D', nargs=1, help="defines an FVCOM " \
            + "directory. Name is in the form YYYY_Mmm_DD_3D", \
            metavar='dirname', type=str, required=True)
    parser.add_argument("--drifter",'-d', nargs='*', help="defines one or " \
            + "more matlab drifter files to be read.", metavar='matfile', \
            type=str, default=False)
    parser.add_argument("--tight", '-t', action='store_true', \
            help="constricts the subdomain region when plotting spatially.")
    parser.add_argument("--bfric", '-B', help='define a bottom friction.', \
            nargs=1, choices=('0.009','0.012','0.015','0.020'), default='0.015', \
            type=str)
    parser.add_argument("--ratio", help='use a ratio adjustment', nargs=1, \
            default=1.0, metavar='ratio', type=float)
    parser._optionals.title = 'optional flag arguments'

    args = parser.parse_args()

    if args.debug:
        # identifies program options
        print '-verbosity turned on.-'

        if args.show:
            print '\tshow set to True, suppressing savepaths...'
        elif args.savepath and not args.show:
            print '\tsavepath identified...'
        if args.loc:
            print '\tlocation tag set to {}...'.format(args.loc)
        if args.dir:
            print '\tfvcom directory set to {}...'.format(args.dir)
        print '\tbottom friction is {}...'.format(args.bfric)

    if not args.loc:
        sys.exit('a location tag is needed. type --help for more info.')

    if args.drifter:
        if args.debug:
            print '\tdrifter files specified...'
        if not args.loc or not args.dir:
            sys.exit('a location tag and fvcom directory is required. ' \
                    + 'type --help for more info.')

    # tags to determine what sort of arguments are returned
    if args.drifter:
        tag = 'yes'
    else:
        tag = 'no'

    if args.ratio:
        print '\tratio selected is {}'.format(args.ratio)

    return args, tag


def setOptions(args):
    """
    Interprets the options returned by argument parser.
    """
    obs_dir = PATH_TO_OBS + args.loc[0] + '/Drifter/'

    if args.drifter:
        # if drifter files are given...
        # use specified drifter files
        if args.debug:
            print 'looking for drifter directory...'
            print '\tgathering all files...'

        matfiles = [obs_dir + file for file in args.drifter]

        for file in matfiles:
            if not osp.exists(file) or not osp.isfile(file):
                sys.exit('problem loading matlab drifter file {}'.format(file))

        if args.debug:
            print '-drifter file(s) successfully identified-'

    elif args.dir:
        # if just an fvcom directory is given...
        # gather all drifter files in region
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
        print 'looking for fvcom directory...'

    if args.bfric:
        path2sim = PATH_TO_SIM + 'BFRIC_' + args.bfric[0] + '/'

    # hacky fix to an odd bug
    if args.dir[0][-1] == ':':
        args.dir[0] = args.dir[0][:-1]

    # locate given fvcom file
    sim_path = path2sim + args.loc[0] + '/' + args.dir[0]

    if not osp.exists(sim_path):
        sys.exit('the directory {} could not be located.'.format(sim_path))
    elif not osp.isdir(sim_path):
        sys.exit('{} is not a directory.'.format(sim_path))
    elif args.debug:
        print '\tfvcom directory found. \n\tloading nc file...'

    sim_path += '/output/subdomain_' + args.loc[0] + '1_0001.nc'
    if not osp.exists(sim_path) or not osp.isfile(sim_path):
         sys.exit('fvcom file not in directory.')
    elif args.debug:
         print '\tfvcom file successfully located.'

    return args.loc[0], sim_path, args.dir[0], obs_dir, matfiles, args.tight


if __name__ == '__main__':

    # parse the command line args and identify parameters
    print '\nparsing command line options...'
    args, opt = parseArgs()
    debug = args.debug

    loc, sim_path, args_dir, obs_dir, obs_files, tight = setOptions(args)

    if debug:
        print '\n--parameters selected--'
        print 'location: ', loc, '\nsim_path: ', sim_path, \
                '\nobs_path: ', obs_dir, '\nmatfiles selected? ', opt

        print '\nloading fvcom object...'

    ncfile = FVCOM(sim_path, debug=False)

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
        sys.exit('drifter files given do not fall within model runtime window.')

    # check if savepath exists, else use default
    savepath = args.savepath

    if not args.show:
        plot=False

        if savepath[-1] != '/':
            savepath += '/'
        if debug:
            print 'savepath selected: ', savepath
            print 'looking for save directory...'

        if not osp.exists(savepath):
            if debug:
                print 'directory not found.'
                print 'creating directories...'
                print 'directory {} succcessfully created.'.format(savepath)
            os.makedirs(savepath)
        elif not osp.isdir(savepath):
            sys.exit('{} is not a directory.'.format(savepath))

    else:
        plot=True

    createPlots(ncfile, files, loc, savepath, args_dir, args.bfric[0], \
            debug=debug, \
            plot=plot, \
            tight=tight, \
            ratio=args.ratio)

    plt.close()
    if debug:
        print '...all done!'
