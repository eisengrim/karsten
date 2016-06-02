#!usr/bin/python2.7
# encoding: utf-8

import sys
import os
from pyseidon import *
import numpy as np
import matplotlib.pyplot as plt
from interpolation_utils import *
import scipy.io as sio
from datetime import datetime, timedelta
import h5py
import netCDF4 as nc

# local imports
from drifterUtils import *


PATH_TO_SIM_FILE="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/" \
   + "drifter_runs/BFRIC_0.015/GP/2013_Aug_08_3D/output/subdomain_GP1_0001.nc"
PATH_TO_OBS_DIR="/EcoII/acadia_uni/workspace/observed/GP/Drifter/"
LOCATION = "GP"
SAVEPATH='/EcoII/acadia_uni/projects/drifters/plots/'

"""
The program plots drifter trajectories taken from a drifter class onto a
spatially varying map of the flow speed (specifically, the flood or ebb
velocity norm). The program relies on a location tag and works in that domain.

Code adopted from Jon Smith's timespan.py.

Program called via command line with optional flags as:
  $ python spatialPlot.py <path_to_fvcom_file> <path_to_obs_dir> <location_tag> -debug
"""


def plotTrajectories(model, dir, matfiles, loc, debug=False):
    """
    Plots the trajectories of the given drifter files against a color map of averaged tidal flow.
    """
    # finds the indices where flood / ebb occurs; computes velocity norm
    fI, eI, pa, pav = model.Util2D.ebb_flood_split_at_point(loc[0], loc[1])
    model.Util3D.velo_norm()

    first = True
    for name in matfiles:
        drift = Drifter(dir+name, debug=debug)
        if debug:
            print 'creating drifter object...'

        if first:
            tide = str(drift.Data['water_level'].tide)
            # averages velocity norm over flood or ebb cycle
            if tide == 'flood':
                tideNorm = np.mean(model.Variables.velo_norm[fI,:,:],0)
            elif tide == 'ebb':
                tideNorm = np.mean(model.Variables.velo_norm[eI,:,:],0)
            # creates spatially varying color map of mean velocity norm
            model.Plots.colormap_var(tideNorm[0,:], mesh=False, hold=True, title='Depth-Averaged Velocity Norm During '+tide.capitalize()+' Tide')
            plt.hold('on')
            first = False

        if not str(drift.Data['water_level'].tide) == tide:
            continue

        x = drift.Variables.lon
        y = drift.Variables.lat
        u = drift.Variables.u
        v = drift.Variables.v

        if debug:
            print 'creating scatter plot..'
        # plt.quiver(x, y, u, v)
        plt.scatter(x,y)
        plt.hold('on')


def compareSpeeds(model, dir, matfiles, loc, savedir='', debug=False):
    """
    Plots a speed v. time graph for each drifter in the directory and saves it
    to a defined directory.
    """

    for i, name in enumerate(matfiles, start=1):
        drift = Drifter(dir+name, debug=debug)
        valid = Validation(drift, model, flow='sf', debug=debug)
        if debug:
            print 'creating drifter and validation object...'

        mTimes = valid.Variables.struct['mod_time'][:-1]

        oU = valid.Variables.struct['obs_timeseries']['u']
        oV = valid.Variables.struct['obs_timeseries']['v']
        mU = valid.Variables.struct['mod_timeseries']['u']
        mV = valid.Variables.struct['mod_timeseries']['v']

        speedS = np.sqrt(mU**2 + mV**2)
        speedO = np.sqrt(oU**2 + oV**2)

        datetimes = [datetime.fromordinal(int(x)) + timedelta(days=x%1) \
                  - timedelta(days = 366) for x in mTimes]

        plt.plot(datetimes, speedS, 'b', label='Simulated')
        plt.hold('on')
        plt.plot(datetimes, speedO, 'ro', label='Observed')
        plt.suptitle('Observed vs. Simulated Speed vs. Time')
        plt.xlabel('Time (HH:MM:SS)')
        plt.ylabel('Speed (m/s)')
        plt.legend()
        plt.gcf().autofmt_xdate()
        if savedir.strip() and name.strip():
            if os.path.exists(savedir) and os.path.isdir(savedir):
                plt.savefig(savedir+LOCATION+'_drifter_'+i+'_speedvtime.png')


if __name__ == '__main__':

    argc = len(sys.argv)

    if '-debug' in sys.argv:
        print '--debug mode on--'
        debug = True
    else:
        debug = False

    if argc >= 4:
        location = sys.argv[3]
    else:
        location = LOCATION

    if argc >= 3:
        simFile = sys.argv[1]
        obsDir = sys.argv[0]
    elif argc < 3:
        if debug:
            print 'missing path to fvcom file and/or drifter directory, using defaults...'
        simFile = PATH_TO_SIM_FILE
        obsDir = PATH_TO_OBS_DIR

    if os.path.exists(simFile) and os.path.isfile(simFile):
        print 'fvcom file succesfully found.'
    else:
        sys.exit('fvcom file not found.')

    if os.path.exists(obsDir) and os.path.isdir(obsDir):
        print 'drifter directory successfully found.'
        print '\tloading all files...'
        dirs = os.listdir(obsDir)

        if len(dirs) == 0:
            sys.exit('no files found in drifter directory.')
    else:
        sys.exit('drifter directory not found.')

    # sets location centre for flood/tide split calculation
    if location == 'GP':
        loc = [-66.33906, 44.26898]
    elif location == 'DG':
        loc = [-65.76000, 44.67751]
    elif location == 'PP':
        loc = [-66.20692, 44.38937]
    elif location == 'MP':
        loc = [-64.40725, 45.34758]
    else: sys.exit("not a valid location tag.")

    # finds relevant time window to work in
    ncfile = FVCOM(simFile, debug=debug)
    mTimes = ncfile.Variables.matlabTime[:]
    mStart, mEnd = float(mTimes[0]), float(mTimes[-1])
    if debug:
        print 'model time is from {} to {}.'.format(mStart, mEnd)

    # finds drifter files in fvcom runtime window
    files = []
    for matfile in dirs:
        if debug:
            print 'gathering all matlab drifter files in model run time...'

        dStart, dEnd = driftTimes(obsDir + matfile, debug=debug)
        dStart, dEnd = float(dStart), float(dEnd)
        if dStart > mStart and mEnd > dEnd:
            files.append(matfile)

    plotTrajectories(ncfile, obsDir, files, loc, debug=debug)
    compareSpeeds(ncfile, obsDir, files, loc, debug=debug, savedir=SAVEPATH)
    plt.show()
