#!/usr/env/python2.7
#! encoding: utf-8

"""
Filename = findStartingLocations.py
Author = Kody Crowell
Version = 1.0

This code is tasked with finding the initial locations
and times (and, subsequently, time steps in the ncfile)
for drifters within an FVCOM run and saving them.
"""

# library imports
import sys, os
import numpy as np
import scipy as sp
import scipy.io as sio
import os.path as osp
from pyseidon import *
from functionsFvcom import *

PATH2SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
PATH2OBS="/EcoII/acadia_uni/workspace/observed/"


def driftTimes(name, obs_dir):
    """
    Identify the timespans for each drifter file. Adapted from Jon Smith's
    timespan.py.

    input:
        - drifter file name
        - drifter directory
    returns:
        - starting and ending time
    """
    # iterate through the drifter files
    start_time = []
    end_time = []

    try:
        print 'examining drifter file...'
        drft = sio.loadmat(obs_dir+name)
        times = drft['gps_observation'][0][0][0][0]
    except KeyError:
        times = drft['time'][0][0][0][0]

    # grab the times and convert them to strings
    start, end = times[0], times[-1]
    return start, end


def findStartConditions(ncfile, files, odir):
    """
    Finds the starting conditions.

    input:
        - ncfile : FVCOM object
        - files : list of matlab filenames in directory
        - odir : directory of observed files
    """

    print '{} drifters will be analysed...'.format(len(files))

    data_count = 0
    names= []
    lon0 = []
    lat0 = []
    x0 = []
    y0 = []
    mtime0 = []
    firststep = []
    laststep = []

    for i, fname in enumerate(files, start=1):
        print 'working with ' + fname + '...'

        drift = Drifter(odir+fname, debug=False)

        print '\tcompiling data...'
        lon0.append(drift.Variables.lon[0])
        lat0.append(drift.Variables.lat[0])
        mtime0.append(drift.Variables.matlabTime[0])
        firststep.append(int(np.argmin(np.abs(ncfile.Variables.matlabTime - \
                    drift.Variables.matlabTime[0]))))
        laststep.append(int(np.argmin(np.abs(ncfile.Variables.matlabTime - \
                    drift.Variables.matlabTime[-1]))))
        names.append(fname)

        idx = closest_points(drift.Variables.lon[0], drift.Variables.lat[0], \
                            ncfile.Grid.lon, ncfile.Grid.lat)

        x0.append(ncfile.Grid.x[idx])
        y0.append(ncfile.Grid.y[idx])

    return names, lon0, lat0, x0, y0, mtime0, firststep, laststep


if __name__ == '__main__':
    """
    usage:
    findStartingLocations.py [GP, DG, PP] [0.015, 0.012, 0.009] YYYY_Mmm_DD_3D
    """

    loc = str(sys.argv[1])
    bfric = str(sys.argv[2])
    sim_name = str(sys.argv[3])

    obs_dir = PATH2OBS + loc + '/Drifter/'
    dir_name = PATH2SIM + 'BFRIC_' + bfric + '/' + loc + '/' + sim_name + \
            '/output/subdomain_' + loc + '1_0001.nc'

    print '\nloading fvcom object...'
    ncfile = FVCOM(dir_name, debug=False)

    # find the relevant time window to work in
    mTimes = ncfile.Variables.matlabTime[:]
    mStart, mEnd = float(mTimes[0]), float(mTimes[-1])

    # from given drifter files, find files in fvcom runtime window
    print 'gathering all drifter files in model runtime window...'

    obs_files = os.listdir(obs_dir)

    files = []
    for matfile in obs_files:
        dStart, dEnd = driftTimes(matfile, obs_dir)
        dStart, dEnd = float(dStart), float(dEnd)
        if dStart > mStart and mEnd > dEnd:
            files.append(matfile)
            print 'file {} is within runtime window.'.format(matfile)
    if not files:
        sys.exit('drifters given are not within model runtime window.')

    nom, lon0, lat0, x0, y0, mtime0, first, last = findStartConditions(ncfile, \
                files, obs_dir)

    # write init loc data to text file
    print 'recording initial positions...'

    outname = loc + '_' + sim_name

    print len(nom), len(lon0), len(first)

    with open('init_locs_'+outname+'.dat', 'w') as f:
        for name, lon, lat in zip(nom, lon0, lat0):
            f.write(str(name) + ' ' + str(lon) + ' ' + str(lat) + '\n')

    with open('init_time_steps_'+outname+'.dat', 'w') as f:
        for name, idx1, idx2 in zip(nom, first, last):
            f.write(str(name) + ' ' + str(idx1) + ' ' + str(idx2) + '\n')

    np.savetxt('pyticle_info_'+outname+'.txt', \
                np.column_stack((nom,lon0,lat0,first,last)), fmt='%s')

    np.savetxt('porpoise_info_'+outname+'.txt', \
                np.column_stack((nom,x0,y0,first,last)), fmt='%s')

    print '...all done!'
