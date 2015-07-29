#!/usr/env/python2.7
#! encoding: utf-8

"""
Filename = analyzeDrifters.py
Author = Kody Crowell
Version = 1.0

This code is an amalgamation of previously written code such as drifterBias.py,
spatialPlot.py, compareSpeeds.py and quickPlotDrifter.py. Several functions to
run various stats / plot data are combined and chosen using command line input.
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
from mpl_toolkits.basemap import Basemap
import seaborn as sns
from pyseidon import *
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
        uspeedS = np.asarray(np.sqrt(mU**2 + mV**2))
        uspeedO = np.asarray(np.sqrt(oU**2 + oV**2))

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
        all_sdev.append(sdev_bias)
        all_mean.append(mean_bias)
        all_bias.extend(diffs)
        all_ubias.extend(udiffs)

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
         mod_uspeed, all_bias, o_lat, o_lon, lon0, lat0, all_ubias


# for ncfile:
#     for drifter:
#         extract information
#         compare
#         save plots
#     add to cumulative data
# add to cumulative data
# calculate stats
# plot results
