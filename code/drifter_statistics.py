#! /usr/env/python2.7

# library imports
import sys, os
import numpy as np
import scipy as sp
import scipy.special as sps
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as tic
import seaborn as sns

# pyseidon
from pyseidon_dvt import *

# local import
from color_map import createColorMap
from utils import *
from drifter_plots import *

def calculate_stats(ncfile, fname, loc, date, tide_opt=None, outpath=None):
    """
    Given a Drifter object from PySeidon, statistics are computed based on
    the closest spatial and temporal points to an FVCOM object.
    """

    mlon = ncfile.Grid.lon
    mlat = ncfile.Grid.lat
    mlonc = ncfile.Grid.lonc
    mlatc = ncfile.Grid.latc

    filename = path_leaf(fname)[:-4]

    # create drifter object
    drift = Drifter(fname, debug=False)
    model = FVCOM(ncfile, debug=False)

    # find the relevant time window to work in
    mTimes = model.Variables.matlabTime[:]
    mStart, mEnd = float(mTimes[0]), float(mTimes[-1])

    # print 'model time is from {} to {}.'.format(mStart, mEnd)
    # from given drifter files, find files in fvcom runtime window
    files = []
    for matfile in obs_files:
        dStart, dEnd = driftTimes(matfile, debug=debug)
        dStart, dEnd = float(dStart), float(dEnd)
        if dStart > mStart and mEnd > dEnd:
            files.append(matfile)
    if not files:
        sys.exit('drifters given are not within model runtime window.')

    # tests tide condition
    try:
        tide = str(drift.Data['water_level'].tide)
    except:
        tide = np.empty(len(drift.Data['tide_time']), dtype="S6")
        tide[drift.Data['tide_time'] > 0.5] = "flood"
        tide[drift.Data['tide_time'] < 0.5] = "ebb"

    if tide_opt is not None:
        tide_idx = np.squeeze(np.argwhere(tide == tide_opt))

    idx_s, idx_t = findClosestInds(drift, ncfile)
    mTimes = mTimes[idx_t]
    oTimes = drift.Variables.matlabTime
    oU = drift.Variables.u
    oV = drift.Variables.v
    mU = model.Variables.u[idx_t, 0, idx_s]
    mV = model.Variables.v[idx_t, 0, idx_s]
    olon = drift.Variables.lon
    olat = drift.Variables.lat

    # calculate relevant statistics
    uspeedS = np.sqrt(mU*mU + mV*mV)
    datetimes = np.asarray([dn2dt(time) for time in mTimes])

    uspeedO = np.sqrt(oU*oU + oV*oV)
    print uspeedS.shape, uspeedO.shape

    speedO = uspeedO * np.sign(oV)
    speedS = uspeedS * np.sign(mV)

    diffs = np.subtract(uspeedS, uspeedO)
    udiffs = np.subtract(speedS, speedO)
    erru = np.subtract(mU, oU)
    errv = np.subtract(mV, oV)
    err_mag = np.sqrt(erru**2 + errv**2)
    err_dir = np.arctan(np.divide(errv, erru))/np.pi * 180
    rmse = np.sqrt(np.mean((uspeedS-uspeedO)*(uspeedS-uspeedO)))

    try:
        beta = drift.Data['water_level'].beta
    except:
        beta = drift.Data['Tr']

    try:
        ind = [np.argmin(np.abs(drift.Data['velocity'].vel_time - x)) \
            for x in valid.Variables.struct['mod_time']]
        alpha = drift.Data['velocity'].alpha[ind]
    except:
        alpha = drift.Data['tide_time']

    # if outpath:
    #     ids=[]
    #     frames=[]
    #     for id, d in drift.iteritems():
    #         ids.append(id)
    #         frames.append(pd.DataFrame(d))
    #     data = pd.concat(frames, keys=ids)
    #     if obs_dir[-1] is not "/":
    #         obs_dir.append("/")
    #     data.to_csv(outpath+"_drifter_data_{}.csv".format(loc))


