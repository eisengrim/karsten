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
from plot_utils import *
from location_bounds import get_bounds

def calculate_stats(ncfile, fname, loc, date, tide_opt=None, plot=False,
        tight=False, outpath=None):
    """
    Given a Drifter object from PySeidon, statistics are computed based on
    the closest spatial and temporal points to an FVCOM object.

    Takes either two filename strings or two PySeidon objects
    """

    if type(ncfile) == str:
        model = FVCOM(ncfile, debug=False)
        drift = Drifter(fname, debug=False)
    else:
        model = ncfile
        drift = fname

    mlon = model.Grid.lon[:]
    mlat = model.Grid.lat[:]
    mlonc = model.Grid.lonc[:]
    mlatc = model.Grid.latc[:]

    mTimes = model.Variables.matlabTime[:]
    oTimes = drift.Variables.matlabTime[:]

    # find the relevant time window to work in
    mStart, mEnd = float(mTimes[0]), float(mTimes[-1])
    oStart, oEnd = float(oTimes[0]), float(oTimes[-1])
    if oStart < mStart or mEnd < oEnd:
        sys.exit('drifter not within model runtime window.')

    # tests tide condition
    try:
        tide = str(drift.Data['water_level'].tide)
    except:
        tide = np.empty(len(drift.Data['tide_time']), dtype="S6")
        tide[drift.Data['tide_time'] > 0.5] = "flood"
        tide[drift.Data['tide_time'] < 0.5] = "ebb"

    if tide_opt is not None:
        tide_idx = np.squeeze(np.argwhere(tide == tide_opt))

    # finds closest points
    idx_s = closest_dist(drift.Variables.lon, drift.Variables.lat, \
                model.Grid.lonc, model.Grid.latc)
    idx_t = closest_vals(drift.Variables.matlabTime, model.Variables.matlabTime)

    # collect indexed variables
    mTimes = mTimes[idx_t]
    oU = drift.Variables.u
    oV = drift.Variables.v
    mU = model.Variables.u[idx_t, 0, idx_s]
    mV = model.Variables.v[idx_t, 0, idx_s]

    olon = drift.Variables.lon
    olat = drift.Variables.lat

    # calculate relevant statistics
    datetimes = np.asarray([dn2dt(time) for time in mTimes])
    speedS = np.sqrt(mU*mU + mV*mV)
    speedO = np.sqrt(oU*oU + oV*oV)

    diffs = np.subtract(speedS, speedO)
    erru = np.subtract(mU, oU)
    errv = np.subtract(mV, oV)
    err_mag = np.sqrt(erru**2 + errv**2)
    err_dir = np.arctan(np.divide(errv, erru))/np.pi * 180
    rmse = np.sqrt(np.mean((speedS-speedO)*(speedS-speedO)))

    try:
        beta = drift.Data['water_level'].beta
    except:
        beta = drift.Data['Tr']

    try:
        ind = [np.argmin(np.abs(drift.Data['velocity'].vel_time - x)) \
            for x in mTimes]
        alpha = drift.Data['velocity'].alpha[ind]
    except:
        alpha = drift.Data['tide_time']

    if outpath:
        ids=[]
        frames=[]
        for id, d in drift.iteritems():
            ids.append(id)
            frames.append(pd.DataFrame(d))
        data = pd.concat(frames, keys=ids)
        if obs_dir[-1] is not "/":
            obs_dir.append("/")
        data.to_csv(outpath+"_drifter_data_{}.csv".format(loc))

    if plot:
        #model.Util3D.velo_norm()
        #vnorm = model.Variables.velo_norm[:]
        vnorm = None

        if tight:
            box = get_bounds(loc)
        else:
            box = []

        fig1 = trajectoryPlot(model.Variables.u, model.Variables.v, \
                mlon, mlat, model.Variables.matlabTime, \
                olon, olat, oTimes, model.Grid.trinodes, vel_norm = vnorm,\
                label = "Mean Velocity Norm (m/s)", bounds=box, \
                title_main = loc + " Trajectory for " + date)

        fig2 = plotTimeSeries(np.reshape(np.tile(datetimes,2),\
                (2, len(datetimes))), np.vstack((speedO, speedS)), \
                loc, label=['Observed', 'Simulated'], where=111, \
                title=loc + ' Drifter Speeds for ' + date,  \
                ylab='Speed (m/s)', style="style_drift.json", lines=['o','^'])

        plt.show()

