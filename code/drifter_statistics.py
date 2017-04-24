#! /usr/env/python2.7

# library imports
import sys, os
import os.path as osp
import numpy as np
import scipy as sp
import scipy.special as sps
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as tic
import seaborn as sns
import pandas as pd

# pyseidon
from pyseidon_dvt import *

# local import
from color_map import createColorMap
from utils import *
from plot_utils import *
from location_bounds import get_bounds

def calculate_stats(ncfile, fname, loc, date, tide_opt=None, plot=False,
        tight=False, outpath=None, multi=False, tt_window=None):
    """
    Given a Drifter object from PySeidon, statistics are computed based on
    the closest spatial and temporal points to an FVCOM object.

    Takes either two filename strings or two PySeidon objects

    If more than one drifter is in a file, deals with it as best it can...
    Enter date variable as YYYYMMDD.
    """
    # checks if filename or object is passed, open if a filename is passed
    if type(ncfile) == str:
        model = FVCOM(ncfile, debug=False)
    else:
        model = ncfile

    if type(fname) == str:
        drift = Drifter(fname, debug=False)
    else:
        drift = fname

    # get model lon and lats
    mlon = model.Grid.lon[:]
    mlat = model.Grid.lat[:]
    mlonc = model.Grid.lonc[:]
    mlatc = model.Grid.latc[:]

    # get times and datetimes
    mTimes = model.Variables.matlabTime[:]
    oTimes = drift.Variables.matlabTime[:]
    odatetimes = np.asarray([dn2dt(time) for time in oTimes])

    # find the relevant time window to work in
    if not multi:
        mStart, mEnd = float(mTimes[0]), float(mTimes[-1])
        oStart, oEnd = float(oTimes[0]), float(oTimes[-1])
        if oStart < mStart or mEnd < oEnd:
            sys.exit('drifter not within model runtime window.')
        idx = np.arange(0, len(oTimes))
    else:
        # sample the indices within the model time window
        # is this better than what is done previously? yes -- deals with nonconsecutive data
        dates = np.asarray([k.strftime("%Y%m%d") for k in odatetimes])
        idx = np.squeeze(np.where(dates == date))
        # if oStart < mStart or mEnd < oEnd:
        #     start = np.where(oTimes > mStart)[0][0]
        #     end = np.where(oTimes < mEnd)[0][-1]
        #     idx = np.arange(start, end)
        # else:
        #     idx = np.arange(0, len(oTimes))
        # this finds the indices to be used for searching beta and other constant vars
        # old treatment
        try:
            look = np.squeeze(np.where([date in str(k) for k in drift.Data['fn']]))
        # new treatment
        except:
            look = np.arange(0, len(drift.Data['Tr']))

    # get beta and alpha
    # something messed up with extracting alpha from multi run in old treatment
    if not multi:
        # old treatment
        try:
            beta = drift.Data['water_level'].beta
            alpha = drift.Data['velocity'].alpha[idx]
        # new treatment
        except:
            alpha = drift.Data['tide_time'][idx]
            beta = drift.Data['Tr'][idx]
    else:
        # old treatment
        try:
            alpha = drift.Data['velocity'].alpha[idx]
            beta = drift.Data['water_level'].beta[look]
        # new treatment
        except:
            alpha = drift.Data['tide_time'][idx]
            beta = drift.Data['Tr'][look]

    # tests tide condition -- doesn't work for multiple drifters in one file
    # if not multi:
        # old treatment
        # try:
        #     tide = str(drift.Data['water_level'].tide)
        # new treatment
        # except:
        #     tide = np.empty(len(alpha), dtype="S6")
        #     tide[alpha > 0.5] = "flood"
        #     tide[alpha < 0.5] = "ebb"
    # else:
        # old treatment
        # try:
        #     tide = [str(k) for k in drift.Data['water_level'].tide[look]]
        # new treatment
        # except:
    tide = np.empty(len(alpha), dtype="S6")
    tide[alpha < 0.5] = "flood"
    tide[alpha > 0.5] = "ebb"
    tide = np.asarray(tide)

    # overrides tt_window if selected
    if tide_opt is "flood":
        tt_window = [0, 0.5]
    elif tide_opt is "ebb":
        tt_window = [0.5, 1]

    # subset based on tide
    # rather than by tide, may subset using a specified alpha window
    if tt_window is not None:
        a1 = np.squeeze(np.argwhere(alpha >= tt_window[0]))
        a2 = np.squeeze(np.argwhere(alpha <= tt_window[1]))
        tide_idx = np.intersect1d(a1, a2)

        if tide_idx == []:
            return False

    oTimes = oTimes[idx]
    odatetimes = odatetimes[idx]
    # finds closest points
    idx_s = closest_dist(drift.Variables.lon[idx], drift.Variables.lat[idx], \
                model.Grid.lonc, model.Grid.latc)
    idx_t = closest_vals(drift.Variables.matlabTime[idx], model.Variables.matlabTime)

    # collect indexed variables
    mTimes = mTimes[idx_t]
    datetimes = np.asarray([dn2dt(time) for time in mTimes])
    oU = drift.Variables.u[idx]
    oV = drift.Variables.v[idx]
    mU = model.Variables.u[idx_t, 0, idx_s]
    mV = model.Variables.v[idx_t, 0, idx_s]
    olon = drift.Variables.lon[idx]
    olat = drift.Variables.lat[idx]

    # calculate relevant statistics
    speedS = np.sqrt(mU*mU + mV*mV)
    speedO = np.sqrt(oU*oU + oV*oV)
    diffs = np.subtract(speedS, speedO)
    erru = np.subtract(mU, oU)
    errv = np.subtract(mV, oV)
    err_mag = np.sqrt(erru**2 + errv**2)
    err_dir = np.arctan(np.divide(errv, erru))/np.pi * 180
    rmse = np.sqrt(np.mean((speedS-speedO)*(speedS-speedO)))

    nrmse = rmse / np.mean(speedS)
    rbias = np.mean(diffs)
    rerru = erru / mU
    rerrv = errv / mV

    # how to deal with beta
    if tt_window is not None:
        mTimes = mTimes[tide_idx]
        datetimes = datetimes[tide_idx]
        oTimes = oTimes[tide_idx]
        odatetimes = odatetimes[tide_idx]
        oU = oU[tide_idx]
        oV = oV[tide_idx]
        mU = mU[tide_idx]
        mV = mV[tide_idx]
        olon = olon[tide_idx]
        olat = olat[tide_idx]
        alpha = alpha[tide_idx]
        tide = tide[tide_idx]

    # new date for printing
    date = datetime(int(date[0:4]), int(date[4:6]), int(date[6:8])).strftime("%Y-%m-%d")

   # call plotting functions!
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

        # if wanting to use untouched obs datetimes, use vstack((odatetimes,datetimes))
        fig2 = plotTimeSeries(np.reshape(np.tile(datetimes,2), \
            (2, len(datetimes))), np.vstack((speedO, speedS)), \
            loc, label=['Observed', 'Simulated'], where=111, \
            title=loc + ' Drifter Speeds for ' + date,  \
            ylab='Speed (m/s)', style="style_drift.json", lines=['o','^'])

        fig3 = spatialError(mlon, mlat, olon, olat, speedO, speedS, \
                model.Grid.trinodes, label='Signed Error', error='sgn',\
                title='Spatially-distributed '+loc+' Signed Error for '+date)

        if not outpath:
            plt.show()
        else:
            if not osp.isdir(outpath):
                sys.exit("outpath does not exist")
            fig1.savefig(outpath + "/_" + loc + "_" + date + "_tr.png", bbox_inches='tight')
            fig2.savefig(outpath + "/_" + loc + "_" + date + "_ts.png", bbox_inches='tight')
            # hacky fix to access plotting function for PolyCollection/AxesSubplot object
            ax = fig3.get_axes()
            ax.figure.savefig(outpath + "/_" + loc + "_" + date + "_er.png", bbox_inches='tight')

        plt.close()

    return (mTimes, oTimes, oU, oV, mU, mV, speedS, speedO, olon, olat, \
            mlon, mlat, alpha, beta, tide, datetimes, odatetimes, nrmse, \
            rbias, rerru, rerrv, diffs)


def call_many_drifts(drift_dir, model, loc, date, tide_opt=None, plot=False, \
        tight=False, outpath=None, export=False, tt_window=None):
    """
    To operate around multi=True, call this function.
    """
    drifters = {}
    mt = []
    ot = []
    osd = []
    msd = []
    olon = []
    olat = []
    tt = []
    tr = []
    tide = []
    ou = []
    ov = []
    mu = []
    mv = []
    dt = []
    odt = []
    nrmse = []
    rmsd = []
    rerru = []
    rerrv = []
    bias = []

    mlon = model.Grid.lon[:]
    mlat = model.Grid.lat[:]
    mlonc = model.Grid.lonc[:]
    mlatc = model.Grid.latc[:]

    # get list of all drifter files in directory
    matfiles = [drift_dir + f for f in os.listdir(drift_dir)]

    mTimes = model.Variables.matlabTime[:]
    mStart, mEnd = float(mTimes[0]), float(mTimes[-1])
    mfiles = []
    for f in matfiles:
        dStart, dEnd = drift_times(f)
        dStart, dEnd = float(dStart), float(dEnd)
        if dStart > mStart and mEnd > dEnd:
            mfiles.append(f)

    for i, mfile in enumerate(mfiles, start=1):
        drift = Drifter(mfile, debug=False)

        try:
            t = str(drift.Data['water_level'].tide)
        except:
            t = np.empty(len(drift.Data['tide_time']), dtype="S6")
            t[drift.Data['tide_time'] < 0.5] = "flood"
            t[drift.Data['tide_time'] > 0.5] = "ebb"

        if tide_opt is not None:
            if isinstance(t, (str, unicode)):
                if t != tide_opt:
                    continue
            else:
                if tide_opt == 'flood':
                    if 'ebb' in t:
                        continue
                else:
                    if 'flood' in t:
                        continue

        drifters[mfile] = {}
        mtt, ott, ouu, ovv, muu, mvv, mss, oss, olonn, olatt, _, _, \
                alpha, beta, ti, dtt, odtt, rmse, rmsdd, erru, errv, diffs \
                = calculate_stats(model, drift, loc, date)

        mt = np.hstack((mt, mtt))
        ot = np.hstack((ot, ott))
        ou = np.hstack((ou, ouu))
        ov = np.hstack((ov, ovv))
        mu = np.hstack((mu, muu))
        mv = np.hstack((mv, mvv))
        osd = np.hstack((osd, oss))
        msd = np.hstack((msd, mss))
        olon = np.hstack((olon, olonn))
        olat = np.hstack((olat, olatt))
        tt = np.hstack((tt, alpha))
        tr = np.hstack((tr, beta))
        tide = np.hstack((tide, ti))
        dt = np.hstack((dt, dtt))
        odt = np.hstack((odt, odtt))
        nrmse = np.hstack((nrmse, rmse))
        rmse = np.hstack((rmsd, rmsdd))
        bias = np.hstack((bias, diffs))
        rerru = np.hstack((rerru, erru))
        rerrv = np.hstack((rerrv, errv))

        # compile drifter data for export
        drifters[mfile]["filename"] = str(mfile)
        drifters[mfile]["mtime"] = pd.Series(mtt.byteswap().newbyteorder())
        drifters[mfile]["otime"] = pd.Series(ott.byteswap().newbyteorder())
        drifters[mfile]["mspeed"] = pd.Series(mss.byteswap().newbyteorder())
        drifters[mfile]["ospeed"] = pd.Series(oss.byteswap().newbyteorder())
        drifters[mfile]["mu"] = pd.Series(muu.byteswap().newbyteorder())
        drifters[mfile]["mv"] = pd.Series(mvv.byteswap().newbyteorder())
        drifters[mfile]["ou"] = pd.Series(ouu.byteswap().newbyteorder())
        drifters[mfile]["ov"] = pd.Series(ovv.byteswap().newbyteorder())
        drifters[mfile]["olon"] = pd.Series(olonn.byteswap().newbyteorder())
        drifters[mfile]["olat"] = pd.Series(olatt.byteswap().newbyteorder())
        drifters[mfile]["tt"] = pd.Series(alpha.byteswap().newbyteorder())
        drifters[mfile]["tr"] = pd.Series(beta)
        drifters[mfile]["tide"] = pd.Series(ti.byteswap().newbyteorder())
        drifters[mfile]["dt"] = pd.Series(dtt.byteswap().newbyteorder())
        drifters[mfile]["dtt"] = pd.Series(odtt.byteswap().newbyteorder())
        drifters[mfile]["nrmse"] = pd.Series(rmse)
        drifters[mfile]["msd"] = pd.Series(rmsdd)
        drifters[mfile]["bias"] = pd.Series(diffs.byteswap().newbyteorder())
        drifters[mfile]["rerru"] = pd.Series(erru.byteswap().newbyteorder())
        drifters[mfile]["rerru"] = pd.Series(errv.byteswap().newbyteorder())

    # new date for printing
    date = datetime(int(date[0:4]), int(date[4:6]), int(date[6:8])).strftime("%Y-%m-%d")

    drift = {}
    drift[date] = drifters
    # export data to csv is desired
    if export:
        if not osp.isdir(outpath):
            sys.exit("outpath does not exist")
        ids=[]
        frames=[]
        for id, d in drifters.iteritems():
            ids.append(id)
            frames.append(pd.DataFrame(d))
        data = pd.concat(frames, keys=ids)
        if outpath[-1] is not "/":
            outpath.append("/")
        data.to_csv(outpath+"_drifter_data_{}.csv".format(loc))

    # call plotting functions!
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
                olon, olat, ot, model.Grid.trinodes, vel_norm = vnorm,\
                label = "Mean Velocity Norm (m/s)", bounds=box, \
                title_main = loc + " Trajectory for " + date)

        # if wanting to use untouched obs datetimes, use vstack((odatetimes,datetimes))
        fig2 = plotTimeSeries(np.reshape(np.tile(dt,2),\
            (2, len(dt))), np.vstack((osd, msd)), \
            loc, label=['Observed', 'Simulated'], where=111, \
            title=loc + ' Drifter Speeds for ' + date,  \
            ylab='Speed (m/s)', style="style_drift.json", lines=['o','^'])

        fig3 = spatialError(mlon, mlat, olon, olat, osd, msd, \
                model.Grid.trinodes, label='Relative Signed Error', error='sgn',\
                title='Spatially-Distributed '+loc+' Signed Error for '+date)

        if not outpath:
            plt.show()
        else:
            if not osp.isdir(outpath):
                sys.exit("outpath does not exist")
            fig1.savefig(outpath + "/_" + loc + "_" + date + "_tr.png", bbox_inches='tight')
            fig2.savefig(outpath + "/_" + loc + "_" + date + "_ts.png", bbox_inches='tight')
            # hacky fix to access plotting function for PolyCollection/AxesSubplot object
            ax = fig3.get_axes()
            ax.figure.savefig(outpath + "/_" + loc + "_" + date + "_er.png", bbox_inches='tight')

        plt.close()

    return drift
