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

def calculateMKE(drift):
    pass

def calculate_stats(ncfile, fname, loc, date, plot=False,
        tight=False, outpath=None):
    """
    Given a Drifter object from PySeidon, statistics are computed based on
    the closest spatial and temporal points to an FVCOM object.

    inputs
    ncfile  : FVCOM object or a netcdf filename
    fname   : Drifter object or a matlab filename
    loc     : location tag, string, one of "GP", "MP", "DG", "PP"
    date    : date string, as YYYMMDD, ie. "20130801"
    plots   : boolean, generate plots
    tight   : boolean, apply tighting bounding boxes to plots
    outpath : directory to save plots in. if none, displays plots
    multi   : apply separate treatments for multiple drifters in one file/object
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

    # get times and datetimes
    mTimes = model.Variables.matlabTime[:]
    oTimes = drift.Variables.matlabTime[:]
    odatetimes = np.asarray([dn2dt(time) for time in oTimes])

    # find the relevant time window to work in
    mStart, mEnd = float(mTimes[0]), float(mTimes[-1])
    oStart, oEnd = float(oTimes[0]), float(oTimes[-1])
    if oStart < mStart or mEnd < oEnd:
        sys.exit('drifter not within model runtime window.')
    idx = np.arange(0, len(oTimes))

    # get beta and alpha
    # something messed up with extracting alpha from multi run in old treatment
    # old treatment
    try:
        beta = drift.Data['water_level'].beta
        alpha = drift.Data['velocity'].alpha[idx]
    # new treatment
    except:
        alpha = drift.Data['tide_time'][idx]
        beta = drift.Data['Tr'][idx]
        tide = np.empty(len(alpha), dtype="S6")

    # get tide info
    tide[alpha < 0.5] = "flood"
    tide[alpha > 0.5] = "ebb"
    tide = np.asarray(tide)

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
    mlonc = model.Grid.lonc[idx_s]
    mlatc = model.Grid.latc[idx_s]
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

    # new date for printing
    date = datetime(int(date[0:4]), int(date[4:6]), int(date[6:8])).strftime("%Y-%m-%d")

   # call plotting functions!
    if plot:
        # uncomment below if your model does not have velocity norm
        # model.Util3D.velo_norm()
        # vnorm = model.Variables.velo_norm[:]
        vnorm = None

        if tight:
            box = get_bounds(loc)
        else:
            box = []

        # call plotting functions
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
            mlonc, mlatc, alpha, beta, tide, datetimes, odatetimes, nrmse, \
            rbias, rerru, rerrv, diffs)


def call_many_drifts(drift_dir, model, loc, date, tide_opt=None, plot=False, \
        tight=False, outpath=None, export=False, hide=False):
    """
    given many drifters, do a bunch of stats. returns a dictionary with all data
    in the form drift['YYYYMMDD']['drifter file name'] ...

    inputs
    drift_dir   : directory of drifter mat files
    model       : fvcom object or netcdf filename
    loc         : location string tag ("GP", "MP", "DG", "PP")
    date        : date string of the model in the form YYYYMMDD
    tide_opt    : subset drifter data based on whether "flood" or "ebb"
    plot        : boolean, generate plots
    tight       : boolean, constrict plotting regions
    outpath     : save directory. if none, displays plots
    export      : export compiled drifter data to csv in outpath directory
    hide        : hide long/lats in spatial plots
    """
    # initialize empty data arrays
    drifters = {}
    mt = []
    ot = []
    osd = []
    msd = []
    olon = []
    olat = []
    mlonc = []
    mlatc = []
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

    # get list of all drifter files in directory
    matfiles = [drift_dir + f for f in os.listdir(drift_dir)]

    # collect drifter files within fvcom time window
    mTimes = model.Variables.matlabTime[:]
    mStart, mEnd = float(mTimes[0]), float(mTimes[-1])
    mfiles = []
    for f in matfiles:
        dStart, dEnd = drift_times(f)
        dStart, dEnd = float(dStart), float(dEnd)
        if dStart > mStart and mEnd > dEnd:
            mfiles.append(f)

    # loop through drifter files
    for i, mfile in enumerate(mfiles, start=1):
        drift = Drifter(mfile, debug=False)

        # get tide of drifter
        try:
            t = str(drift.Data['water_level'].tide)
        except:
            t = np.empty(len(drift.Data['tide_time']), dtype="S6")
            t[drift.Data['tide_time'] < 0.5] = "flood"
            t[drift.Data['tide_time'] > 0.5] = "ebb"

        # test if drifter is in tide option, if not, skip
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
        # call stats function
        mtt, ott, ouu, ovv, muu, mvv, mss, oss, olonn, olatt, mlon, mlat, \
                alpha, beta, ti, dtt, odtt, rmse, rmsdd, erru, errv, diffs \
                = calculate_stats(model, drift, loc, date)

        # append to existing array for stats
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
        mlonc = np.hstack((mlonc, mlon))
        mlatc = np.hstack((mlatc, mlat))
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
        drifters[mfile]["mlon"] = pd.Series(mlon.byteswap().newbyteorder())
        drifters[mfile]["mlat"] = pd.Series(mlat.byteswap().newbyteorder())
        drifters[mfile]["tt"] = pd.Series(alpha.byteswap().newbyteorder())
        drifters[mfile]["tr"] = pd.Series(beta)
        drifters[mfile]["tide"] = pd.Series(ti.byteswap().newbyteorder())
        drifters[mfile]["dt"] = pd.Series(dtt.byteswap().newbyteorder())
        drifters[mfile]["odt"] = pd.Series(odtt.byteswap().newbyteorder())
        drifters[mfile]["nrmse"] = pd.Series(rmse)
        drifters[mfile]["msd"] = pd.Series(rmsdd)
        drifters[mfile]["bias"] = pd.Series(diffs.byteswap().newbyteorder())
        drifters[mfile]["rerru"] = pd.Series(erru.byteswap().newbyteorder())
        drifters[mfile]["rerru"] = pd.Series(errv.byteswap().newbyteorder())

    # final drifter structure for export
    drift = {}
    drift[date] = drifters

    # new date for printing
    date = datetime(int(date[0:4]), int(date[4:6]), int(date[6:8])).strftime("%Y-%m-%d")

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
                title_main = loc + " Trajectory for " + date, hide=hide)

        # if wanting to use untouched obs datetimes, use vstack((odatetimes,datetimes))
        fig2 = plotTimeSeries(np.reshape(np.tile(dt,2),\
            (2, len(dt))), np.vstack((osd, msd)), \
            loc, label=['Observed', 'Simulated'], where=111, \
            title=loc + ' Drifter Speeds for ' + date,  \
            ylab='Speed (m/s)', style="style_drift.json", lines=['o','^'])

        fig3 = spatialError(mlon, mlat, olon, olat, osd, msd, \
                model.Grid.trinodes, label='Relative Signed Error', error='sgn',\
                title='Spatially-Distributed '+loc+' Signed Error for '+date, hide=hide)

        # MEASURE OF SIMILARITY BETWEEN TRAJECTORY SPEEDS
        # meanm = []
        # meano = []
        # sdm = []
        # sdo = []
        # _, idx_ll, cntl = np.unique(mlonc, return_index=True, return_counts=True)
        # _, idx_tt, cntt = np.unique(mt, return_index=True, return_counts=True)
        # for i, ll in enumerate(idx_ll):
        #     meanm = np.hstack((meanm, np.mean(osd[ll:ll+cnt[i]])))



        # fig4 = plotTimeSeries(np.reshape(np.tile(dt,2),\
        #     (2, len(dt))), np.vstack((osd, msd)), \
        #     loc, label=['Observed', 'Simulated'], where=111, \
        #     title=loc + ' Drifter Speeds for ' + date,  \
        #     ylab='Speed (m/s)', style="style_drift.json", lines=['o','^'])

        # display or save plots to outpath
        if not outpath:
            plt.show()
        else:
            if not osp.isdir(outpath):
                sys.exit("outpath does not exist")
            if tide_opt is None:
                tide_opt = ""
            fig1.savefig(outpath + "/_" + loc + tide_opt + "_" + date + "_tr.png", bbox_inches='tight')
            fig2.savefig(outpath + "/_" + loc + tide_opt + "_" + date + "_ts.png", bbox_inches='tight')
            # hacky fix to access plotting function for PolyCollection/AxesSubplot object
            ax = fig3.get_axes()
            ax.figure.savefig(outpath + "/_" + loc + tide_opt + "_" + date + "_er.png", bbox_inches='tight')

        plt.close()

    return drift
