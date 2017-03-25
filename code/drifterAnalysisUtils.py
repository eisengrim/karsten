#! /usr/env/python2.7

import numpy as np
import pandas as pd
import scipy.io as sio
#from pyseidon import *
from pyseidon_dvt import *
import sys
import scipy.special as sps
# local imports
from drifterUtils import dn2dt, path_leaf, checkIDs

def computeTideNorms(model, pytkl):
    """
    *Run after initializeFVCOM.py*

    **add drifter computation

    model :: FVCOM model
    pytkl :: pytkl object
    """
    model.Util3D.velo_norm()
    tmodel = model.Variables.julianTime
    tpytkl = pytkl.variables['time'][:]
    win1 = (np.abs(tmodel-tpytkl.min())).argmin()
    win2 = (np.abs(tmodel-tpytkl.max())).argmin()
    if win1 == win2:
        tideNorm = np.mean(model.Variables.velo_norm[win1,:,:], 0)
    else:
        tideNorm = np.mean(model.Variables.velo_norm[win1:win2,:,:], 0)

    return tideNorm


def cubeRatio(uspdO, uspdS, debug=False, plot=False):
    """
    Estimates the cubed-ratio and returns it.
    """
    ratio = sps.cbrt(np.mean(np.power(uspdO,3))/np.mean(np.power(uspdS,3)))

    if debug:
        print 'speed ratio is {}'.format(ratio)

    return ratio


def evalWindSpeed(drift, debug=False):
    """
    Compares the direction and magnitude of the wind speed for a particular
    drift with the calculated average error magnitude and direction, using the
    cosine similarity (-1 is exact opposite, 1 is exact same, 0 is orthogonal).
    """
    wind = drift['wind_speed']
    if wind != '0':
        mag = float(''.join([k for k in wind if k.isdigit()]))
        dir_brng = ''.join([k for k in wind if k.isalpha()])

        headings = {'N' : 360.0,
                    'NNE' : 22.5,
                    'NE' : 45.0,
                    'ENE' : 67.5,
                    'E' : 90.0,
                    'ESE' : 112.5,
                    'SE' : 135.0,
                    'SSE' : 160.0,
                    'S' : 180.0,
                    'SSW' : 202.5,
                    'SW' : 225.0,
                    'WSW' : 247.5,
                    'W' : 270.0,
                    'WNW' : 292.5,
                    'NW' : 315.0,
                    'NNW' : 337.5}

        dir = np.deg2rad(headings[dir_brng])
        if debug:
            print 'bearing is {}'.format(dir_brng)
            print 'tide is {}'.format(drift['tide'])

        e_dir = np.deg2rad(np.mean(drift['err_dir']))
        e_mag = np.mean(drift['err_mag'])

        wind_vec = np.array([mag*np.cos(dir), mag*np.sin(dir)])
        err_vec = np.array([e_mag*np.cos(e_dir), e_mag*np.sin(e_dir)])

        similarity = np.dot(wind_vec, err_vec) / (e_mag * mag)
        if debug:
            print 'cosine similarity is {}'.format(similarity)
        return similarity


def compareBFRIC(drift, debug=False):
    """
    For a single drift, a comparative plot is generated with all the model runs
    that have a different bottom friction.
    """
    pass


def drifterParse(ncfile, drift, idx=None, debug=False):
    """
    Does all the laborious compiations and computations for calculateBias.
    """


def calculateBias(ncfile, files, loc, date, tide_opt=None, debug=False,
        self_build=False):
    """
    Compiles necessary data, calculates individual biases and
    takes the mean and standard deviation. Also calculates the
    bias for ALL drifter files. Returns the mean biases, the drifter
    numbers, the model and observed speeds the the std deviations.
    Checks for an array of indices for which unique drifters exist assuming
    the drifters are compiled into one matlab file.

    input:
        - ncfile : FVCOM object
        - files : list of matlab filenames in directory
        - loc : location tag
        - date : simulation date
        - tide_opt : 'ebb' or 'flood' ... returns only drifter data for tide
            This option is not yet available for compiled drifter data
        - sim_name : name of simulation directory used
        - self_build : don't use PySeidon Drifter object
    """

    if debug:
        print '{} drifters will be analysed...'.format(len(files))
    #    print 'calculating depth...'
    # ncfile.Util3D.depth()

    data_count = 0
    drifters = {}
    dates = {}
    all_diff = []
    all_udiff = []
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
    depths = []
    all_erru = []
    all_errv = []
    all_err_mag = []
    all_err_dir = []
    mean_depth = []

    mlon = ncfile.Grid.lon
    mlat = ncfile.Grid.lat
    mlonc = ncfile.Grid.lonc
    mlatc = ncfile.Grid.latc

    for i, fname in enumerate(files, start=1):

        # drifters[fname[48:-4]] = {}
        filename = path_leaf(fname)[:-4]
        drifters[filename] = {}
        if debug:
            print 'creating drifter object {}...'.format(i)
            print fname

        if not self_build:
            drift = Drifter(fname, debug=False)

        check = checkIDs(fname, debug=True)
         # tests tide condition
        try:
            tide = str(drift.Data['water_level'].tide)
        except:
            tide = np.empty(len(drift.Data['tide_time']), dtype="S6")
            tide[drift.Data['tide_time'] > 0.5] = "flood"
            tide[drift.Data['tide_time'] < 0.5] = "ebb"

        if tide_opt is not None:
            if isinstance(tide, (str, unicode)):
                if tide != tide_opt:
                    if debug:
                        print fname + ' not in tide selection'
                    continue
            else:
                tide_idx = np.squeeze(np.argwhere(tide == tide_opt))

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
        oTimes = valid.Variables.struct['obs_time']
        oU = valid.Variables.struct['obs_timeseries']['u']
        oV = valid.Variables.struct['obs_timeseries']['v']
        mU = valid.Variables.struct['mod_timeseries']['u']
        mV = valid.Variables.struct['mod_timeseries']['v']
        olon = valid.Variables.struct['lon']
        olat = valid.Variables.struct['lat']
        # what if len(modtime) =/= len(lon)?
        uspeedS = np.sqrt(mU**2 + mV**2)
        uspeedO = np.sqrt(oU**2 + oV**2)

        # index finder is only in dvt branch
        depth = np.asarray([ncfile.Util3D.depth_at_point(j,k) for j,k in zip(olon,olat)])
        datetimes = np.asarray([dn2dt(time) for time in mTimes])

        if debug:
            print '\tcalculating signed speeds...'
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
        mean_diff = np.mean(diffs)
        sdev_diff = np.std(diffs)
        erru = np.subtract(mU, oU)
        errv = np.subtract(mV, oV)
        err_mag = np.sqrt(erru**2 + errv**2)
        err_dir = np.arctan(np.divide(errv, erru))/np.pi * 180

        # info on # of data this drifter has and other individual data
        drifters[filename]['filename'] = fname
        fname = filename
        drifters[fname]['date'] = date
        drifters[fname]['num'] = i
        drifters[fname]['start'] = data_count
        drifters[fname]['stop'] = data_count + len(olon) - 1
        data_count += len(olon)
        try:
            drifters[fname]['beta'] = drift.Data['water_level'].beta
        except:
            drifters[fname]['beta'] = drift.Data['Tr']
        try:
            ind = [np.argmin(np.abs(drift.Data['velocity'].vel_time - x)) \
                for x in valid.Variables.struct['mod_time']]
            drifters[fname]['alpha'] = pd.Series(drift.Data['velocity'].alpha[ind])
        except:
            drifters[fname]['alpha'] = pd.Series(drift.Data['tide_time'])
        if "ID" in drift.Data.keys():
            drifters[fname]["ID"] = drift.Data["ID"]
        drifters[fname]['tide'] = tide
        drifters[fname]['obs_speed'] = pd.Series(speedO)
        drifters[fname]['sim_speed'] = pd.Series(speedS)
        drifters[fname]['diff'] = pd.Series(diffs)
        drifters[fname]['mean_diff'] = pd.Series(mean_diff)
        drifters[fname]['sdev_diff'] = pd.Series(sdev_diff)
        drifters[fname]['err_u'] = pd.Series(erru)
        drifters[fname]['err_v'] = pd.Series(errv)
        drifters[fname]['err_mag'] = pd.Series(err_mag)
        drifters[fname]['err_dir'] = pd.Series(err_dir)
        # why was mTimes being used here???
        idx = np.asarray([np.argmin(np.abs(ncfile.Variables.matlabTime-x))\
                for x in oTimes])
        depth = [depth[k][i][-1] for k, i in enumerate(idx, start=0)]

        drifters[fname]['depth'] = pd.Series(depth)
        drifters[fname]['u_sim_speed'] = pd.Series(uspeedS)
        drifters[fname]['u_obs_speed'] = pd.Series(uspeedO)
        drifters[fname]['u_diff'] = pd.Series(udiffs)
        drifters[fname]['lon'] = pd.Series(olon)
        drifters[fname]['lat'] = pd.Series(olat)
        drifters[fname]['uobs'] = pd.Series(oU)
        drifters[fname]['vobs'] = pd.Series(oV)
        drifters[fname]['usim'] = pd.Series(mU)
        drifters[fname]['vsim'] = pd.Series(mV)

        # drifters[fname]['obs_timeseries'] = {'u' : oU, 'v' : oV}
        # drifters[fname]['mod_timeseries'] = {'u' : mU, 'v' : mV}
        # drifters[fname]['datetimes'] = datetimes
        drifters[fname]['mat_times'] = pd.Series(mTimes)
        drifters[fname]['obs_times'] = pd.Series(oTimes)

        # find wind speedi if in name
        if fname.split('.')[0].split('_')[-1].startswith('S'):
            wind = fname.split('.')[0].split('_')[-1]
        elif fname.split('.')[0].split('_')[-1].startswith('N'):
            wind = fname.split('.')[0].split('_')[-1]
        elif fname.split('.')[0].split('_')[-1].startswith('W'):
            wind = fname.split('.')[0].split('_')[-1]
        elif fname.split('.')[0].split('_')[-1].startswith('E'):
            wind = fname.split('.')[0].split('_')[-1]
        else:
            if debug:
                print '\tno recorded wind speed...'
            wind = '0'
        drifters[fname]['wind_speed'] = wind

        all_sdev.append(sdev_diff)
        all_mean.append(mean_diff)
        all_diff.extend(diffs)
        all_udiff.extend(udiffs)
        all_erru.extend(erru)
        all_errv.extend(errv)
        all_err_mag.extend(err_mag)
        all_err_dir.extend(err_dir)
        depths.extend(depth)
        mean_depth.append(np.mean(depth))

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

    if tide_opt is not None:
        print 'all drifters for ' + tide_opt + ' selected...'

    return drifters, all_mean, all_sdev, obs_speed, mod_speed, obs_uspeed, \
         mod_uspeed, all_diff, o_lat, o_lon, lon0, lat0, all_udiff, depths, \
         all_erru, all_errv, all_err_mag, all_err_dir, mean_depth


