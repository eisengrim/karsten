#! /usr/env/python2.7

import numpy as np
from pyseidon import *
import scipy.special as sps

# local imports
from drifterUtils import dn2dt

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


def calculateBias(ncfile, files, loc, tide_opt=None, debug=False):
    """
    Compiles necessary data, calculates individual biases and
    takes the mean and standard deviation. Also calculates the
    bias for ALL drifter files. Returns the mean biases, the drifter
    numbers, the model and observed speeds the the std deviations.

    input:
        - ncfile : FVCOM object
        - files : list of matlab filenames in directory
        - loc : location tag
        - tide_opt : 'ebb' or 'flood' ... returns only drifter data for tide
        - sim_name : name of simulation directory used
    """

    if debug:
        print '{} drifters will be analysed...'.format(len(files))
    #    print 'calculating depth...'
    # ncfile.Util3D.depth()

    data_count = 0
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
        drifters[fname[48:-4]] = {}
        if debug:
            print 'creating drifter object {}...'.format(i)
            print fname

        drift = Drifter(fname, debug=False)

        # tests tide condition
        if tide_opt is not None:
            tide = str(drift.Data['water_level'].tide)
            if tide != tide_opt:
                if debug:
                    print fname + ' not in tide selection'
                continue

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

        uspeedS = np.sqrt(mU**2 + mV**2)
        uspeedO = np.sqrt(oU**2 + oV**2)

        # index finder is only in dvt branch
        depth = [ncfile.Util3D.depth_at_point(j,k) for j,k in zip(olon,olat)]
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
        mean_bias = np.mean(diffs)
        sdev_bias = np.std(diffs)
        erru = np.subtract(mU, oU)
        errv = np.subtract(mV, oV)
        err_mag = np.sqrt(erru**2 + errv**2)
        err_dir = np.arctan(np.divide(errv, erru))/np.pi * 180

        # info on # of data this drifter has and other individual data
        drifters[fname[48:-4]]['filename'] = fname[48:-4]
        fname = fname[48:-4]
        drifters[fname]['num'] = i
        drifters[fname]['(start, stop)'] = (data_count,data_count+len(olon)-1)
        data_count += len(olon)
        drifters[fname]['tide'] = str(drift.Data['water_level'].tide)
        drifters[fname]['beta'] = drift.Data['water_level'].beta
        ind = [np.argmin(np.abs(drift.Data['velocity'].vel_time - x)) \
                for x in valid.Variables.struct['mod_time']]
        drifters[fname]['alpha'] = drift.Data['velocity'].alpha[ind]
        drifters[fname]['obs_speed'] = speedO
        drifters[fname]['sim_speed'] = speedS
        drifters[fname]['bias'] = diffs
        drifters[fname]['mean_bias'] = mean_bias
        drifters[fname]['sdev_bias'] = sdev_bias
        drifters[fname]['err_u'] = erru
        drifters[fname]['err_v'] = errv
        drifters[fname]['err_mag'] = err_mag
        drifters[fname]['err_dir'] = err_dir
        idx = [np.argmin(np.abs(ncfile.Variables.matlabTime-x)) for x in mTimes]
        depth = [depth[k][i][-1] for k, i in enumerate(idx)]
        drifters[fname]['depth'] = depth
        drifters[fname]['u_sim_speed'] = uspeedS
        drifters[fname]['u_obs_speed'] = uspeedO
        drifters[fname]['u_bias'] = udiffs
        drifters[fname]['lon'] = olon
        drifters[fname]['lat'] = olat
        drifters[fname]['obs_timeseries'] = {'u' : oU, 'v' : oV}
        drifters[fname]['mod_timeseries'] = {'u' : mU, 'v' : mV}
        drifters[fname]['datetimes'] = datetimes
        drifters[fname]['mat_times'] = mTimes

        # find wind speed
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
                print '\tno wind speed...'
            wind = '0'
        drifters[fname]['wind_speed'] = wind

        all_sdev.append(sdev_bias)
        all_mean.append(mean_bias)
        all_bias.extend(diffs)
        all_ubias.extend(udiffs)
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
         mod_uspeed, all_bias, o_lat, o_lon, lon0, lat0, all_ubias, depths, \
         all_erru, all_errv, all_err_mag, all_err_dir, mean_depth


