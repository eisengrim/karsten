#! /usr/env/python2.7

from datetime import datetime, timedelta
import numpy as np
import scipy.io as sio
import netCDF4 as nc


def broadcast(arr1, arr2, debug=False):
    '''
    Finds the closest elements of two arrays for broadcasting.
    Returns a list of indices. arr1 is compared to arr2.

    arr1, arr2 :: arrays of different length
    len(arr1) > len(arr2)
    '''
    return [np.abs(arr1 - t).argmin() for t in arr2]


def haversine(lon, lat):
    '''
    Computes the distance between consecutive lat/lon coordinates in m.

    '''
    R = 6371000
    dLon = np.deg2rad(np.diff(lon))
    dLat = np.deg2rad(np.diff(lat))

    a = np.square(np.sin(dLat/2)) + np.square(np.sin(dLon/2)) \
            + (np.cos(np.deg2rad(lat[:-1]))*np.cos(np.deg2rad(lat[1:])))
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))

    return R*c


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
            print 'loading file {}...'.format(name)
        drft = sio.loadmat(name)
        times = drft['gps_observation'][0][0][0][0]
    except KeyError:
        times = drft['time'][0][0][0][0]

    # grab the times and convert them to strings
    start, end = times[0], times[-1]
    return start, end


def fvcomTime(filename, debug=False):
    '''
    Identify the timespan for a given FVCOM file.
    '''
    # iterate through the fvcom files

    fvcom = nc.Dataset(filename)
    if debug:
        print 'nc fvcom file loaded....\n reading variables...'

    try:
        time = fvcom.variables['time'][:]
    except KeyError:
        time = fvcom.variables['time_JD'][:]

    # grab the times and convert them to strings
    start, end = mjd2num(time[0]), mjd2num(time[-1])
    return start, end


def mjd2num(x, debug=False):
    '''
    Convert FVCOM time to matlab datenum.
    '''
    y = x + 678942
    if debug:
        print 'converting from modified julian date to matlab datenum...'
    return y


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


