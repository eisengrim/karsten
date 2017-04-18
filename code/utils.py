#! /usr/env/python2.7

from datetime import datetime, timedelta
import numpy as np
import scipy.io as sio
import netCDF4 as nc
import ntpath as ntp


def closest_vals(arr1, arr2):
    """
    Given two arrays, returns an array of length arr1 with
    the indices corresponding to the closest values in arr2.
    """
    arr1 = np.array(arr1).T
    arr1 = arr1[:, np.newaxis]
    arr2 = np.array(arr2).T

    return np.argmin(abs(arr2 - arr1), axis=1)


def closest_dist(x, y, x_list, y_list):
    """
    Given two pairs of arrays (x,y), finds the closest corresponding
    values in a different pair, (x_list,y_list) and returns their indices.
    """
    points = np.array([x, y]).T
    points_list = np.array([x_list, y_list]).T

    dpt0 = points_list[:, 0] - points[:, 0, np.newaxis]
    dpt1 = points_list[:, 1] - points[:, 1, np.newaxis]

    return np.argmin((dpt0*dpt0 + dpt1*dpt1), axis=1)


def path_leaf(path):
    head, tail = ntp.split(path)
    return tail or ntp.basename(head)


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


def drift_times(name, debug=False):
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
        drft = sio.loadmat(name, struct_as_record=False, squeeze_me=True)
        times = drft['gps_observation'][0][0][0][0]
    except KeyError:
        try:
            times = drft['time'][0][0][0][0]
        except IndexError:
            times = drft['time']
    # grab the times, sort them, and convert them to strings

    start, end = np.sort(times)[0], np.sort(times)[-1]
    return start, end


def fvcom_time(filename, debug=False):
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
    Convert modified Julian datetime to matlab datenum.
    '''
    y = x + 678942
    if debug:
        print 'converting from modified julian date to matlab datenum...'
    return y


def dn2dt(datenum, debug=False):
    """
    Convert matlab datenum to python datetime.

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


