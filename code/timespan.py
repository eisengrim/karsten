import scipy.io as sio
import sys
import numpy as np
from datetime import datetime, timedelta
import h5py
import netCDF4 as nc


# USAGE:
# python timespan.py TYPE FILE1 FILE2 FILE3 ...
# TYPE can be one of the following: FVCOM, ADCP, DRIFTER
# all files have to be the same type


def mjd2num(x):
    '''
    Convert FVCOM time to matlab datenum.
    '''
    y = x + 678942

    return y


def dn2dt(datenum):
    '''
    Convert matlab datenum to python datetime.
    '''
    return datetime.fromordinal(int(datenum)) + \
        timedelta(days=datenum % 1) - \
        timedelta(days=366)


def fvcomTimes(filenames):
    '''Identify and print the list of timespans for each FVCOM file
    '''
    # iterate through the fvcom files
    for name in filenames:
        fvcom = nc.Dataset(name)

        try:
            time = fvcom.variables['time'][:]
        except KeyError:
            time = fvcom.variables['time_JD'][:]

        # grab the times and convert them to strings
        start, end = mjd2num(time[0]), mjd2num(time[-1])
        start, end = dn2dt(start), dn2dt(end)
        start_str, end_str = str(start), str(end)
        start_str = start_str.split('.')[0]
        end_str = end_str.split('.')[0]
        name = name.split('/')[-1]
        print "{} runs from {} to {}".format(name, start_str, end_str)


def adcpTimes(filenames):
    '''Identify and print the list of timespans for each ADCP
    '''
    # iterate through adcp files
    for name in filenames:
        # don't open raw files or station files
        if 'raw' in name.lower() or 'station' in name.lower():
            continue

        try:
            adcp = sio.loadmat(name)
            times = adcp['time'][0][0][0][0]
        except NotImplementedError:
            adcp = h5py.File(name)
            times = np.rot90(adcp['time']['mtime'][:])[0]

        # grab the times and convert them to strings
        start, end = dn2dt(times[0]), dn2dt(times[-1])
        start_str, end_str = str(start), str(end)
        start_str = start_str.split('.')[0]
        end_str = end_str.split('.')[0]
        name = name.split('/')[-1]
        print "{} runs from {} to {}".format(name, start_str, end_str)


def driftTimes(filenames):
    '''Identify and print the list of timespans for each ADCP
    '''
    # iterate through the drifter files
    for name in filenames:
        try:
            drft = sio.loadmat(name)
            times = drft['gps_observation'][0][0][0][0]
        except KeyError:
            times = drft['time'][0][0][0][0]

        # grab the times and convert them to strings
        start, end = dn2dt(times[0]), dn2dt(times[-1])
        start_str, end_str = str(start).split('.')[0], str(end).split('.')[0]
        name = name.split('/')[-1]
        print "{} runs from {} to {}".format(name, start_str, end_str)


if __name__ == '__main__':

    argc = len(sys.argv)

    if argc < 3:
        print "Insufficient command line arguments! Must add a type and" + \
            "at least one file."
        sys.exit(1)

    data_type = sys.argv[1]

    # find the data timespan based on the type of data
    if data_type == 'ADCP':
        adcpTimes(sys.argv[2:])
    elif data_type == 'FVCOM':
        fvcomTimes(sys.argv[2:])
    elif data_type == 'DRIFTER':
        driftTimes(sys.argv[2:])
