#! /usr/bin/python2.7

"""
Parses netCDF4 data from FVCOM in a similar manner to tawe-telemac-utils.

usage: python netCDF4.py path/to/file.nc path/to/vel/simdate path/to/output
"""

from __future__ import division
import numpy as np
import scipy as sp
from scipy.io import netcdf
import netCDF4 as nc
import sys, os
import os.path as osp


if __name__ == '__main__':

    filename = sys.argv[1]
    dir_vels = sys.argv[2]
    outpath = sys.argv[3]

    # loading netcdf file
    if filename.endswith('.nc'):
        print 'retrieving data from ' + filename + '...'

    try:
        data = netcdf.netcdf_file(filename, 'r', mmap=True)
    except (OverflowError, TypeError, ValueError) as e:
        data = nc.Dataset(filename, 'r', format='NETCDF4_CLASSIC')

    print 'loading timestamps...'
    try:
        jtime = data.variables['time'].data
    except AttributeError as e:
        # exception due to nc.Dataset type data
        if e == AttributeError: jtime = data.variables['time'][:]
    print 'full temporal domain loaded...'

    # deal with time data
    # converts julian time to seconds elapsed
    times = (jtime - jtime[0])*8640/0.1
    nt - len(times)

    print 'loading spatial domain and connectivity...'
    # deal with nodal coordinates
    x = data.variables['x'].data
    y = data.variables['y'].data
    npoin = len(x)

    # deal with elemental coordinates and connectivity
    nelem = len(data.variables['xc'].data)
    ndp = 6

    try:
        nv = np.transpose(data.variables['nv'].data) - 1
    except AttributeError:
        # exception due to nc.Dataset type data
        nv = np.transpose(data.variables['nv'][:]) - 1
    # since neighbouring nodes are preserved throughout layers,
    # triangular prisms may be constructed by copying nv...
    nv = np.hstack((nv,nv))

    print 'loading velocity files...'
    # deal with velocities
    u = np.load(dir_vels+'_u.pkl', mmap_mode = 'r')
    v = np.load(dir_vels+'_v.pkl', mmap_mode = 'r')
    w = np.load(dir_vels+'_w.pkl', mmap_mode = 'r')

    nvars = 4
    var0 = 'Z'
    var1 = 'U'
    var2 = 'V'
    var3 = 'W'

    # question: how to get z coordinate? how to print out velocities?
    print 'loading z coordinate variable...'
    h = data.variables['h'].data
    zeta = data.variables['zeta'].data
    z = h + zeta

    # save to file
    print 'saving to file...'


