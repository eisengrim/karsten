#! /usr/bin/python2.7


from __future__ import division
import numpy as np
import scipy as sp
from scipy.io import netcdf
import netCDF4 as nc
import sys, os
import os.path as osp




if __name__ == '__main__':


    filename = sys.argv[1]

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
    # creates matlab time
    mtime = self.jtime[:] + 678942.0
    print 'full temporal domain loaded...'
    print 'loading velocities...'

    u = data.variables['u'].data
    v = data.variables['v'].data
    w = data.variables['w'].data




