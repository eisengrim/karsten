"""
Generates a CSV file containing x,y and lon,lat for the generation of attracting sources in Thomas Lake's Swansea code. Lon/Lat is given in WGS84. Outputs text file file suitable for a field var input.

Parameters:
    generateFood.py -f fvcom -b min_lon max_lon min_lat max_lat -d depth -o out -r radius -s str
"""

import argparse
import numpy as np
from pyseidon import *
import sys

if  __name__ == '__main__':

    parser = argparse.ArgumentParser(description='usage: ' \
          + 'generateFood.py -f fvcom -b min_lon,max_lon,min_lat,max_lat ' \
          + '-d depth -o out -r radius -s str')
    parser.add_argument('-f', action='store', default=None, dest='path2fvcom', \
            help='path to fvcom nc file', type=str)
    parser.add_argument('-b', action='store', default=None, dest='box', \
            help='bounding box of for a region within the grid', \
            metavar='minlon maxlon minlat maxlat',type=float, nargs=4)
    parser.add_argument('-o', action='store', default=None, dest='outpath', \
            help='path to outfile basename', type=str)
    parser.add_argument('-r', action='store', default=None, dest='radius', \
            help='radius of constant strength', type=float)
    parser.add_argument('-s', action='store', default=None, dest='strength', \
            help='strength of food source', type=float)
    parser.add_argument('-d', action='store', default=None, dest='depth', \
            help='maximum depth tolerated', type=float)

    args = parser.parse_args()

    # load netcdf file
    if not args.path2fvcom.endswith('.nc'):
        sys.exit('not a valid netcdf4 file.')

    model = FVCOM(args.path2fvcom, ax=args.box)
    print 'h is of size ' + str(model.Grid.h.shape[0]) + '...'

    idx1 = np.where(model.Grid.h > 0)
    idx2 = np.where(model.Grid.h < args.depth)
    idx = np.intersect1d(idx1, idx2)

    lon = model.Grid.lon[idx]
    lat = model.Grid.lat[idx]
    rads = np.tile(args.radius, len(lon))
    stren = np.tile(args.strength, len(lon))

    print 'lon, lat of size ' + str(lon.shape[0]) + ' ' + str(lat.shape[0]) + '...'
    print 'saving lon/lat info...'

    with open(args.outpath + '_ll.txt', 'wb') as f:
        f.write(b'{}\t2\n'.format(len(lon)))
        np.savetxt(f, np.vstack((lat, lon, np.zeros(len(lon)), rads, stren)).T, \
                fmt='%f,%f,%f\t%f\t%f')

    # np.savetxt(args.outpath + '_ll.txt', np.vstack((lat, lon)).T, fmt='%f, %f')
    # dat = np.vstack((lat,lon)).T

    # nd = []
    # from ll2utm import ll2utm
    # for row in dat:
    #     a = ll2utm(row[0], row[1])
    #     a = list(a)
    #     a.extend(np.delete(row, [0,1], 0))
    #     nd.append(a)
