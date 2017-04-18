"""
Generates a CSV file containing x,y and lon,lat for the generation of particles in Thomas Lake's Swansea code. Lon/Lat is given in WGS84. Outputs text file file suitable for a particle input.

**Unfinished

Parameters:
    generateFood.py -f fvcom -c x_c y_c -d min_depth -o out -l x_r y_r
"""

import argparse
import numpy as np
from pyseidon import *
import sys

if  __name__ == '__main__':

    parser = argparse.ArgumentParser(description='usage: ' \
          + '')
    parser.add_argument('-f', action='store', default=None, dest='path2fvcom', \
            help='path to fvcom nc file', type=str, required=True)
    parser.add_argument('-o', action='store', default=None, dest='outpath', \
            help='path to outfile basename', type=str, required=True)
    parser.add_argument('-r', action='store', nargs=3, dest='range', \
            help='range for x, y, z rectangle', type=float)
    parser.add_argument('-c', action='store', nargs=3, dest='centre', \
            help='centre coordinates for bounding box in x,y,z')
    parser.add_argument('-n', action='store', type=float, dest='np', \
            help='number of particles', default=2000)
    parser.add_argument('-v', action="store", nargs=3, default=[1.0,1.0,0.75], \
            dest='vel', help='velocity range')
    parser.add_argument('--orient', action='store', nargs=3, \
            default=[np.pi, np.pi, np.pi], dest='orient', help='range of orientation')
    parser.add_argument('-p', action='store', required=True, dest='pfile', \
            help='input particle file')
    args = parser.parse_args()

    # load netcdf file
    if not args.path2fvcom.endswith('.nc'):
        sys.exit('not a valid netcdf4 file.')
    model = FVCOM(args.path2fvcom)

    idx = np.where(model.Grid.h > 0)

    print 'saving particle information...'

    with open(args.outpath, 'wb') as f:
        f.write(b"\nposition = ({:+.10d}, {+:.10d}, {:+.10d})\n\n".format(args.centre[0], \
                args.centre[1], args.centre[2]))

        fd = ""
        posx = np.random.uniform(args.centre[0]-args.range[0], \
                                 args.centre[0]+args.range[0], args.np)
        posy = np.random.uniform(args.centre[1]-args.range[1], \
                                 args.centre[1]+args.range[1], args.np)
        posz = np.random.uniform(args.centre[2]-args.range[2], \
                                 args.centre[2]+args.range[2], args.np)

        velx = np.random.uniform(-args.vel[0], \
                                 +args.vel[0], args.np)
        vely = np.random.uniform(-args.vel[1], \
                                 +args.vel[1], args.np)
        velz = np.random.uniform(-args.vel[2], \
                                 +args.vel[2], args.np)

        orix = np.random.uniform(-args.orient[0], \
                                 +args.orient[0], args.np)
        oriy = np.random.uniform(-args.orient[1], \
                                 +args.orient[1], args.np)
        oriz = np.random.uniform(-args.orient[2], \
                                 +args.orient[2], args.np)

        for i in xrange(args.np):







