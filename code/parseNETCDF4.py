#! /usr/bin/python2.7

"""
Parses netCDF4 data from FVCOM in a similar manner to tawe-telemac-utils.

usage: python netCDF4.py loc bfric sim path/to/out/dir/
"""

from __future__ import division
import numpy as np
import scipy as sp
from scipy.io import netcdf
import netCDF4 as nc
import sys, os
import os.path as osp

PATH2SIM='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/'
PATH2VEL='/home/array/119865c/karsten/swansea/vel_interp/bfric_'


if __name__ == '__main__':

    loc = sys.argv[1]
    sim = sys.argv[3]
    bfric = sys.argv[2]
    outpath = sys.argv[4]

    if loc not in ['GP', 'PP', 'DG']:
        sys.exit('not a valid location tag.')
    if bfric not in ['0.015', '0.012', '0.009']:
        sys.exit('not a valid bottom friction.')

    mesh = loc + '_' + sim
    filename = PATH2SIM+'BFRIC_'+bfric+'/'+loc+'/'+sim+ \
                '/output/subdomain_'+loc+'1_0001.nc'
    dir_vel = PATH2VEL+bfric+'/'+mesh+'/'

    # loading netcdf file
    if not filename.endswith('.nc'):
        sys.exit("not a valid netcdf4 file.")

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
    nt = len(times)
    times = np.vstack((np.arange(nt), times)).T

    print 'loading spatial domain and connectivity...'
    # deal with nodal coordinates
    x = data.variables['x'].data
    y = data.variables['y'].data
    nnode = len(x)

    x = np.hstack((x,x,x,x,x,x,x,x,x,x)).astype(object)
    y = np.hstack((y,y,y,y,y,y,y,y,y,y)).astype(object)
    npoin = nnode*10

    # deal with elemental coordinates and connectivity
    nelem = len(data.variables['xc'].data)*9
    ndp = 6

    try:
        nv = np.transpose(data.variables['nv'].data) - 1
    except AttributeError:
        # exception due to nc.Dataset type data
        nv = np.transpose(data.variables['nv'][:]) - 1
    # since neighbouring nodes are preserved throughout layers,
    # triangular prisms may be constructed by copying nv...
    nv = np.hstack((nv,nv+nnode))
    nv = np.vstack((nv, nv+nnode, nv+2*nnode, nv+3*nnode, nv+4*nnode, \
                    nv+5*nnode, nv+6*nnode, nv+7*nnode, nv+8*nnode)).reshape(-1)

    # create element and node lists
    elems = np.arange(nelem)
    elems = np.vstack((elems,elems,elems,elems,elems,elems)).reshape(-1,order='F')
    nodes = np.arange(npoin).astype(object)

    print 'loading velocity files...'
    # deal with velocities
    u = np.load(dir_vels+'_u.pkl', mmap_mode = 'r')
    v = np.load(dir_vels+'_v.pkl', mmap_mode = 'r')
    w = np.load(dir_vels+'_w.pkl', mmap_mode = 'r')

    # flatten velocity arrays in all but time dim
    u = u.reshape(u.shape[0], -1).astype(object)
    v = v.reshape(v.shape[0], -1).astype(object)
    w = w.reshape(w.shape[0], -1).astype(object)

    nvars = 4
    var0 = 'Z'
    var1 = 'U'
    var2 = 'V'
    var3 = 'W'

    # deal with height coordinate
    print 'loading z coordinate variable...'
    h = data.variables['h'].data
    siglay = data.variables['siglay'].data
    zeta = data.variables['zeta'].data

    # reshape data to be broadcast together
    depth = h + zeta
    depth = np.hstack((depth,depth,depth,depth,depth, \
                            depth,depth,depth,depth,depth))

    z = np.multiply(depth,siglay.reshape(-1))

    # save to file
    print 'saving to files...'
    if not osp.exists(outpath):
        os.makedirs(outpath)

    if outpath[-1] != '/':
        outpath.append('/')

    outfile = outpath+mesh

    np.savetxt(outfile + '.x.txt', np.vstack((nodes, x)).T, header=str(npoin), \
               delimiter='\t')
    np.savetxt(outfile + '.y.txt', np.vstack((nodes, y)).T, header=str(npoin), \
               delimiter='\t')
    np.savetxt(outfile + '.conn.txt', np.vstack((elems, nv)).T, delimiter='\t', \
               header='{}\t{}'.format(nelem, ndp))
    np.savetxt(outfile + '.times.txt', times, delimiter='\t'. header=str(nt))

    with open(outfile + '.vars.txt', 'w') as f:
        out_str='{}\n{}\t{}\n{}\t{}\n{}\t{}\n{}\t{}'\
                .format(nvars, 0, var0, 1, var1, 2, var2, 3, var3)
        f.write(out_str)

    for t in np.xrange(nt):
        np.savetxt(outfile+.'var0.t'+str(t)+'.txt', np.vstack((nodes, z[t])).T, \
                   delimiter = '\t')
        np.savetxt(outfile+.'var1.t'+str(t)+'.txt', np.vstack((nodes, u[t])).T, \
                   delimiter = '\t')
        np.savetxt(outfile+.'var2.t'+str(t)+'.txt', np.vstack((nodes, v[t])).T, \
                   delimiter = '\t')
        np.savetxt(outfile+.'var3.t'+str(t)+'.txt', np.vstack((nodes, w[t])).T, \
                   delimiter = '\t')

    print 'creating .ini file...'
    with open(outfile + '.ini', 'w') as f:
        to_write = """[paths]
        basename = %s
        output = %s

        [settings]
        dimensions = 3
        particle_steps = 1

        [indexes]
        fish = 0
        noise = 0
        z = 0
        u = 1
        v = 2
        w = 3
        """

        f.write(to_write % (outfile, outfile + '_output/'))

    if not osp.exists(outfile + '_output/'):
        os.makedirs(outfile + '_output/')

    print 'all done!'
