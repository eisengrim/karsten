#! /usr/bin/python2.7
"""
Linearly interpolates FVCOM grid velocities to nodes.

Usage:
    python runInterpAvg.py BF LOC PATH2SIM.nc SAVEDIR/ ELEMS/ MESH

Example Inputs:
PATH2SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
LOC = 'GP'
BF = '0.015'
SAVEDIR = "/EcoII/acadia_uni/projects/drifters/swansea/vel_interp/"
ELEMS = "/array/home/119865c/karsten/nearest_elems/"
MESH = "2013_Aug_08_3D" <- whatever your mesh is named, usually date of run
"""


from pyseidon import *
from interpolation_utils import *
import numpy as np
import scipy as sp
from scipy.io import netcdf
import sys, os
import os.path as osp
import h5py

if __name__ == '__main__':
    if len(sys.argv) == 7:
        bf = str(sys.argv[1])
        loc = sys.argv[2]
        sim_path = sys.argv[3]
        savedir = sys.argv[4]
        elemdir = sys.argv[5]
        mesh_name = sys.argv[6]
    else:
        sys.exit('insufficient number of command line args.')

    print 'locating files for bf={} and loc={}...'.format(bf, loc)

    if not osp.exists(sim_path) or not osp.isfile(sim_path):
        print sys.exit('fvcom file {} not found.'.format(sim_path))

    dir_name = str(sim_path)
    print 'loading fvcom object...'
    ncfile = FVCOM(dir_name, debug=False)
    print 'ncfile for {} loaded.'.format(dir_name)
    print '\tshape of speed is: {}.'.format(ncfile.Variables.u.shape)
    print '\tshape of trinodes is: {}.'.format(ncfile.Grid.trinodes.shape)
    print '\tnnode: {}, nele: {}'.format(ncfile.Grid.nnode, \
                ncfile.Grid.nele)

    trinodes = ncfile.Grid.trinodes
    nnode = ncfile.Grid.nnode

    print '\tfinding nearest element indices for {}...'.format(loc)
    fname = elemdir + '{}_ne_idx.dat'.format(loc)
    ne_idx = find_nearest_elems(trinodes, nnode, masked = False)
    if not osp.exists(fname):
        to_dump = np.ma.masked_invalid(np.asarray([np.append(\
                        ne_idx[i], [np.nan]*(7-len(ne_idx[i]))) \
                        for i in xrange(len(ne_idx))]))
        np.ma.dump(to_dump, fname)

    print '\tswapping axes...'
    # numpy is row major ordered - quicker to have elements as first dim
    ncfile.Variables.u = np.swapaxes(ncfile.Variables.u,0,2)
    ncfile.Variables.v = np.swapaxes(ncfile.Variables.v,0,2)
    ncfile.Variables.w = np.swapaxes(ncfile.Variables.w,0,2)

    print '\tinitializing empty arrays...'
    uvarE = np.empty(nnode, dtype=object)
    vvarE = np.empty(nnode, dtype=object)
    wvarE = np.empty(nnode, dtype=object)
    print 'shape of arrays: ', uvarE.shape

    savepath = savedir + 'bfric_' + bf + '/' + loc + '_' + mesh_name + '/'

    if osp.exists(savepath):
        sys.exit('{} savedir already created.'.format(mesh_name))

    print '\tcreating new subdirectory...'
    print '\t\tfiles will be placed in {}...'.format(mesh_name)

    os.makedirs(savepath)
    savepath = savepath + loc + '_' + mesh_name

    print '\taveraging for u...'
    for node in xrange(nnode):
        uvarE[node] = np.mean(ncfile.Variables.u[ne_idx[node],:,:], axis=0)
    uvarE = np.squeeze([np.vstack(uvarE[i]) for i in xrange(nnode)])
    print '\t\trestoring axes...'
    uvarE = np.swapaxes(uvarE, 0, 2)
    print '\t\tsaving as pickle...'
    h5f = h5py.File(savepath+'_u.h5', 'w')
    h5f.create_dataset('u_interp', data=uvarE)
    h5f.close()
    # uvarE.dump(savepath + '_u.pkl')
    print '\tsize of mean speed(s): {}'.format(uvarE.shape)
    del uvarE

    print '\taveraging for v...'
    for node in xrange(nnode):
        vvarE[node] = np.mean(ncfile.Variables.v[ne_idx[node],:,:], axis=0)
    vvarE = np.squeeze([np.vstack(vvarE[i]) for i in xrange(nnode)])
    print '\t\trestoring axes...'
    vvarE = np.swapaxes(vvarE, 0, 2)
    print '\t\tsaving as pickle...'
    h5f = h5py.File(savepath+'_v.h5', 'w')
    h5f.create_dataset('v_interp', data=vvarE)
    h5f.close()
    # vvarE.dump(savepath + '_v.pkl')
    print '\t\tsize of mean speed(s): {}'.format(vvarE.shape)
    del vvarE

    for node in xrange(nnode):
        wvarE[node] = np.mean(ncfile.Variables.w[ne_idx[node],:,:], axis=0)
    wvarE = np.squeeze([np.vstack(wvarE[i]) for i in xrange(nnode)])
    print '\t\trestoring axes...'
    wvarE = np.swapaxes(wvarE, 0, 2)
    h5f = h5py.File(savepath+'_w.h5', 'w')
    h5f.create_dataset('w_interp', data=wvarE)
    h5f.close()
    print '\t\tsaving as pickle...'
    # wvarE.dump(savepath + '_w.pkl')
    print '\t\tsize of mean speed(s): {}'.format(wvarE.shape)
    del wvarE

    # clear structure to save some space
    del ncfile

    # benefit of saving in matlab: you don't need to know original shape
    # try in pickle, numpy.save or netcdf4??
    # sp.io.savemat(savepath + '_u.mat', mdict={'out': uvarE}, \
    #               oned_as='row')
    # sp.io.savemat(savepath + '_v.mat', mdict={'out': vvarE}, \
    #               oned_as='row')
    # sp.io.savemat(savepath + '_w.mat', mdict={'out': wvarE}, \
    #               oned_as='row')

    # save as netcdf
    # u = netcdf.netcdf_file(savepath + '_u.nc', 'w')
    # u.history = 'BF/LOC/FILE//VEL: ' + bf + loc + sim_name + 'u'
    # u.createDimension('u', len(

    # v = netcdf.netcdf_file(savepath + '_v.nc', 'w')
    # v.history = 'BF/LOC/FILE//VEL: ' + bf + loc + sim_name + 'v'
    # w = netcdf.netcdf_file(savepath + '_w.nc', 'w')
    # w.history = 'BF/LOC/FILE//VEL: ' + bf + loc + sim_name + 'w'


    # to load the data:
    #     > matdata = scipy.io.loadmat(savepath)
    # to check if data is unchanged
    #     > assert np.all(x == matdata['out'])
    # or as a pickle, load with numpy.load or pickle.load...


    # var = interp_to_nodes_avg(ncfile.Variables.u, ncfile.Grid.nnode, \
    #        ncfile.Grid.trinodes)
