#! /usr/bin/python2.7

from pyseidon import *
from interpolation_utils import *
import numpy as np
import scipy as sp
from scipy.io import netcdf
import sys, os
import os.path as osp

PATH2SIM = "/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
LOC = ['GP']#,'DG','PP']
#BF = ['0.015', '0.012', '0.009']
BF = ['0.015']
SAVEDIR = "/array/home/119865c/karsten/vel_interp/"
ELEMS = "/array/home/119865c/karsten/nearest_elems/"

if __name__ == '__main__':

    print 'locating all directories...'
    for bf in BF:
        for loc in LOC:
            print 'locating files for bf={} and loc={}...'.format(bf, loc)
            path2sim = PATH2SIM + 'BFRIC_' + bf + '/' + loc + '/'

            dirs = os.listdir(path2sim)

            sim_path = [path2sim + file + '/output/subdomain_' + loc \
                        + '1_0001.nc' for file in dirs]

            for path in sim_path:
                if not osp.exists(path) or not osp.isfile(path):
                    print '\tfvcom file {} not found. removing...'.format(path)
                    sim_path.remove(path)

            if not osp.exists(SAVEDIR):
                sys.exit('save directory not found. exiting...')

    if len(sim_path) == 0:
        sys.exit('no ncfiles found')

    done = []
    for dir_name in sim_path:
        print 'loading fvcom object...'
        ncfile = FVCOM(dir_name, debug=False)
        print 'ncfile for {} loaded.'.format(dir_name)
        print '\tshape of speed is: {}.'.format(ncfile.Variables.u.shape)
        print '\tshape of trinodes is: {}.'.format(ncfile.Grid.trinodes.shape)
        print '\tnnode: {}, nele: {}'.format(ncfile.Grid.nnode, ncfile.Grid.nele)

        trinodes = ncfile.Grid.trinodes
        nnode = ncfile.Grid.nnode

        if 'GP' not in done and 'GP' in dir_name:
            done.append('GP')
            print '\tfinding nearest element indices for GP...'
            fname = ELEMS + 'GP_ne_idx.dat'
            ne_idx_gp = find_nearest_elems(trinodes, nnode, masked = False)
            if not osp.exists(fname):
                to_dump = np.ma.masked_invalid(np.asarray([np.append(\
                        ne_idx_gp[i], [np.nan]*(7-len(ne_idx_gp[i]))) \
                        for i in xrange(len(ne_idx_gp))]))
                np.ma.dump(to_dump, fname)

        elif 'DG' not in done and 'DG' in dir_name:
            done.append('DG')
            print '\tfinding nearest element indices for DG...'
            fname = ELEMS + 'DG_ne_idx.dat'
            ne_idx_dg = find_nearest_elems(trinodes, nnode, masked = False)
            if not osp.exists(fname):
                to_dump = np.ma.masked_invalid(np.asarray([np.append( \
                        ne_idx_dg[i], [np.nan]*(7-len(ne_idx_dg[i]))) \
                        for i in xrange(len(ne_idx_dg[i]))]))
                np.ma.dump(to_dump, fname)

        elif 'PP' not in done and 'PP' in dir_name:
            done.append('PP')
            print '\tfinding nearest element indices for PP...'
            fname = ELEMS + 'PP_ne_idx.dat'
            ne_idx_pp = find_nearest_elems(trinodes, nnode, masked = False)
            if not osp.exists(fname):
                to_dump = np.ma.masked_invalid(np.asarray([np.append( \
                        ne_idx_pp[i], [np.nan]*(7-len(ne_idx_pp[i]))) \
                        for i in xrange(len(ne_idx_pp[i]))]))
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

        print '\tcreating new subdirectory...'
        sim_name = dir_name[-43:-29]
        print '\t\tfiles will be placed in {}...'.format(sim_name)
        savepath = SAVEDIR + 'bfric_' + bf + '/' + loc + '_' + sim_name
        os.makedirs(savepath)
        savepath = savepath + loc + '_' + sim_name

        if 'DG' in dir_name:
            ne_idx = ne_idx_dg
        elif 'GP' in dir_name:
            ne_idx = ne_idx_gp
        elif 'PP' in dir_name:
            ne_idx = ne_idx_pp

        print '\taveraging for u...'
        for node in xrange(nnode):
            uvarE[node] = np.mean(ncfile.Variables.u[ne_idx[node],:,:], axis=0)
        uvarE = np.squeeze([np.vstack(uvarE[i]) for i in xrange(nnode)])
        print '\t\trestoring axes...'
        uvarE = np.swapaxes(uvarE, 0, 2)
        print '\t\tsaving as pickle...'
        uvarE.dump(savepath + '_u.pkl')
        print '\tsize of mean speed(s): {}'.format(uvarE.shape)
        del uvarE

        print '\taveraging for v...'
        for node in xrange(nnode):
            vvarE[node] = np.mean(ncfile.Variables.v[ne_idx[node],:,:], axis=0)
        vvarE = np.squeeze([np.vstack(vvarE[i]) for i in xrange(nnode)])
        print '\t\trestoring axes...'
        vvarE = np.swapaxes(vvarE, 0, 2)
        print '\t\tsaving as pickle...'
        vvarE.dump(savepath + '_v.pkl')
        print '\t\tsize of mean speed(s): {}'.format(vvarE.shape)
        del vvarE

        for node in xrange(nnode):
            wvarE[node] = np.mean(ncfile.Variables.w[ne_idx[node],:,:], axis=0)
        wvarE = np.squeeze([np.vstack(wvarE[i]) for i in xrange(nnode)])
        print '\t\trestoring axes...'
        wvarE = np.swapaxes(wvarE, 0, 2)
        print '\t\tsaving as pickle...'
        wvarE.dump(savepath + '_w.pkl')
        print '\t\tsize of mean speed(s): {}'.format(wvarE.shape)
        del wvarE

        # clear structure to save some space
        del ncfile

        # benefit of saving in matlab: you don't need to know original shape
        # try in pickle, numpy.save or netcdf4??
        # sp.io.savemat(savepath + '_u.mat', mdict={'out': uvarE}, oned_as='row')
        # sp.io.savemat(savepath + '_v.mat', mdict={'out': vvarE}, oned_as='row')
        # sp.io.savemat(savepath + '_w.mat', mdict={'out': wvarE}, oned_as='row')

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




