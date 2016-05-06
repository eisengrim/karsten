from pyseidon import *
from interpolation_utils import *
import numpy as np
import scipy as sp
import sys, os
import os.path as osp
from datetime import datetime, timedelta

PATH2SIM = "/EcoII/acadiau_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/"
LOC = ['GP', 'DG', 'PP']
#BF = ['0.015', '0.012', '0.009']
BF = ['0.015']
SAVEDIR = "/array/home/119865c/karsten/vel_interp/"


if __name__ == 'main':

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
            ne_idx_gp = find_nearest_elems(trinodes, nnode)
        elif 'DG' not in done and 'DG' in dir_name:
            done.append('DG')
            print '\tfinding nearest element indices for DG...'
            ne_idx_dg = find_nearest_elems(trinodes, nnode)
        elif 'PP' not in done and 'PP' in dir_name:
            done.append('PP')
            print '\tfinding nearest element indices for PP...'
            ne_idx_pp = find_nearest_elems(trinodes, nnode)

        print '\tswapping axes...'
        # numpy is row major ordered - quicker to have elements as first dim
        u=np.swapaxes(ncfile.Variables.u,0,1)
        v=np.swapaxes(ncfile.Variables.v,0,1)
        w=np.swapaxes(ncfile.Variables.w,0,1)

        print '\tinitializing empty array...'
        varE = np.empty(len(nnode), dtype=object)

        if 'DG' in dir_name:
            for node in xrange(len(nnode)):
                varE[i] = np.mean(u[ne_idx_dg,:,:], axis=0)

        elif 'GP' in dir_name:
            for node in xrange(len(nnode)):
                varE[i] = np.mean(u[ne_idx_gp,:,:], axis=0)

        elif 'PP' in dir_name:
            for node in xrange(len(nnode)):
                varE[i] = np.mean(u[ne_idx_pp,:,:], axis=0)

        x = np.squeeze([np.vstack(varE[i]) for i in xrange(len(varEE))])

        print '\tcreating new subdirectory...'

        sim_name = dir_name[15:29]
        savepath = SAVEDIR + 'bfric_' + bf + '/' + loc + '_' + sim_name
        os.makedirs(savepath)
        savepath = savepath + loc + '_' + sim_name + '.mat'

        # benefit of saving in matlab: you don't need to know original shape
        # try in pickle, numpy.save or netcdf4??
        sp.io.savemat(savepath, mdict={'out': x}, oned_as='row')

        # to load the data:
        #     > matdata = scipy.io.loadmat(savepath)
        # to check if data is unchanged
        #     > assert np.all(x == matdata['out'])

        # var = interp_to_nodes_avg(ncfile.Variables.u, ncfile.Grid.nnode, \
        #        ncfile.Grid.trinodes)




