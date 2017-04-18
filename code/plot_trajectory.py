#! /usr/env/python2.7

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from matplotlib.dates import DateFormatter
import numexpr as ne
import seaborn as sns
from datetime import datetime
import os.path as osp
import os
import sys

from utils import *
from color_map import *
from drifter_plots import *

def trajectoryPlots(ncfile, drift, loc, date, \
                fname=None, tight=False, debug=False, \
                plot=False, bounds=[]):
    """
    Compiles necessary data, creates plots and saves / shows them all. The plot is
    a spatially varying color map of speed with a drifter trajectory superposed and
    a comparitive speed-time plot for each drifter in the given directory.

    input:
        - ncfile : PySeidon FVCOM object
        - drift : PySeidon Drifter object
        - loc : location tag (str)
        - savepath : savepath to be used (if plot=False)
        - date : date of simulation used (str)
        - bfric : bottom friction value (str)
        - fname : drifter filename
        - plot : True if plot is shown and not passed
        - tight : boolean, True if subdomain region is to be constricted
        - ratio : ratio to adjust model data
        - bounds : array of min max lon lat *priority over tight*
    returns:
        - fig1, a trajectory, fig2, a timeseries
    """

    sns.set(font="serif")

    # find the location centre for flood/tide split calculation
    # not yet working...
    if loc == 'GP':
        centre = [-66.33906, 44.26898]
        if tight and not bounds:
            bounds = [-66.355, -66.31, 44.245, 44.2925]
        elif not bounds:
            bounds = []
    elif loc == 'DG':
        centre = [-65.76000, 44.67751]
        if tight and not bounds:
            bounds = [-65.775, -65.77, 44.665, 44.69]
        elif not bounds:
            bounds = []
    elif loc == 'PP':
        centre = [-66.206924, 44.389368]
        # find out the tightness required for PP
        if tight and not bounds:
            bounds = [-66.225, -66.195, 44.37, 44.41]
        elif not bounds:
            bounds = []
    elif 'MP' in loc:
        centre = [-64.418942, 45.3499696]
        if tight and not bounds:
            bounds = [-64.51, -64.31, 45.31, 45.38]
        elif not bounds:
            bounds = []

    # computes variable norm
    fI, eI, _, _ = ncfile.Util2D.ebb_flood_split_at_point(centre[0], centre[1])
    ncfile.Util3D.velo_norm()

    if debug:
        print 'creating plots...'

    if savepath:
        savepath = savepath + 'bfric_' + bfric + '/' + loc + '_' + date
        if ratio != 1.0:
            savepath = savepath + '/with_ratio_{}'.format(str(ratio))

        # creates a subdirectory, so as not to overwrite existing files
        # if debug:
        #     print 'creating new subdirectory...'
        # now = datetime.now()
        # now = now.strftime("%Y%m%d")
        if not osp.exists(savepath):
             os.makedirs(savepath)
        # else:
        #     savepath = savepath + '/_' + now
        #     os.makedirs(savepath)
        savepath = savepath + '/'

    # check = checkIDs(fname, debug=True)
    # if check is not None:
    #     print "warning: more than one drifter in this file!"

    # creates drifter object window for flow map
    if debug:
        print 'creating drifter object window...'
    tModel = ncfile.Variables.matlabTime
    tDrift = drift.Variables.matlabTime
    win1 = (np.abs(tModel-tDrift.min())).argmin()
    win2 = (np.abs(tModel-tDrift.max())).argmin()

    # tests tide condition
    if debug:
        print 'testing tide conditions...'
    try:
        tide = str(drift.Data['water_level'].tide)
    except:
        tide = np.empty(len(drift.Data['tide_time']), dtype="S6")
        tide[drift.Data['tide_time'] > 0.5] = "flood"
        tide[drift.Data['tide_time'] < 0.5] = "ebb"

    # averages velocity norm over flood or ebb cycle within drifter window
    if debug:
        print 'calculating tide norm...'

    if isinstance(tide, (str, unicode)):
        if bounds:
            pos1=np.where((model.Grid.lon >= bounds[0]) & \
                          (model.Grid.lon <= bounds[1]))[0]
            pos2=np.where((model.Grid.lat >= bounds[2]) & \
                          (model.Grid.lat <= bounds[3]))[0]
            box = np.intersect1d(pos1, pos2)

            if tide == 'flood':
                tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,box,:], 0)
            elif tide == 'ebb':
                tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,box,:], 0)
        else:
            if tide == 'flood':
                tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
            elif tide == 'ebb':
                tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
        label='Mean Velocity Norm during '+tide.capitalize()+' Tide (m/s)'
    else:
        tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
        label='Mean Velocity Norm (m/s)'

   # create spatially varying color map of mean velocity norm
    if debug:
        print 'preparing to create colormap...'
    if not fname:
        title_main = "Trajectory for {} on {}".format(loc, date)
    else:
        title_main = "Trajectory for {}".format(fname)

    if len(tide) > 1:
        label='Mean Velocity Norm (m/s)'
    else:
        label='Mean Velocity Norm during '+tide[0]+' Tide (m/s)'

    fig = createColorMap(ncfile, tideNorm[0,:], mesh=False, \
            bounds=bounds, title=title_main, debug=debug, \
            label=label)

    # create title
    # fig.suptitle('Data from ' + fname[:-4], fontsize=14)
    x = drift.Variables.lon
    y = drift.Variables.lat
    u = drift.Variables.u
    v = drift.Variables.v

    if debug:
        print 'creating scatter plot...'
    plt.scatter(x,y)

    if debug:
        print 'preparing to plot time series...'

    # create validation structure
    if debug:
        print 'creating validation object...'

    try:
        valid = Validation(drift, ncfile, flow='sf', debug=False)
    except IndexError:
        print 'cannot create validation object for drifter %i.' % i
        return 0

    # calculate speed from interpolated and observed date
    mTimes = valid.Variables.struct['mod_time']
    oU = valid.Variables.struct['obs_timeseries']['u']
    oV = valid.Variables.struct['obs_timeseries']['v']
    mU = valid.Variables.struct['mod_timeseries']['u']
    mV = valid.Variables.struct['mod_timeseries']['v']

    if debug:
        print '\tcalculating speeds...'
    speedS = np.asarray(np.sqrt(mU**2 + mV**2))
    speedO = np.asarray(np.sqrt(oU**2 + oV**2))

    # ratio addition
    if debug:
        print '\tadding ratio adjustments...'
    speedS = speedS * ratio
    datetimes = np.asarray([dn2dt(time) for time in mTimes])

    # For now, separate the two plots.
    fig2=plt.figure()
    result = plotTimeSeries(fig2, np.reshape(np.tile(datetimes,2),\
            (2, len(datetimes))), np.vstack((speedS, speedO)), \
            loc, label=['Simulated','Observed'], where=111, \
            title=loc + ' Drifter Speeds for ' + date,  \
            ylab='Speed (m/s)')

    # if not result:
    #     if debug:
    #         print '...error encountered with drifter {}.'.format(i)
    #         print 'continuing...'
    #     plt.close()
    #     continue

    if plot:
        if debug:
            print 'displaying plot...'
        plt.show()
    if savepath:
        if debug:
            print 'saving plot...'
        fig.savefig(savepath + '_traj.png')
        result.savefig(savepath + '_ts.png')
        if debug:
            print '...plot saved to: ', savepath

    # clear the figure window
    # plt.close()
    return fig, result
