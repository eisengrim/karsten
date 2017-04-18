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
from velo_norm import *

def trajectoryPlots(u, v, lon, lat, mtime, x, y, otime, drift, savepath, \
                title_main=None, plot=False, bounds=[]):
    """
    Compiles necessary data, creates plots and saves / shows them all. The plot is
    a spatially varying color map of speed with a drifter trajectory superposed and
    a comparitive speed-time plot for each drifter in the given directory.

    input:
        - bounds : array of min max lon lat *priority over tight*
    returns:
        - fig1, a trajectory, fig2, a timeseries
    """
    sns.set(font="serif")

    # computes variable norm
    vel_norm = velo_norm(u, v)

    if debug:
        print 'creating plots...'

    if savepath:
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
    tModel = mtime
    tDrift = otime
    win1 = (np.abs(tModel-tDrift.min())).argmin()
    win2 = (np.abs(tModel-tDrift.max())).argmin()

    # averages velocity norm over flood or ebb cycle within drifter window
    if bounds:
        pos1=np.where((lon >= bounds[0]) & \
                          (lon <= bounds[1]))[0]
        pos2=np.where((lat >= bounds[2]) & \
                          (lat <= bounds[3]))[0]
        box = np.intersect1d(pos1, pos2)

        tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,box,:], 0)
    else:
        tideNorm = np.mean(ncfile.Variables.velo_norm[win1:win2,:,:], 0)
    label='Mean Velocity Norm (m/s)'

    # create spatially varying color map of mean velocity norm
    fig = createColorMap(ncfile, tideNorm[0,:], mesh=False, \
            bounds=bounds, title=title_main, debug=debug, \
            label=label, where=111)

    plt.scatter(x,y)

    return fig
