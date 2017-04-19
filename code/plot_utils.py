#! /usr/env/python2.7

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from matplotlib.dates import DateFormatter
import matplotlib as mpl
import seaborn as sns
from datetime import datetime
import os.path as osp
import os
import sys
import json

# local imports
from utils import *
from color_map import createColorMap
from velo_norm import velo_norm


def varCorr(var1, var2, xlabel='', ylabel='', title='', where=111, \
        style=None, plot=False):
    """
    Calculates the correlation between two given variables, so long as they
    are the same size (Uses the Pearson product-moment correlation coefficient).
    If plot==True, a scatter plot of the two variables will be generated.

    Returns the Pearson coefficient and the two-tailed p-value.
    """
    fig = plt.Figure()
    if style is not Nont:
        s = json.load(open(style))
        mpl.rcParams.update(s)

    ax = fig.add_subplot(where)

    ax.scatter(var1, var2)
    if title:
        plt.title(title)
    if ylabel:
        plt.ylabel(ylabel)
    if xlabel:
        plt.xlabel(xlabel)
    plt.grid(True)

    print 'type \'plt.show()\' to display...'

    return pearsonr(var1, var2)


def spatialError(glon, glat, lon, lat, obs, sim, trinodes):
    """
    Creates a spatially-varying plot of the errors.
    """
    bounds = [np.min(glon), np.max(glon), np.min(glat), np.max(glat)]

    err = obs - sim

    f=plt.figure()
    ax = f.add_axes([.125,.1,.775,.8])
    ax.triplot(glon, glat, trinodes, zorder=10, lw=10)
    clim=np.percentile(err,[5,95])

    cb = ax.scatter(lon, lat, c=err, s=10, edgecolor='None', \
        vmin=clim[0], vmax=clim[1], zorder=20)
    plt.colorbar(cb)

    return f


def plotTimeSeries(dt, var, loc, label=['',''], title='', \
                style=None, ylab='', where=111, lines=['--','-'],\
                axx=None, axy=None, legend=True):
    """
    Creates a comparative var vs. time graph from simulated and observed
    arrays. This function is also passed an existing figure object to plot on.

    input:
        - datetimes : two datetimes
        - vars : two variables to plot
        - loc : location tag
        - axx : x axis to share
        - axy : y axis to share
        - label : label for legend
        - title : title of plot
        - style : a json filename for mpl styling
        - lines : linestyles
        - ylab : label of var axis
        - legend : boolean. add legend.
        - where : location of plot on fig (3 digit int)
    return:
        - fig, ax
    """
    fig = plt.figure()

    if style is not None:
        s = json.load(open(style))
        mpl.rcParams.update(s)
    # add subplot and configure axes
    if axx or axy:
        if axx:
            ax = fig.add_subplot(where, sharex=axx)
        if axy:
            ax = fig.add_subplot(where, sharey=axy)
    else:
        ax = fig.add_subplot(where)

    try:
        ax.plot(dt[0], var[0], lines[0], label=label[0])
        ax.plot(dt[1], var[1], lines[1], label=label[1])
    except ValueError:
        return False

    ax.set_ylabel(ylab)
    ax.set_xlabel('Time (HH:MM:SS)')
    ax.set_title(title)
    if legend:
        leg = ax.legend(loc='best')
        for lab in leg.get_texts():
            lab.set_fontsize(12)
        for lab in leg.get_lines():
            lab.set_linewidth(1)

    # set the axis limits (hardcoded for consistency)
    # ax2.set_ylim(0.0, 3.0)
    fomt = DateFormatter('%H:%M:%S')
    plt.gca().xaxis_date()
    plt.gca().get_xaxis().set_major_formatter(fomt)
    plt.gcf().autofmt_xdate()
    plt.grid(True)
    plt.tight_layout()

    plt.hold('on')

    return fig


def trajectoryPlot(u, v, lon, lat, mtime, x, y, otime, trinodes, \
                vel_norm=None, label=None, title_main=None, bounds=[]):
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
    if not vel_norm:
        vel_norm = velo_norm(u, v)

    # creates drifter object window for flow map
    win1 = (np.abs(mtime-otime.min())).argmin()
    win2 = (np.abs(mtime-otime.max())).argmin()

    # averages velocity norm over flood or ebb cycle within drifter window
    tideNorm = np.mean(vel_norm[win1:win2,:,:], 0)

    # create spatially varying color map of mean velocity norm
    fig = createColorMap(tideNorm[0,:], lon, lat, trinodes, \
            mesh=False, bounds=bounds, title=title_main, \
            label=label, where=111)

    plt.scatter(x,y)

    return fig
