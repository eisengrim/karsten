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

def varCorr(var1, var2, xlabel='', ylabel='', title='', debug=False, plot=False):
    """
    Calculates the correlation between two given variables, so long as they
    are the same size (Uses the Pearson product-moment correlation coefficient).
    If plot==True, a scatter plot of the two variables will be generated.

    Returns the Pearson coefficient and the two-tailed p-value.
    """
    fig = plt.Figure()
    ax = fig.add_subplot(111)

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


def spatialError(glon, glat, lon, lat, obs, sim, trinodes, debug=False):
    """
    Creates a spatially-varying plot of the errors.
    """
    if debug:
        print '\tstarting spatial plot...'
    if debug:
        print '\tcomputing bounding boxes...'
    bounds = [np.min(glon), np.max(glon), np.min(glat), np.max(glat)]

    if debug:
        print '\tcomputing colors...'

    err = obs - sim

    if debug:
        print '\tcreating map...'
        print '\tcreating scatter plot...'
	f=plt.figure()
	ax = f.add_axes([.125,.1,.775,.8])
	ax.triplot(glon, glat, trinodes, zorder=10, lw=10)
	clim=np.percentile(err,[5,95])
    if debug:
        print '\tcreating color bar...'

	cb = ax.scatter(lon, lat, c=err, s=10, edgecolor='None', \
            vmin=clim[0], vmax=clim[1], zorder=20)
    plt.colorbar(cb)
    if debug:
        print '\tcreating color bar...'

    return f


def plotTimeSeries(fig, dt, var, loc, label=['',''], title='', bg='#f5deb3', \
                styles=['b--', '#8B0000'], ylab='', where=121, axx=None, \
                axy=None, legend=True, debug=False):
    """
    Creates a comparative var vs. time graph from simulated and observed
    arrays. This function is also passed an existing figure object to plot on.

    input:
        - fig : figure object
        - datetimes : two datetimes
        - vars : two variables to plot
        - loc : location tag
        - axx : x axis to share
        - axy : y axis to share
        - label : label for legend
        - title : title of plot
        - bg : hex code for background
        - styles : line styles
        - ylab : label of var axis
        - legend : boolean. add legend.
        - where : location of plot on fig (3 digit int)
    return:
        - fig, ax
    """
    if debug:
        print '\tcreating subplot...'
        print '\tconfiguring axes...'

    # add subplot and configure axes
    if axx or axy:
        if axx:
            ax = fig.add_subplot(where, axisbg=bg, sharex=axx)
        if axy:
            ax = fig.add_subplot(where, axisbg=bg, sharey=axy)
    else:
        ax = fig.add_subplot(where, axisbg=bg)

    if debug:
        print 'shapes are var: ({}, {}) t: ({}, {})...'.format(len(var[0]), \
                    len(var[1]), len(dt[0]), len(dt[1]))

    if len(var[0]) < 5 or len(var[1]) < 5:
        return False
    # if a value error is encountered due to the data in pyseidon,
    # do not plot, ignore and move on...
    try:
        ax.plot(dt[0], var[0], styles[0], label=label[0], linewidth=2)
        ax.plot(dt[1], var[1], styles[1], label=label[1], linewidth=2)
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

    if debug:
        print '...time series successfully plotted.'

    return fig



