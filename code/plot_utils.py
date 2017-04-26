#! /usr/env/python2.7

import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import matplotlib.ticker as tic
import matplotlib.tri as Tri
import seaborn as sns
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


def varCorr(var1, var2, xlabel='', ylabel='', title='', style=None, where=111):
    """
    Calculates the correlation between two given variables, so long as they
    are the same size (uses the Pearson product-moment correlation coefficient).
    If plot==True, a scatter plot of the two variables will be generated.

    Returns the Pearson coefficient and the two-tailed p-value.

    var1, var2  : two variable arrays, same length
    xlabel, ylabel, title : self-explanatory
    style : json style filename
    where : where on figure to place plot
    """
    sns.set(font="serif")
    fig = plt.Figure()
    if style is not None:
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


def spatialError(glon, glat, lon, lat, obs, sim, trinodes, where=111, \
        label=None, title=None, figsize=(18,10), hide=False, error="sgn"):
    """
    Creates a spatially-varying plot of the errors between two spatially
    varying data.

    Type of error measurement: "abs", "sgn", "rel", "pct",
    """
    bounds = [np.min(glon), np.max(glon), np.min(glat), np.max(glat)]

    if error is "abs":
        err = np.abs(sim - obs)
    elif error is "bias":
        err = sim - obs
    elif error is "rel":
        err = np.abs(sim - obs) / sim
    elif error is "sgn":
        err = (sim - obs) / sim
    elif error is "pct":
        err = 100 * np.abs(sim - obs) / sim

    # define figure window
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(where, aspect=(1.0/np.cos(np.mean(glat) * np.pi/180.0)))

    # create triangulation object
    tri = Tri.Triangulation(glon, glat, triangles=trinodes)

    # setting limits and levels of colormap
    cmin = err[:].min()
    cmax = err[:].max()
    step = (cmax-cmin) / 50.0

    # depth contours to plot
    levels = np.arange(cmin, (cmax+step), step)

    # triangular grid
    f = ax.tripcolor(tri, np.zeros(len(glon[:])), cmap=plt.cm.PRGn)
    plt.triplot(glon, glat, trinodes, zorder=10, lw=0.25, color='white')

    # scatter plot and color bar
    clim = np.percentile(err,[5,95])
    cb = ax.scatter(lon, lat, c=err, s=20, edgecolor='None', \
        vmin=clim[0], vmax=clim[1], zorder=20, cmap=plt.cm.hot)
    cbar = fig.colorbar(cb, ax=ax)
    cbar.set_label(label, rotation=-90, labelpad=30)

    # label and axis parameters
    if hide:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    else:
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    ax.patch.set_facecolor('0.5')
    scale = 1

    # ticker for coordinate degree axis
    if hide:
        ax.set_yticklabels([])
        ax.set_xticklabels([])
    else:
        ticks = tic.FuncFormatter(lambda lon, pos: '{0:g}'.format(lon/scale))
        ax.xaxis.set_major_formatter(ticks)
        ax.yaxis.set_major_formatter(ticks)

    ax.set_xlim([bounds[0], bounds[1]])
    ax.set_ylim([bounds[2], bounds[3]])
    ax.grid()
    plt.title(title)

    plt.hold('on')

    return f


def plotTimeSeries(dt, var, loc, sd=None, iqr=False, label=['',''], title='', \
                style=None, ylab='', where=111, lines=['--','-'],\
                axx=None, axy=None, legend=True, multi=False):
    """
    Creates a comparative var vs. time graph from simulated and observed
    arrays. This function is also passed an existing figure object to plot on.

    input:
        - datetimes : two datetimes
        - vars : two variables to plot
        - sd : vector of sdevs for first var (optional)
        - iqr : boolean. plots boxplots instead for var 1 rather than errorbars
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
    # format date, set style
    fomt = DateFormatter('%H:%M:%S')

    if style is not None:
        s = json.load(open(style))
        mpl.rcParams.update(s)

    # deal with breaks in x axis
    if multi:
        offset = np.nonzero(dt[0] - dt[0][0])[0][0]
        step = (dt[0][offset] - dt[0][0]).seconds

        # find the breaks in time
        dtt = np.asarray([d.seconds for d in np.diff(dt[0])])
        breaks = np.squeeze(np.where(dtt > 4*step))
        num = len(breaks) + 1
        breaks = np.insert(breaks, 0, -1)
        breaks = np.append(breaks, -1)

        # generate subplots based on # breaks, plots them
        fig, ax = plt.subplots(1, num, sharey=True, tight_layout=True)
        fig.subplots_adjust(wspace=0.005)
        for k in xrange(num):
            try:
                if sd is not None:
                    ax[k].errorbar(dt[0], var[0], yerr=sd, fmt=lines[0], label=label[0])
                elif iqr:
                    pass
                else:
                    ax[k].plot(dt[0], var[0], lines[0], label=label[0])
                ax[k].plot(dt[1], var[1], lines[1], label=label[1])
            except ValueError:
                return False

            # format axes to show dates
            ax[k].set_xlim(dt[0][breaks[k]+1], dt[0][breaks[k+1]])
            ax[k].xaxis_date()
            ax[k].get_xaxis().set_major_formatter(fomt)

            if k == 0:
                ax[k].spines['right'].set_visible(False)
                ax[k].tick_params(right='off')
            # elif k == (num - 1):
            #     ax[k].spines['left'].set_visible(False)
            else:
                ax[k].spines['right'].set_visible(False)
                ax[k].spines['left'].set_visible(False)
                ax[k].tick_params(right='off')
                ax[k].tick_params(left='off')

            # rotate dates 30 degreres
            ticks = ax[k].get_xticklabels()
            for tick in ticks:
                tick.set_rotation(30)

        # set labels and legends
        ax[0].set_ylabel(ylab)
        ax[num/2].set_xlabel('Time (HH:MM:SS)')
        ax[num/2].set_title(title)
        if legend:
            leg = ax[num - 1].legend(loc='best')
            for lab in leg.get_texts():
                lab.set_fontsize(12)
            for lab in leg.get_lines():
                lab.set_linewidth(1)

        plt.gcf().autofmt_xdate()
        plt.grid(True)
        plt.tight_layout()

    # if not multi, dont bother with temporal breaks
    else:
        # add subplot and configure axes
        fig = plt.figure()

        if axx or axy:
            if axx:
                ax = fig.add_subplot(where, sharex=axx)
            if axy:
                ax = fig.add_subplot(where, sharey=axy)
        else:
            ax = fig.add_subplot(where)

        # generate plots
        try:
            if sd is not None:
                ax[k].errorbar(dt[0], var[0], yerr=sd, fmt=lines[0], label=label[0])
            elif iqr:
                pass
            else:
                ax.plot(dt[0], var[0], lines[0], label=label[0])
            ax.plot(dt[1], var[1], lines[1], label=label[1])
        except ValueError:
            return False

        # set labels, legend
        ax.set_ylabel(ylab)
        ax.set_xlabel('Time (HH:MM:SS)')
        ax.set_title(title)
        if legend:
            leg = ax.legend(loc='best')
            for lab in leg.get_texts():
                lab.set_fontsize(12)
            for lab in leg.get_lines():
                lab.set_linewidth(1)

        # set the axis limits
        fomt = DateFormatter('%H:%M:%S')
        plt.gca().xaxis_date()
        plt.gca().get_xaxis().set_major_formatter(fomt)
        plt.gcf().autofmt_xdate()
        plt.grid(True)
        plt.tight_layout()

    plt.hold('on')

    return fig


def trajectoryPlot(u, v, lon, lat, mtime, x, y, otime, trinodes, hide=False, \
                vel_norm=None, label=None, title_main=None, bounds=[]):
    """
    Compiles necessary data, creates plots and saves / shows them all. The plot is
    a spatially varying color map of speed with a drifter trajectory superposed and
    a comparitive speed-time plot for each drifter in the given directory.

    input:
        - u, v      : u and v model velocity arrays
        - lon, lat  : model lon and lat
        - mtime     : model date times
        - x, y      : drifter lon, lat arrays
        - otime     : observed date times
        - trinodes  : model trinodes
        - hide      : hide long/lat on plots, use hidden location
        - vel_norm  : use a particular vel_norm
        - label     : label for color bar
        - title_main: title for plot
        - bounds    : array of min max lon lat *priority over tight*
    returns:
        - fig
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
            label=label, where=111, hide=hide)

    plt.scatter(x,y,c='#9b009b')

    return fig

