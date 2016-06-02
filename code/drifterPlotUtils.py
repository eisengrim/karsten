#! /usr/env/python2.7

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

from drifterUtils import dn2dt


def plotBias(speedS, ubias, bfric, debug=False):
    """
    Plots the bias information. Returns the R-squared.
    """
    if debug:
        print 'creating plot...'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(speedS, ubias, alpha=0.25)
    ax.set_xlabel('Model Speed (m/s)')
    ax.set_ylabel('Bias')
    ax.set_title('Bias vs. Model Speed for BFRIC={}'.format(bfric))

    # determine line of best fit
    if debug:
        print '\tdetermining line of best fit...'

    par = np.polyfit(speedS, ubias, 1)
    m = par[-2]
    b = par [-1]

    variance = np.var(ubias)
    residuals = np.var([(m*xx + b - yy)  for xx,yy \
            in zip(speedS, ubias)])
    Rsqr = np.round(1-residuals/variance, decimals=5)
    if debug:
        print '\tR^2 value for bias plot is {}...'.format(Rsqr)
    plt.hold('on')
    ax.plot(speedS, m*speedS+b, 'r-')
    plt.grid(True)
    print 'type \'plt.show()\' to display...'

    return Rsqr


def plotCubeSpeeds(speedS, speedO, bfric, debug=False):
    """
    Creates a plot of the speeds cubed, and returns the R-squared.
    """
    fi = plt.figure()
    ax1 = fi.add_subplot(111)
    speedS3 = np.power(speedS, 3)
    speedO3 = np.power(speedO, 3)
    ax1.scatter(speedS3, speedO3, alpha=0.25)
    ax1.set_xlabel('Model Speed (m/s)')
    ax1.set_ylabel('Drifter Speed (m/s)')
    ax1.set_title('Model and Drifter Speed Comparison for BFRIC={}'.format(bfric))
    coeff = np.polyfit(speedS3, speedO3, 1)
    m = coeff[-2]
    b = coeff[-1]
    if debug:
        print '\tcoeffs for cube plot are: \n\t\tm={}\n\t\tb={}'.format(m,b)
    plt.hold('on')
    ax1.plot(speedS3, m*speedS3+b)
    variance = np.var(speedO3)
    residuals = np.var([(m*xx + b - yy)  for xx,yy \
            in zip(speedS3, speedO3)])
    Rsqr = np.round(1-residuals/variance, decimals=5)
    if debug:
        print '\tR^2 for cube plot is {}'.format(Rsqr)
    print 'type \'plt.show()\' to display...'

    return Rsqr


def plotBiasvDrifter(bias, num_drift, mean, bfric, debug=False):
    """
    Creates a plot of the bias vs. the drifters.
    """
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    if debug:
        print '\tnum drift is ', np.arange(1,num_drift+1).size
        print '\tlen of mean is ', len(mean)

    # hacky fix for differing lengths
    try:
        ax2.plot(np.arange(1,num_drift+1), mean, 'go')
    except:
        ax2.plot(np.arange(1,num_drift), mean, 'go')
    ax2.set_ylabel('Bias')
    ax2.set_xlabel('Drifter')
    ax2.set_title('Individual Bias vs Drifter for BFRIC={}'.format(bfric))
    plt.hold('on')
    ax2.axhline(y=np.mean(bias), linewidth=2)
    plt.grid(True)

    print 'type \'plt.show()\' to display...'


def plotBiasvDepth(drift, debug=False):
    """
    Plots error vs. depth for a single drifter.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if debug:
        print 'creating plot...'

    bias = drift['bias']
    depth = drift['depth']

    ax.plot(depth, bias, 'bo')
    ax.set_ylabel('Bias')
    ax.set_xlabel('Depth along Trajectory (m)')
    ax.set_title('Bias vs. Depth for {}'.format(drift['filename']))
    plt.grid(True)
    print 'type \'plt.show()\' to display...'


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


def spatialRatios(model, lon, lat, uspdO, uspdS, debug=False):
    """
    Creates a spatially-varying plot of the power ratios.
    """
    if debug:
        print '\tstarting spatial plot...'
    glon = model.Grid.lon
    glat = model.Grid.lat
    if debug:
        print '\tcomputing bounding boxes...'
    bounds = [np.min(glon), np.max(glon), np.min(glat), np.max(glat)]

    uspeedO3 = np.power(uspdO, 3)
    uspeedS3 = np.power(uspdS, 3)
    if debug:
        print '\tcomputing colors...'
    var3 = np.divide(uspeedO3, uspeedS3)
    # var3 = np.divide(uspeedO3, uspeedS3) # does this produce a different plot?
    # color = np.subtract(var3, np.min(var3)) / (np.max(var3) - np.max(var3))

    if debug:
        print '\tcreating map...'
        print '\tcreating scatter plot...'
	f=plt.figure()
	ax = f.add_axes([.125,.1,.775,.8])
	ax.triplot(glon, glat, model.Grid.trinodes, zorder=10, lw=10)
	clim=np.percentile(var3,[5,95])
	cb = ax.scatter(lon, lat, c=var3, s=10, edgecolor='None', \
            vmin=clim[0], vmax=clim[1], zorder=20)
    plt.colorbar(cb)
    if debug:
        print '\tcreating color bar...'


def plotTimeSeries(fig, dt, var, loc, label=['',''], title='', bg='#f5deb3', \
                styles=['b--', '#8B0000'], axis_label='', where=121, axx=None, \
                axy=None, legend=True, debug=False):
    """
    Creates a comparative var vs. time graph from a model object and a drifter
    object. This function is also passed an existing figure object to plot on.
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
        - axis_label : label of var axis
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

    ax.set_ylabel(axis_label)
    ax.set_xlabel('Time (HH:MM:SS)')
    ax.set_title(title)
    if legend:
        leg = ax.legend(loc='best')
        for lab in leg.get_texts():
            lab.set_fontsize('small')
        for lab in leg.get_lines():
            lab.set_linewidth(1)

    # set the axis limits (hardcoded for consistency)
    # ax2.set_ylim(0.0, 3.0)
    plt.gcf().autofmt_xdate()
    plt.grid(True)
    plt.tight_layout()

    plt.hold('on')

    if debug:
        print '...time series successfully plotted.'

    return fig, ax
