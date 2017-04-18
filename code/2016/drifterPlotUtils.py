#! /usr/env/python2.7

import numpy as np
from pyseidon_dvt import *
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from matplotlib.dates import DateFormatter
from drifterUtils import dn2dt, checkIDs
from createColorMap import *
import seaborn as sns
from datetime import datetime
import os.path as osp
import os
import sys

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


def spatialError(model, lon, lat, uspdO, uspdS, debug=False):
    """
    Creates a spatially-varying plot of the errors.
    """
    if debug:
        print '\tstarting spatial plot...'
    glon = model.Grid.lon
    glat = model.Grid.lat
    if debug:
        print '\tcomputing bounding boxes...'
    bounds = [np.min(glon), np.max(glon), np.min(glat), np.max(glat)]

    if debug:
        print '\tcomputing colors...'

    speedErr = uspdO - uspdS

    if debug:
        print '\tcreating map...'
        print '\tcreating scatter plot...'
	f=plt.figure()
	ax = f.add_axes([.125,.1,.775,.8])
	ax.triplot(glon, glat, model.Grid.trinodes, zorder=10, lw=10)
	clim=np.percentile(speedErr,[5,95])
    if debug:
        print '\tcreating color bar...'

	cb = ax.scatter(lon, lat, c=speedErr, s=10, edgecolor='None', \
            vmin=clim[0], vmax=clim[1], zorder=20)
    plt.colorbar(cb)
    if debug:
        print '\tcreating color bar...'

    return f


def plotTimeSeries(fig, dt, var, loc, label=['',''], title='', bg='#f5deb3', \
                styles=['b--', '#8B0000'], ylab='', where=121, axx=None, \
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


def trajectoryPlots(ncfile, drift, loc, date, savepath=None, bfric='0.015', \
                fname=None, tight=False, ratio=1.0, debug=False, \
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

    if debug:
        print 'calculating ebb/flood split at centre of location...'
        print 'calculating model velocity norm...'

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
