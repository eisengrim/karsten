#!/usr/env/python2.7

# lib imports
from __future__ import division
from pyticle_tracker import pyticle
import sys, os
import os.path as osp
import scipy as sp
import numpy as np
import time
import scipy.io as sio
import netCDF4 as nc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import matplotlib.tri as Tri
from matplotlib.dates import DateFormatter
from matplotlib import rc
from datetime import datetime, timedelta
from pyseidon import *
import argparse as arp
from sklearn.preprocessing import normalize
import pandas as pd
import pylab

# local imports
from createColorMap import createColorMap
from drifterUtils import *
from drifterPlotUtils import plotTimeSeries
from plotParticles import *


def plotTracks(pytkl, drift, ncfile, dname, tideNorm, sim, loc, save=False, \
        write=False, verbose=False, saveplot=None):
    """
    Plots the tracks of the two objects.
    """
    # check whether or not to overwrite
    savename = dname + '_plot.png'

    if save and not write:
        if osp.exists(saveplot + savename):
             return

    if verbose:
        print 'preparing to create color map...'

    # hacky fix for odd tide norm behaviour when window splicing is the same
    if len(tideNorm.shape) > 1:
        tideNorm = tideNorm[0,:]

    fig = createColorMap(ncfile, tideNorm[:], mesh=False, \
            title = 'Particle Track for ' + dname,
                        label = 'Mean Velocity Norm (m/s)')

    # add scatter plots of data
    lon = pytkl.variables['lon'][:]
    lat = pytkl.variables['lat'][:]
    if verbose:
        print 'adding scatter plot...'
    pytkl_path = plt.scatter(lon, lat, c='m', lw=0, alpha=0.6, marker="^", s=10)

    lonD = drift.Variables.lon
    latD = drift.Variables.lat
    drift_path = plt.scatter(lonD, latD, c='k', lw=0, s=25)

    plt.legend([drift_path, pytkl_path],['Drifter', 'Model'])

    if save:
        if verbose:
            print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def geographicError(pytkl, drift, ncfile, dname, sim, loc, save=False, \
            write=False, verbose=False, saveplot=None, tn=None):
    """
    Calculates errors in lon/lat coordinates and velocities graphically.

    pytkl :: pyticle object
    ncfile :: FVCOM mode
    drift :: drifter data
    sim :: sim date
    loc :: sim location
    dname :: drifter file name
    save :: bool, save plot
    write :: bool, overwrite plots
    verbose :: bool, verbosity
    saveplot :: name of save plot dir
    tn :: optional tideNorm, prints side by side
    """
    # check whether or not to overwrite
    savename = dname + '_lonlatErr.png'

    if save and not write:
        if osp.exists(saveplot + savename):
             return

    if verbose:
        print 'calculating lon/lat error...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat
    uD = drift.Variables.u
    vD = drift.Variables.v
    u = np.array(pytkl.variables['u'][:])
    v = np.array(pytkl.variables['v'][:])

    dtime = drift.Variables.matlabTime
    ptime = mjd2num(pytkl.variables['time'][:])

    if verbose:
        print str(lon.shape[1]) + ' pyticles...'
        print 'no. pts pytkl: ' + str(len(lon))
        print 'no. pts drift: ' + str(len(lonD))
        print 'broadcasting...'

    # truncate longer list to broadcast
    if len(lon) > len(lonD):
        idx = broadcast(ptime, dtime)
        lon = lon[idx]
        lat = lat[idx]
        u = u[idx]
        v = v[idx]
        ptime = ptime[idx]
    elif len(lonD) > len(lon):
        idx = broadcast(dtime, ptime)
        dtime = dtime[idx]
        lonD = lonD[idx]
        latD = latD[idx]
        uD = uD[idx]
        vD = vD[idx]

    dtime = [dn2dt(x) for x in dtime]
    ptime = [dn2dt(x) for x in ptime]

    # if multiple particles, take the mean
    if len(lon.shape) > 1:
        lonErr = np.squeeze(np.array([np.abs(row-lonD) for row in lon.T])).T
        latErr = np.squeeze(np.array([np.abs(row-latD) for row in lat.T])).T
        uErr = np.squeeze(np.array([np.abs(row-uD) for row in u.T])).T
        vErr = np.squeeze(np.array([np.abs(row-vD) for row in v.T])).T

        lonDiff = np.mean(lonErr, axis=1)
        latDiff = np.mean(latErr, axis=1)
        uDiff = np.mean(uErr, axis=1)
        vDiff = np.mean(vErr, axis=1)
        lab1 = ['Mean Lon. Diff.', 'Mean Lat. Diff.']
        lab2 = ['Mean u Diff.', 'Mean v Diff.']
    else:
        lonDiff = np.abs(lonD - lon)
        latDiff = np.abs(latD - lat)
        uDiff = np.abs(uD - u)
        vDiff = np.abs(vD - v)
        lab1 = ['Lon. Diff.', 'Lat. Diff.']
        lab2 = ['u Diff.', 'v Diff.']

    if tn != None:
        ## **Needs work done to ensure proper scaling of subplots on figure!
        where1 = 221
        where2 = 223
        fig = createColorMap(ncfile, tn[0,:], mesh=False, \
                title='Trajectory for {}'.format(dname), \
                label='Mean Velocity Norm (m/s)')
        fig.suptitle('Positional Timeseries :: ' + \
                 'BFRIC={} | Filename={}'.format(str(bfric), dname),
                 fontsize=10)
        plt.scatter(lon, lat)
        plt.scatter(lonD, latD)
    else:
        fig = plt.figure()
        where1 = 211
        where2 = 212

    if verbose:
        print 'creating figures...'
    result, axis = plotTimeSeries(fig, (ptime, dtime), (lonDiff, latDiff), loc, \
            label=lab1, where=111, \
            title='Longitudal/Latitudal Timeseries Error', \
            axis_label='Lon/Lat Error ($^\circ$)', styles=['#FFC300','#581845'],\
            debug=True, legend=True)

    # styles=['#1DA742','#900C3F'],

    if not result:
        sys.exit('error plotting longitudal/latitudal data.')

    # result, axis = plotTimeSeries(fig, (ptime, dtime), (uDiff, vDiff), loc, \
    #         label=lab2, where=where2, \
    #         title='u/v-Velocity Timeseries Error', \
    #         axis_label='u/v-Velocity Error (m/s)', styles=['#FFC300','#581845'], \
    #         axx=axis, debug=True, legend=True)

    # if not result:
    #     sys.exit('error plotting velocity data.')

##    result, axis = plotTimeSeries(fig, (ptime, dtime), (lon, lonD), loc, \
##            label=lab, where=where1, \
##            title='Longitudal Timeseries Comparison', \
##            axis_label='Longitude ($^\circ$)', styles=['#1DA742','#900C3F'], \
##            debug=True, legend=leg)
##
##    if not result:
##        sys.exit('error plotting longitudal data.')
##
##    result, axis = plotTimeSeries(fig, (ptime, dtime), (lat, latD), loc, \
##            label=['Model Drifter','Observed Drifter'], where=where2, \
##            title='Latitudal Timeseries Comparison', \
##            axis_label='Latitude ($^\circ$)', styles=['#FFC300','#581845'], \
##            axx=axis, debug=True, legend=leg)
##
##    if not result:
##        sys.exit('error plotting latitudal data.')

    if save:
        if verbose:
            print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def boxplotError(pytkl, drift, ncfile, dname, sim, loc, save=False, write=False, \
            verbose=False, saveplot=None):
    """
    Plots a box plot for each model-drifter time series error.
    """
    # check whether or not to overwrite
    savename = dname + '_lonlatErr.png'

    if save and not write:
        if osp.exists(saveplot + savename):
             return

    if verbose:
        print 'calculating lon/lat error...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat
    uD = drift.Variables.u
    vD = drift.Variables.v
    u = np.array(pytkl.variables['u'][:])
    v = np.array(pytkl.variables['v'][:])

    dtime = drift.Variables.matlabTime
    ptime = mjd2num(pytkl.variables['time'][:])

    if verbose:
        print str(lon.shape[1]) + ' pyticles...'
        print 'no. pts pytkl: ' + str(len(lon))
        print 'no. pts drift: ' + str(len(lonD))
        print 'broadcasting...'

    # truncate longer list to broadcast
    if len(lon) > len(lonD):
        idx = broadcast(ptime, dtime)
        lon = lon[idx]
        lat = lat[idx]
        u = u[idx]
        v = v[idx]
        ptime = ptime[idx]
    elif len(lonD) > len(lon):
        idx = broadcast(dtime, ptime)
        dtime = dtime[idx]
        lonD = lonD[idx]
        latD = latD[idx]
        uD = uD[idx]
        vD = vD[idx]

    dtime = [dn2dt(x) for x in dtime]
    ptime = [dn2dt(x) for x in ptime]

    if len(lon.shape) > 1:
        lonErr = np.squeeze(np.array([np.abs(row-lonD) for row in lon.T])).T
        latErr = np.squeeze(np.array([np.abs(row-latD) for row in lat.T])).T
        uErr = np.squeeze(np.array([np.abs(row-uD) for row in u.T])).T
        vErr = np.squeeze(np.array([np.abs(row-vD) for row in v.T])).T

    if verbose:
        print 'creating figure...'
    # PANDAS? PYLAB?
    # Is the error normal??
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    bp = ax1.boxplot(lonErr, showmeans=True, patch_artist=True)
    for box in bp['boxes']:
        box.set(color='#1C2833',linewidth=2)
        box.set(facecolor='lightblue', alpha=1)
    for whisker in bp['whiskers']:
        whisker.set(color='#1C2833', linewidth=2)
    for cap in bp['caps']:
        cap.set(color='#1C2833', linewidth=2)
    for median in bp['medians']:
        median.set(color='#1C2833')
    for flier in bp['fliers']:
        flier.set(marker='+')
        flier.set(color='k', alpha=0.5)
    for mean in bp['means']:
        mean.set(marker='o')
        mean.set(color='k')
    ax1.set_ylabel('Longitudal Error')
    ax1.set_title('Particle-Drifter Error')

    ax2 = fig.add_subplot(212)
    bp = ax2.boxplot(lonErr, showmeans=True, patch_artist=True)
    for box in bp['boxes']:
        box.set(color='#1C2833',linewidth=2)
        box.set(facecolor='lightblue', alpha=1)
    for whisker in bp['whiskers']:
        whisker.set(color='#1C2833', linewidth=2)
    for cap in bp['caps']:
        cap.set(color='#1C2833', linewidth=2)
    for median in bp['medians']:
        median.set(color='#1C2833')
    for flier in bp['fliers']:
        flier.set(marker='+')
        flier.set(color='k', alpha=0.5)
    for mean in bp['means']:
        mean.set(marker='o')
        mean.set(color='k')
    ax2.set_xlabel('Particle Label')
    ax2.set_ylabel('Latitudal Error')

    if save:
        if verbose:
            print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def spatialProbability(pytkl, drift, ncfile, dname, sim, loc, save=False, \
                bfric=0.015, write=False, verbose=False, saveplot=None):
    """
    Creates a color-varying spatial plot of the model drifters.
    """
    # check whether or not to overwrite
    savename = dname + '_prob.png'

    if save and not write:
        if osp.exists(saveplot + savename):
             return

    if verbose:
        print 'calculating domain frequencies...'
    elems = np.array(pytkl.variables['indomain'][:])

    # calculate spatial probability density
    # count the number of occurrences of elements
    freq = np.zeros(ncfile.Grid.nele)
    idx, cnt = np.unique(elems, return_counts=True)
    freq[idx] = np.log(cnt)

    # freq = freq.normalize(freq[:,np.newaxis], axis=0).ravel()
    fig = createColorMap(ncfile, freq/np.amax(freq), mesh=False, \
            label = "Log Probability of Model Drifter Transit",
            title = "Spatially-Varying Probability Plot of Model Drifter\n" + \
            'BFRIC={} | Filename={}'.format(str(bfric), dname))

    plt.scatter(drift.Variables.lon, drift.Variables.lat, c='k', lw=0)

    if save:
        if verbose:
            print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)

    else:
        plt.show()

    plt.close(fig)


def pyticleStats(pytkl, drift, dname, verbose=False):
    """
    Generates statistics of the error for run.
    """

    if verbose:
        print 'collecting lon/lats...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat

    if verbose:
        print str(lon.shape[1]) + ' pyticles...'
        print 'no. pts pytkl: ' + str(len(lon))
        print 'no. pts drift: ' + str(len(lonD))

    if verbose:
        print 'calculating statistics...'
    lonM = np.mean(lon, axis=1)
    latM = np.mean(lat, axis=1)
    lonSD = np.std(lon, axis=1)
    latSD = np.std(lat, axis=1)

    dlonM = lonM - lonM[0]
    dlatM = latM - latM[0]

    disp = np.sqrt(np.square(dlonM*79868) + np.square(dlatM*111117))
    dispSD = np.sqrt(np.square(lonSD*79868) + np.square(latSD*111117))

    t = sp.stats.t.ppf(0.95, len(lonM)-1)
    CI = t*dispSD

    error_lon = np.sqrt(((lonD - lonM) ** 2).mean())
    error_lat = np.sqrt(((latD - latM) ** 2).mean())

    print error_lon, error_lat


def dispersionPlots(pytkl, drift, ncfile, dname, sim, loc, save=False, \
                write=False, verbose=False, saveplot=None):
    """
    Plots mean long./lat. and relative displacement of the pyticles
    against time, with a 95% confidence interval.
    """
    # check whether or not to overwrite
    savename = dname + '_disp.png'

    if save and not write:
        if osp.exists(saveplot + savename):
             return

    if verbose:
        print 'collecting lon/lats...'
    lon = np.array(pytkl.variables['lon'][:])
    lat = np.array(pytkl.variables['lat'][:])
    lonD = drift.Variables.lon
    latD = drift.Variables.lat

    dtime = drift.Variables.matlabTime
    ptime = mjd2num(pytkl.variables['time'][:])
    dtime = [dn2dt(x) for x in dtime]
    ptime = [dn2dt(x) for x in ptime]

    if verbose:
        print str(lon.shape[1]) + ' pyticles...'
        print 'no. pts pytkl: ' + str(len(lon))
        print 'no. pts drift: ' + str(len(lonD))

    if verbose:
        print 'calculating statistics...'
    lonM = np.mean(lon, axis=1)
    latM = np.mean(lat, axis=1)
    lonSD = np.std(lon, axis=1)
    latSD = np.std(lat, axis=1)

    dlonM = lonM - lonM[0]
    dlatM = latM - latM[0]

    disp = np.sqrt(np.square(dlonM*79868) + np.square(dlatM*111117))
    dispSD = np.sqrt(np.square(lonSD*79868) + np.square(latSD*111117))

    t = sp.stats.t.ppf(0.95, len(lonM)-1)
    CI = t*dispSD

    # plotting
    if verbose:
        print 'creating plot...'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ptime, disp,'o', color='#34495E', markersize=5, markeredgewidth=1,\
            markeredgecolor='#2E4053', markerfacecolor='#34495E')
    ax.fill_between(ptime, disp+CI, disp-CI, color='#F39C12', edgecolor='')
    # minor hack for labelling CI
    ax.plot(ptime, disp+CI, color='#F39C12', label='95% Confidence Interval')
    ax.set_ylabel('Displacement (m)')
    ax.set_xlabel('Time (HH:MM:SS)')
    ax.set_title('Dispersion of Particles')
    fomt = DateFormatter('%H:%M:%S')
    plt.gca().xaxis_date()
    plt.gca().get_xaxis().set_major_formatter(fomt)
    plt.gcf().autofmt_xdate()
    plt.grid(True)
    plt.tight_layout()
    plt.ylim(0, np.max(disp)+np.max(dispSD))

    if save:
        if verbose:
            print 'creating save directory...'
        if not osp.exists(saveplot):
             os.makedirs(saveplot)
        plt.savefig(saveplot + savename)
    else:
        plt.show()

    plt.close(fig)


