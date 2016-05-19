#!/usr/bin/env python2.7

import os
import sys
import argparse
import numpy as np
from matplotlib import use
use("GDK")
from matplotlib import pyplot
import pandas as pd
import seaborn as sb


def plot_orientations(data, graphfile):
    try:
        axes=data.plot(y=[u'Yaw', u'Pitch', u'Roll'])
        axes.legend(["Yaw", "Pitch", "Roll"])
        figure = axes.get_figure()
        figure.suptitle("Particle Orientation angles")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_orientations]"
        print e

def plot_forces(data, graphfile):
    try:
        axes=data.plot(y=[u'Force_X', u'Force_Y', u'Force_Z'])
        axes.legend(("F_x", "F_y", "F_z"))
        figure = axes.get_figure()
        figure.suptitle("Particle forces")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_forces]"
        print e

def plot_coords(data, graphfile):
    try:
        axes = data.plot(y=[u'Pos_X', u'Pos_Y', u'Pos_Z'])
        axes.legend(("x", "y", "z"))
        figure = axes.get_figure()
        figure.suptitle("Particle Coordinates")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_coords]"
        print e

def plot_velocities(data, graphfile):
    try:
        axes = data.plot(y=[u'Vel_X', u'Vel_Y', u'Vel_Z'])
        axes.legend(("V_x", "V_y", "V_z"))
        figure = axes.get_figure()
        figure.suptitle("Particle velocities")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_velocities]"
        print e

def plot_2dpos(data, graphfile, xlim=None, ylim=None):
    try:
        axes = data.plot(x='Pos_X', y='Pos_Y')

        if xlim is None:
            axes.set_xlim(xmin=data.Pos_X.min(), xmax=data.Pos_X.max(), auto=0)
        else:
            axes.set_xlim(xmin=xlim[0], xmax=xlim[1], auto=0)

        if ylim is None:
            axes.set_ylim(ymin=data.Pos_Y.min(), ymax=data.Pos_Y.max(), auto=0)
        else:
            axes.set_ylim(ymin=ylim[0], ymax=ylim[1], auto=0)

        axes.legend(("Path", "Start", "End"))
        figure = axes.get_figure()
        figure.suptitle("Particle position - horizontal movements")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_2dpos]"
        print e

def plot_vertpos(data, graphfile):
    try:
        axes = data.plot(y=u'Pos_Z')
        axes.legend("Vertical Position")
        figure = axes.get_figure()
        figure.suptitle("Particle vertical position")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_vertpos]"
        print e

def plot_targetvec(data, graphfile):
    try:
        axes = data.plot(y=['Target_X', 'Target_Y', 'Target_Z'])
        axes.legend(("T_x", "T_y", "T_z"))
        figure = axes.get_figure()
        figure.suptitle("Target Vector Components")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_targetvec]"
        print e

def plot_orand(data, graphfile):
    try:
        axes = data.plot(y='O_Rand')
        figure = axes.get_figure()
        figure.suptitle("Orientation Random factor")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_orand]"
        print e

def plot_vrand(data, graphfile):
    try:
        axes = data.plot(y='V_Rand')
        figure = axes.get_figure()
        figure.suptitle("Velocity Random factor")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_vrand]"
        print e

def plot_targetvec_xy(data, graphfile):
    try:
        cmap = sb.cubehelix_palette(as_cmap=True, reverse=True)
        if (data['Target_X'].min() == data['Target_X'].max() or data['Target_Y'].min() == data['Target_Y'].max()):
            print "There is no variation in Target vector - not plotting"
            return -1
        axes = pyplot.quiver(data['Pos_X'], data['Pos_Y'], data['Target_X'], data['Target_Y'], data['Timestamp'], cmap=cmap)
        figure = axes.get_figure()
        figure.suptitle("Target Vector in 2D")
        figure.savefig(graphfile, dpi=200)
        pyplot.close(figure)
    except Exception as e:
        print "Unable to plot data [in plot_targetvec_xy]"
        print e

def plot_all(data, base):
    plot_orientations(data, "%s.orientations.png" % base)
    plot_forces(data, "%s.forces.png" % base)
    plot_coords(data, "%s.coords.png" % base)
    plot_velocities(data, "%s.velocities.png" % base)
    plot_2dpos(data, "%s.2dpos.png" % base)
    plot_vertpos(data, "%s.vertpos.png" % base)
    plot_targetvec(data, "%s.targetvec.png" % base)
    plot_orand(data, "%s.orand.png" % base)
    plot_vrand(data, "%s.vrand.png" % base)
    plot_targetvec_xy(data, "%s.targetvec_xy.png" % base)

if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="Plot details from a CSV file produced by extractparticle.\nSeems to work better using the GDK backend")
    cmdparser.add_argument("filename", help="CSV file produced by extractparticle")
    cmdparser.add_argument("-o", "--output", action="store", default=None, dest="output_dir", help="Output to directory 'foo'", metavar="foo")
    cmdparser.add_argument("-x", "--xkcdify", action="store_true", default=False, dest="xkcd", help="Output XKCD style plots")
    cmdparser.add_argument("-t", "--trim", action="store_true", default=False, dest="trim", help="Drop duplicates from data")

    args = cmdparser.parse_args()
    args.filename = os.path.abspath(args.filename)

    if args.xkcd:
        pyplot.xkcd()

    if args.output_dir:
        os.chdir(args.output_dir)

    print "Outputting plots to %s" % os.getcwd()

    np.seterr('raise')
    data = pd.read_csv(args.filename, quotechar="'", skipinitialspace=True, index_col="Time")
    if args.trim:
        data = data.drop_duplicates(subset=data.columns.drop("Timestamp"))

    base = os.path.basename(args.filename)
    base = os.path.splitext(base)[0]
    plot_all(data, base)
