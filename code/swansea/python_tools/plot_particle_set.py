#!/usr/bin/env python2.7

import os
import sys
import argparse
import pandas as pd
import plot_extracted_particle as pep

cmdparser = argparse.ArgumentParser(description="Plot a collection of particle tracks")
cmdparser.add_argument("filename", help="CSV file produced by extractparticle", nargs='+')
cmdparser.add_argument("-t", "--trim", action="store_true", default=False, dest="trim",
        help="Drop duplicated results from tracks\nNB: Will result in plots representing variable durations")
cmdparser.add_argument("-k", "--keep-together", action="store_true", default=False, dest="keeptogether",
        help="Save plots adjacent to original CSV files.")
cmdparser.add_argument("-gx", "--xlimit", dest="xlim", type=float, nargs=2, default=None, help="X axis limits for all plots")
cmdparser.add_argument("-gy", "--ylimit", dest="ylim", type=float, nargs=2, default=None, help="Y axis limits for all plots")

args = cmdparser.parse_args()

if args.xlim and not len(args.xlim) == 2:
    print "X limit must be specified as min max"
    sys.exit(1)

if args.ylim and not len(args.ylim) == 2:
    print "Y limit must be specified as min max"
    sys.exit(1)

for p in args.filename:
    print "Reading file %s...." % p
    data = pd.read_csv(p, quotechar="'", skipinitialspace=True, index_col="Time")
    if args.trim:
        data = data.drop_duplicates(subset=data.columns.drop("Timestamp"))
    filename = p

    if not args.keeptogether:
        filename = os.path.normpath(filename).replace(os.sep, "_")

    filename = os.path.splitext(filename)[0]
    print "\tPlotting data to %s..." % filename
    pep.plot_2dpos(data, "%s.2dpos.png" % filename, args.xlim, args.ylim)
    del data
