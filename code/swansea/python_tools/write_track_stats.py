#!/bin/env python
# coding: utf-8
import os
import sys
import argparse
import pandas as pd
from read_track_file import *


cmdparser = argparse.ArgumentParser(description="Extract statistics for a collection of particle tracks")
cmdparser.add_argument("filename", help="CSV file produced by extractparticle", nargs='+')
cmdparser.add_argument("-k", "--keep-together", action="store_true", default=False, dest="keeptogether",
        help="Save plots adjacent to original track files.")

args = cmdparser.parse_args()

for p in args.filename:
    print "Reading file %s...." % p
    filename = p

    if not args.keeptogether:
        filename = os.path.normpath(filename).replace(os.sep, "_")

    filename = "%s.stats.txt" % os.path.splitext(filename)[0]
    print "\tOutputting stats to %s..." % filename
    tracks = pd.read_csv(p, delimiter="\t", skipinitialspace=True, skiprows=1, names=track_cols)#, usecols=[0,2,3,4,5,6,7])
    stats = tracks.describe()
    with open(filename, "w") as text_file:
            text_file.write(stats.to_string())
    del tracks
    del stats
