#!/usr/bin/env python
"""
VTKExport.py

Export mesh inputs used by floaters/prepmesh in a format usable by ParaView
"""

import argparse

cmdparser = argparse.ArgumentParser(description="Export mesh to VTK format")
cmdparser.add_argument("config", help="Case/Settings file used with floaters/prepmesh")
cmdparser.add_argument("-s", "--suffix", action="store", default="export", dest="suffix", help="Add suffix to generated file names")
cmdparser.add_argument("-t", "--time", action="store", default=None, dest="times", nargs="+", help="Only include specified time steps", metavar="t")

args=cmdparser.parse_args()


from CaseSettings import CaseSettings
conf = CaseSettings(args.config)

if args.times is None:
    args.times = conf.getTimes()
try:
    _ = iter(args.times)
except TypeError:
    args.times = [args.times]

import numpy
from evtk.hl import pointsToVTK

x = numpy.loadtxt("%s.x.txt" % conf.basename, usecols=[1], skiprows=1)
y = numpy.loadtxt("%s.y.txt" % conf.basename, usecols=[1], skiprows=1)
for time in [t[0] for t in args.times]:
    data = {}
    for v in conf.getVariables():
        d = numpy.loadtxt("%s.var%d.t%d.txt" % (conf.basename, v[0], time), usecols=[1])
        data[v[1]] = d

    if conf.mode == 2:
        z = numpy.zeros_like(x)
    elif conf.mode == 3:
        z = numpy.loadtxt("%s.var%d.t%d.txt" % (conf.basename, conf.z, time), usecols=[1])

    pointsToVTK("%s/%s-%s.t%d" % (conf.outputpath, conf.basefilename, args.suffix, time), x, y, z, data=data)
