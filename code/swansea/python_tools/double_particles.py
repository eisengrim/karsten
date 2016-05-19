#!/usr/bin/python
""" Read an existing particle file and write out a new version with twice the number of particles"""
import ParticleFile as PF

from random import uniform

import argparse
import copy

cmdparser = argparse.ArgumentParser(description="Double up an existing particle file, adding a random offset to existing positions")
cmdparser.add_argument("-r", help="Range", default=50, dest="range")
cmdparser.add_argument("-o", help="Output file", default=None, required=True, dest="output")
cmdparser.add_argument("input", help="Existing particle definition file")

args = cmdparser.parse_args()

pfile = PF.ParticleFile()
pfile.read(args.input)

inparts = len(pfile.particles)
outparts = 2*inparts

print "Reading %d particles and preparing %d for output..." % (inparts, outparts)
pcopy = copy.deepcopy(pfile.particles)
for i in pcopy:
    p = i
    for x in [0, 1]:
        p.position[x] = uniform(i.position[x] - args.range/2.0, i.position[x] + args.range/2.0)
    p.id = p.id + inparts
    pfile.particles.append(p)

print "Particles generated, writing to %s..." % args.output

pfile.comment = "[Doubled] %s" % pfile.comment
pfile.np = outparts
pfile.write(args.output)

print "Done"
