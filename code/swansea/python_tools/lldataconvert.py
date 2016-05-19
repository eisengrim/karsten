#!/usr/bin/env python
import random
import sys
import math


def generateData(conf, varnum, varname, force, datagen):
    timesteps = 0
    print "Reading times file..."
    with open("%s.times.txt" % conf.basename, "r") as timestepfile:
        timesteps = int(timestepfile.readline())

    print "%d timesteps to generate" % timesteps
    tsoutfreq = pow(10, max(1, int(math.log(timesteps, 10)) - 1))

    if varname is None:
        print "Variable name required."
        return

    with open("%s.vars.txt" % conf.basename, "r") as varfile:
        numvars = int(varfile.readline())
        varnames = varfile.readlines()
        varnames = [x.split("\t")[1] for x in varnames]
        print "%d names read from vars file (%d variables declared)" % (len(varnames), numvars)
        if len(varnames) != numvars:
            print "Vars file appears to be inconsistent - aborting."
            sys.exit(1)

        if varnum > numvars:
            print "Specified varnum (%d) too large, setting to next valid number (%d)" % (varnum, numvars)
            varnum = None


        if varnum is None:
            varnum = numvars

        if varnum == numvars:
            numvars += 1
            varnames.append(args.varname)

        elif varnum < numvars:
            if not force:
                print "Varnum < number of variables and force not specified - use force flag if you really want to clobber the existing data (%s)" % varnames[varnum].rstrip()
                sys.exit(1)
            else:
                varnames[varnum] = varname

    print "Updating variables file to include %s..." % varname

    with open("%s.vars.txt" % conf.basename, "w+") as varfile:
        varfile.write("%d\n" % numvars)
        varfile.writelines(["%d\t%s\n" % (x, varnames[x].rstrip()) for x in range(numvars)])


    print "Reading coordinate data..."
    coords = None 
    with open("%s.x.txt" % conf.basename, "r") as meshx:
        numnodes = int(meshx.readline())
        print "%d nodes" % numnodes
        coords = [[0, 0] for x in range(numnodes)]
        for line in meshx:
            (node, sep, xc) = line.rpartition("\t")
            coords[int(node)][0] = float(xc)

    with open("%s.y.txt" % conf.basename, "r") as meshy:
        nn2 = int(meshy.readline())
        if not nn2 == numnodes:
            print "Inconsistent mesh!"
            sys.exit(1)
        for line in meshy:
            (node, sep, yc) = line.rpartition("\t")
            coords[int(node)][1] = float(yc)



    print "Coordinates read. Calculating values at nodes..."
    values = []
    noutfreq = pow(10, max(0, int(math.log(numnodes, 10)) - 1))
    for nid in range(numnodes):
        if (nid % noutfreq == 0):
            print "Node %d / %d" % (nid, numnodes)
        values.append(datagen(coords, nid))

    print "Outputting data..."
    for ts in range(timesteps):
        if (ts % tsoutfreq == 0):
            print "Writing file %d / %d" % (ts, timesteps)

        with open("%s.var%d.t%d.txt" % (conf.basename, int(varnum), ts), "w") as tsfile:
            tsfile.writelines(["%d\t%f\n" % (nid, values[nid]) for nid in range(numnodes)])

def linear(points, values, coordinates, nid):
    v = 0
    for i in range(len(points)):
        dist = math.sqrt((coordinates[nid][0] - points[i][0]) ** 2 + (coordinates[nid][1] - points[i][1]) ** 2)
        v += float(values[i])/dist
    return v

def sublinear(points, values, coordinates, nid):
    v = 0
    for i in range(len(points)):
        dist = math.sqrt((coordinates[nid][0] - points[i][0]) ** 2 + (coordinates[nid][1] - points[i][1]) ** 2)
        v += float(values[i])/pow(dist, 0.25)
    return v

def invsq(points, values, coordinates, nid):
    v = 0
    for i in range(len(points)):
        dist = math.sqrt((coordinates[nid][0] - points[i][0]) ** 2 + (coordinates[nid][1] - points[i][1]) ** 2)
        v += float(values[i])/(dist ** 2)
    return v

if __name__ == "__main__":
    import argparse

    cmdparser = argparse.ArgumentParser(description="Populate random values at designated nodes and add to simulation")
    cmdparser.add_argument("-f", "--force", action="store_true", default=False, dest="force", help="Overwrite existing data.")
    cmdparser.add_argument("-n", "--number", action="store", default=None, dest="varnum", help="Write to specified variable number.", metavar="vnum")
    cmdparser.add_argument("-l", "--label", action="store", default=None, dest="varname", help="Variable name", metavar="vname")

    cmdparser.add_argument("config", action="store", help="Case/Settings file")
    cmdparser.add_argument("csv", action="store", help="CSV file (lat, lon, val)")

    args = cmdparser.parse_args()

    if not args.config:
        print "Case/Settings file must be specified"
        parser.print_help()
        sys.exit(1)

    if not args.csv:
        print "CSV file must be specified"
        parser.print_help()
        sys.exit(1)

    from CaseSettings import CaseSettings
    conf = CaseSettings(args.config)

    import numpy
    dat = numpy.loadtxt(args.csv, delimiter=",", skiprows=1)

    from ll2utm import ll2utm
    ps=[]
    vs=[]
    for row in dat:
        E,N,LZ = ll2utm(row[0], row[1])
        ps.append([E,N])
        vs.append(row[2])

    def generator(coords, nid):
        return sublinear(ps, vs, coords, nid)
    generateData(conf, int(args.varnum), args.varname, args.force, generator)
