#!/usr/bin/env python
import random
import sys
import math

def genRandom(coordinates, nid):
    return random.random()

def linear(points, values, coordinates, nid):
    v = 0
    for i in range(len(points)):
        dist = math.sqrt((coordinates[nid][0] - points[i][0]) ** 2 + (coordinates[nid][1] - points[i][1]) ** 2)
        v += float(values[i])/dist
    return v

def generateData(conf, nodefile, varnum, varname, force, datagen=genRandom):
    timesteps = 0
    with open("%s.times.txt" % conf.basename, "r") as timestepfile:
        timesteps = int(timestepfile.readline())

    if varname is None:
        varname = os.path.splitext(os.path.basename(nodefile))[0]

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

    with open("%s.vars.txt" % conf.basename, "w+") as varfile:
        varfile.write("%d\n" % numvars)
        varfile.writelines(["%d\t%s\n" % (x, varnames[x].rstrip()) for x in range(numvars)])


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



    values = []
    if nodefile:
        for x in range(numnodes):
            values.append(0)

        with open(nodefile, "r") as nodes:
            for nid in nodes.readlines():
                values[int(nid)] = datagen(coords, nid)
    else:
        for nid in range(numnodes):
            values.append(datagen(coords, nid))

    for ts in range(timesteps):
        with open("%s.var%d.t%d.txt" % (conf.basename, int(varnum), ts), "w") as tsfile:
            tsfile.writelines(["%d\t%f\n" % (nid, values[nid]) for nid in range(numnodes)])

if __name__ == "__main__":
    import argparse

    cmdparser = argparse.ArgumentParser(description="Populate random values at designated nodes and add to simulation")
    cmdparser.add_argument("-f", "--force", action="store_true", default=False, dest="force", help="Overwrite existing data.")
    cmdparser.add_argument("-n", "--number", action="store", default=None, dest="varnum", help="Write to specified variable number.", metavar="vnum")
    cmdparser.add_argument("-l", "--label", action="store", default=None, dest="varname", help="Variable name", metavar="vname")

    nodeselectors = cmdparser.add_mutually_exclusive_group(required=True)
    nodeselectors.add_argument("-i", "--ids", action="store", dest="nodefile", help="File containing node IDs")
    nodeselectors.add_argument("-a", "--all", action="store_true", dest="allnodes", help="Generate values at all nodes")

    cmdparser.add_argument("-m", "--mode", action="store", default="genRandom", dest="mode", help="Specify function for data generation")
    cmdparser.add_argument("-p", "--point", action="append", dest="points", help="Specify point sources (for applicable generators)")
    cmdparser.add_argument("-v", "--value", action="append", dest="values", help="Values for point sources (for applicable generators)")

    cmdparser.add_argument("config", action="store", help="Case/Settings file")

    args = cmdparser.parse_args()

    if not args.config:
        print "Case/Settings file must be specified"
        parser.print_help()
        sys.exit(1)

    from CaseSettings import CaseSettings
    conf = CaseSettings(args.config)
    if args.mode == "random" or args.mode == "genRandom":
        mode = genRandom
    else:
        rmode = locals()[args.mode]
        points = []
        values = []
        from ParticleFile import s2vector
        if args.points is None:
            args.points = ["0,0,0"]
            args.values = [1]

        for p in args.points:
            points.append(s2vector(p))

        if len(args.values) == 0:
            for i in points:
                values.append(1)
        elif len(args.values) == 1:
            for i in points:
                values.append(args.values[0])
        else :
            for v in args.values:
                values.append(v)

        def psource(coords, nid):
            return rmode(points, values, coords, nid)

        mode = psource

    if args.allnodes:
        args.nodefile = None
    generateData(conf, args.nodefile, int(args.varnum), args.varname, args.force, mode)
