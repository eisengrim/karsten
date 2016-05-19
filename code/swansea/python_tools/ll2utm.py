#!/usr/bin/env python2

from math import pi, sqrt, sin, cos, tan, floor

def ll2utm(lat, lon):
    """ Convert WGS84 Latitude/Longitude pairs to UTM Eastings/Northings """
    a = 6378137
    b = 6356752.3
    e = sqrt(1 - pow(b/a, 2))
    e2 = pow(e, 2) / (1 - pow(e, 2))
    k0 = 0.9996
    n = (a - b) / (a + b)

    LonZone = 31 + floor(lon/6)
    LatRad = (lat * pi) / 180
    LonRad = (lon - (6 * LonZone - 183)) * pi / 180

    A0 = a*(1-n+(5*pow(n, 2)/4)*(1-n)+(81*pow(n, 4)/64)*(1-n))
    B0 = (3*a*n/2)*(1-n-(7*pow(n, 2)/8)*(1-n)+55*pow(n, 2)/64)
    C0 = (15*a*pow(n, 2)/16)*(1-n+(3*pow(n, 2)/4)*(1-n))
    D0 = (35*a*pow(n, 3)/48)*(1-n+11*pow(n, 2)/16)
    E0 = (315*a*pow(n, 4)/51)*(1-n)
    MeridonalArc = A0 * LatRad - B0 * sin(2 * LatRad) + C0 * sin(4 * LatRad) \
            - D0 * sin(6 * LatRad)+E0 * sin(8 * LatRad)

    T = a / sqrt(1 - e * pow(sin(LatRad), 2))
    ScaleMArc = MeridonalArc * k0

    Y = T * sin(LatRad)*cos(LatRad/2)
    Z = (T * sin(LatRad)*pow(cos(LatRad), 3)) / 24 \
            * (5 - pow(tan(LatRad), 2) + 9 * e2 * pow(cos(LatRad), 2) \
            + 4 * pow(e2, 2) * pow(cos(LatRad), 4)) * k0

    Northing = ScaleMArc + Y * pow(LonRad, 2) + Z * pow(LonRad, 4)

    AA = T * cos(LatRad) * k0
    AB = pow(cos(LatRad), 3) * T / 6 * k0 \
            * (1 - pow(tan(LatRad), 2) + e2 * pow(cos(LatRad), 2))
    Easting = 500000 + AA * LonRad + AB * pow(LonRad, 3)

    return (Easting, Northing, LonZone)

if __name__ == "__main__":
    import argparse
    import sys

    cmdparser = argparse.ArgumentParser(description=ll2utm.__doc__)
    cmdparser.add_argument("-f", "--filename", action="store", default=None, dest="filename", help="CSV filename")
    cmdparser.add_argument("-s", "--skip-header", action="store_true", dest="skiphead", help="CSV file contains header row")
    cmdparser.add_argument("-p", "--point", action="store", default=None, dest="point", help="Convert single point from decimal lat lon format", metavar="(lat, lon)")

    args = cmdparser.parse_args()
    if (args.filename and args.point) or not (args.filename or args.point):
        print "Either a filename or point (not both) must be provided"
        sys.exit(1)

    if args.point:
        import numpy
        ar = numpy.fromstring(args.point, sep=",")
        E, N, LZ = ll2utm(ar[0], ar[1])

        print "Easting:\t%f\nNorthing:\t%f\nLon. Zone:\t%d" % (E, N, LZ)
        sys.exit(0)

    if args.filename:
        import os
        import numpy

        header = ""
        if args.skiphead:
            f = open(args.filename, "r")
            h = f.readline()
            hr = h.split(",")
            hr = hr[2:]
            nh = ["Easting", "Northing", "LonZone"]
            nh.extend(hr)
            header = ",".join(nh).rstrip()
            f.close()

        dat = numpy.loadtxt(args.filename, delimiter=',', skiprows=(1 if args.skiphead else 0))

        nd = []
        for row in dat:
            a = ll2utm(row[0], row[1])
            a = list(a)
            a.extend(numpy.delete(row, [0, 1], 0))
            nd.append(a)

        numpy.savetxt("%s/UTM-%s" % (os.path.dirname(args.filename), os.path.basename(args.filename)), nd, fmt="%10f", delimiter=",", header=(header if args.skiphead else ""))
