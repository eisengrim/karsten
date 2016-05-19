#!/usr/bin/env python2.7
import os
import ParticleFile as PF
from random import uniform

## @file
## @brief Generate particle definition files
## @ingroup executable

## @brief Accept user input in vector form
def vector_input(prompt):
    """
    Prompt for user input in vector form (%f, %f, %f)
    Retry on bad data
    """
    while True:
        try:
            vec = PF.s2vector(raw_input("%s [vector]: " % str(prompt)))
            return vec
        except ValueError as e:
            print "Error: %s" % e

def runInteractive():
    """
    Generate a particle definition file

    Creates a number of particles at a given location, with uniformly distributed
    positions and velocities within given limits
    """
    import readline
    import argparse

    cmdparser = argparse.ArgumentParser(description="Generate particle definition file")
    cmdparser.add_argument("-o", help="Output file", default=None, dest="filename")
    cmdparser.add_argument("-v", "--verbose", action="store_true", default=False, dest="verbose", help="More verbose output")

    args = cmdparser.parse_args()

    filename = args.filename

    if (not filename):
        filename = raw_input("Output file: ")
        filename = os.path.expandvars(os.path.expanduser(filename))

    filename = os.path.normpath(filename)
    if os.path.isfile(filename):
        confirm = raw_input("File '%s' already exists. Overwrite? [Y/N]: " % filename)
        confirm = confirm.lower()
        if not confirm in ["y", "yes"]:
            print "Output file exists and should not be overwritten. Exiting"
            return False

    pfile = PF.ParticleFile()

    while True:
        try:
            pfile.np = int(raw_input("Number of particles: "))
            break
        except ValueError as e:
            print "Error: %s" % e

    pfile.default = PF.Particle()
    pfile.comment = raw_input("Comment/Description: ")

    print ""
    print "[vector] in the prompt indicates a vector quantity"
    print "Input these as comma seperated decimals - e.g (0.8, -0.2, 2300.1)"
    print ""
    print "Particles will be distributed randomly in an area centred on the release coordinates"
    print "The release range determines the distance to the boundary of this area from the release coordinates in each direction"
    print "i.e a range of (5, 8, 2) from coordinates (0, 0, 0) would allow particles to be distributed between +5 and -5 in the x-direction, +8 and -8 in the y-direction and +2 and -2 in the z-direction"

    try:
        releaseCoords = vector_input("Release Coordinates")
        releaseRange = vector_input("Release Range")
        pfile.default.position = releaseCoords

        releaseVelocity = vector_input("Release Velocity")
        releaseVRange = vector_input("Velocity Range")
        #print "The centre of pressure is defined as an offset relative to the centre of mass, measured in metres"
        #CoP = vector_input("Centre of Pressure")
        pfile.default.cent_pressure = [0, 0, 0]

        pfile.default.orientation = vector_input("Orientation")
        orientationRange = vector_input("Orientation range")

        pfile.default.coeff_drag = vector_input("Drag Coefficient")
        #pfile.default.coeff_lift = float(raw_input("Lift Coefficient: "))
        pfile.default.coeff_lift = 0
        pfile.default.area = vector_input("Projected Area")

        while True:
            try:
                pfile.default.mass = float(raw_input("Mass: "))
                break
            except ValueError as e:
                print "Error: %s" % e

    except (EOFError, KeyboardInterrupt):
        print "Abandoning input at user request"
        return False

    pfile.particles = []
    for i in range(0, pfile.np):
        t = PF.Particle(pfile.default)
        t.id = i
        for x in [0, 1, 2]:
            t.position[x] = uniform(releaseCoords[x]-releaseRange[x], releaseCoords[x]+releaseRange[x])
            t.velocity[x] = uniform(releaseVelocity[x]-releaseVRange[x], releaseVelocity[x]+releaseVRange[x])
            t.orientation[x] = uniform(pfile.default.orientation[x]-orientationRange[x], pfile.default.orientation[x]+orientationRange[x])

        pfile.particles.append(t)

    fd = open(filename, "w")
    pfile.write(fd)
    fd.close()

if __name__ == "__main__":
    runInteractive()
