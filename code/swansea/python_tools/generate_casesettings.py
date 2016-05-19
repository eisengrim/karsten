#!/usr/bin/env python2.7
import os
import CaseSettings
from vecutils import s2vector, v2string

## @file
## @brief Generate Case Settings file
## @ingroup executable

## @brief Accept user input in vector form
def vector_input(prompt):
    """
    Prompt for user input in vector form (%f, %f, %f)
    Retry on bad data
    """
    while True:
        try:
            vec = s2vector(raw_input("%s [vector]: " % str(prompt)))
            return vec
        except ValueError as e:
            print "Error: %s" % e

def runInteractive():
    defnthreshold = 3

    cs = CaseSettings.CaseSettings()
    cs.mode = 3
    cs.ptsteps = 400
    cs.zthreshold = 9
    cs.fish = 23
    cs.depth = 15
    cs.noise = 16
    cs.scaling = (100,100,0.000001)
    cs.basename = "/home/tom/NorthSea/processed/r3snsNk001Cor.slf"
    
    distributions = ["HP4","HP5"]
    foodweights = {"A": 0, "B": 0.1, "C": 1, "D": 10, "E": 100}
    noisemultis = {"A": 0, "B": 0.1, "C": 1, "D": 10, "E": 100}
    for dist in distributions:
        for fw in foodweights.keys():
            for nm in noisemultis.keys():
                cs.particles = "/home/tom/NorthSea/Thames%s.txt" % dist
                cs.outputpath = "/home/tom/NorthSea/parametric-study/results/%s-%s-%s" % (fw, nm, dist)
                try:
                    os.mkdir("/home/tom/NorthSea/parametric-study/results/%s-%s-%s" % (fw, nm, dist))
                except OSError as e:
                    pass

                cs.nthreshold = 3 * noisemultis[nm]
                cs.foodweight = foodweights[fw]
                cs.write("%s-%s-%s.ini" % (fw, nm, dist))

if __name__ == "__main__":
    runInteractive()
