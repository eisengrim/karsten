import numpy as np
import sys, re
import fileinput as fin

if __name__ == '__main__':

    pfile = sys.argv[1]
    plog = sys.argv[2]
    outfile = sys.argv[3]

    if "Particle #" in open(plog).read():
        parts = re.findall(r'#(\w+)', open(plog).read())

    for i in parts:
        for line in open(plog):
