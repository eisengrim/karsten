import re

findVector = re.compile("\(?([+-]?[\d.]+),\s*([+-]?[\d.]+),\s*([+-]?[\d.]+)\)?")

def s2vector(string):
    """Convert a string representation of a vector - '(%f, %f, %f)' - to a list of floats"""
    v = findVector.match(string)
    if v:
        components = [float(x) for x in v.group(1, 2, 3)]
        return components
    raise ValueError("Input string (%s) is not in the format '(%%f, %%f, %%f)'" % string)

def v2string(vector):
    """Convert a vector into a string representation  - '(%f, %f, %f)'"""
    return "(%+.10f, %+.10f, %+.10f)" % (vector[0], vector[1], vector[2])


