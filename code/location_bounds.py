#! /usr/env/python2.7

def get_centre(loc):
    """
    get centres of certain domains (used for velo_norm)
    """
    # find the location centre for flood/tide split calculation
    if loc == 'GP':
        return [-66.33906, 44.26898]
    elif loc == 'DG':
        return [-65.76000, 44.67751]
    elif loc == 'PP':
        return [-66.206924, 44.389368]
    elif loc == 'MP':
        return [-64.418942, 45.3499696]


def get_bounds(loc):
    """
    get location bounds for tight plotting
    """
    if loc == 'GP':
        return [-66.355, -66.31, 44.245, 44.2925]
    elif loc == 'DG':
        return [-65.775, -65.77, 44.665, 44.69]
    elif loc == 'PP':
        return [-66.225, -66.195, 44.37, 44.41]
    elif 'MP' in loc:
        return [-64.51, -64.31, 45.31, 45.38]
