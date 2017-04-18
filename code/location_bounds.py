#! /usr/env/python2.7

def get_bounds(loc):
    # find the location centre for flood/tide split calculation
    if loc == 'GP':
        centre = [-66.33906, 44.26898]
        if tight and not bounds:
            bounds = [-66.355, -66.31, 44.245, 44.2925]
        elif not bounds:
            bounds = []
    elif loc == 'DG':
        centre = [-65.76000, 44.67751]
        if tight and not bounds:
            bounds = [-65.775, -65.77, 44.665, 44.69]
        elif not bounds:
            bounds = []
    elif loc == 'PP':
        centre = [-66.206924, 44.389368]
        # find out the tightness required for PP
        if tight and not bounds:
            bounds = [-66.225, -66.195, 44.37, 44.41]
        elif not bounds:
            bounds = []
    elif 'MP' in loc:
        centre = [-64.418942, 45.3499696]
        if tight and not bounds:
            bounds = [-64.51, -64.31, 45.31, 45.38]
        elif not bounds:
            bounds = []

return centre, bounds
