#! /usr/env/python2.7
import numexpr as ne

def velo_norm(u, v, w=None):
    """
    This method computes a velocity norm given 3 cartesian velocities.
    """
    #Check if w if there
    if w is not None:
        try:
            vel = ne.evaluate('sqrt(u**2 + v**2 + w**2)').squeeze()
        except MemoryError as e:
            print 'data too large for machine memory or server'
            raise
    else:
        try:
            # computing velocity norm
            vel = ne.evaluate('sqrt(u**2 + v**2)').squeeze()
        except MemoryError as e:
            print 'data too large for machine memory'
            raise

    return vel

