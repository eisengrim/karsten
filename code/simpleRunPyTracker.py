from pyticle_tracker import pyticle
import os
import numpy as np
import time
import scipy.io as sio
import matplotlib.pyplot as plt

loc = 'HH'
sim = '2011_Oct_22-29_AF_3D'
filename = '/EcoII/acadia_uni/projects/acadia_force_numerical_model_r20/2011-10-22_2011-10-29/output/acadia_force_3d_0001.nc'
outpath = '/EcoII/acadia_uni/projects/drifters/pyticle_tracker/' \
        + 'HH_20111022-29/'+loc+'_'+sim+'_output.nc'

inlocs = None

options={}
options['starttime']=0
options['endtime']=-2
options['interpolationratio']=10
options['outputratio']=2
options['ncformat']='NETCDF4_CLASSIC'
options['useLL']=True
options['randomStart']=True
options['centre']=(-64.625361, 45.206072)
options['radius']=500
options['np'] = 10
options['layer']=0
options['gridDim']='3D'

options['diffusion']=False
# projection for dngridCSR
#options['projstr']='lcc +lon_0=-64.55880 +lat_0=41.78504 ' + \
#            '+lat_1=39.69152 +lat_2=43.87856'

# projection for acadia force
options['projstr']='lcc +lon_0=-64.55880 +lat_0=41.84493 +lat_1=39.72147 '\
        + '+lat_2=43.96838'

start = time.clock()
mypy=pyticle(filename, inlocs, outpath, options=options)
mypy.run()
print('run in: %f' % (time.clock() - start))



