from pyticle_tracker import pyticle
import os
import numpy as np
import time
import scipy.io as sio
import matplotlib.pyplot as plt

loc = 'PP'
bfric = '0.015'
sim = '2014_Aug_12_3D'

filename = '/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/' + \
            'drifter_runs/BFRIC_'+bfric+'/'+loc+'/'+sim+'/output/subdomain_'+ \
            loc+'1_0001.nc'
outpath = '/array/home/119865c/karsten/pytrkr/'+loc+'_'+sim+'_output.nc'
start_locs = '/array/home/119865c/karsten/drifters/start_info/init_locs_'+\
            loc+'_'+sim+'.dat'

inlocs=np.genfromtxt(start_locs)

print(inlocs.shape)

options={}
options['starttime']=10
options['endtime']=35
options['interpolationratio']=2
options['outputratio']=1
options['ncformat']='NETCDF4_CLASSIC'
options['useLL']=True
options['layer']=0
options['gridDim']='2D'
options['projstr']='lcc +lon_0=-64.55880 +lat_0=41.78504 +lat_1=39.69152 +lat_2=43.87856'

start = time.clock()
mypy=pyticle(filename, inlocs, outpath, options=options)
mypy.run()
print('run in: %f' % (time.clock() - start))



