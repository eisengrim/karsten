from pyticle_tracker import pyticle
import numpy as np
import time



ncin='/media/moe46/runs/vhhigh_v2/vhhigh_v2_2012-02-01_2012-03-01/output/vhhigh_v2_0001.nc'
lldict=np.load('vhhigh_v2_600k_12kmx5km.npy')
lldict=lldict[()]

locations=np.vstack([lldict['lon'],lldict['lat']]).T
print(locations.shape)

options={}
options['starttime']=1349
options['endtime']=480+options['starttime']
options['interpolationratio']=15
options['outputratio']=10
options['ncformat']='NETCDF3_64BIT'
options['zlib']=False
options['useLL']=True
options['layer']='da'
options['gridDim']='2D'


#start = time.clock()
#mypy=pyticle(ncin,locations,'output/starttime_1349.nc',options=options)
#mypy.run()
#print('run in: %f' % (time.clock() - start))


start = time.clock()
options['starttime']=1353
mypy=pyticle(ncin,locations,'output/starttime_1353.nc',options=options)
mypy.run()
print('run in: %f' % (time.clock() - start))


start = time.clock()
options['starttime']=1357
mypy=pyticle(ncin,locations,'output/starttime_1357.nc',options=options)
mypy.run()
print('run in: %f' % (time.clock() - start))


start = time.clock()
options['starttime']=1361
mypy=pyticle(ncin,locations,'output/starttime_1361.nc',options=options)
mypy.run()
print('run in: %f' % (time.clock() - start))


start = time.clock()
options['starttime']=1365
mypy=pyticle(ncin,locations,'output/starttime_1365.nc',options=options)
mypy.run()
print('run in: %f' % (time.clock() - start))


start = time.clock()
options['starttime']=1369
mypy=pyticle(ncin,locations,'output/starttime_1369.nc',options=options)
mypy.run()
print('run in: %f' % (time.clock() - start))
