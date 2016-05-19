from pyticle_tracker import pyticle
import os
import numpy as np
import time
import scipy.io as sio
import matplotlib.pyplot as plt

os.system('rm test.nc')



#x=np.linspace(-66,-67,100)
#y=np.linspace(44,45,100)

#locx,locy=np.meshgrid(x,y)
#locx=locx.flatten()
#locy=locy.flatten()

#locations=np.vstack([locx,locy]).T

#print("Loading savelag1")
#fileload=h5.File('/home/moe46/workspace_matlab/lagtracker/savedir/vh_high/2012-02-01_2012-02-05/vh_high_test_large_20000pp_s0.mat')
#savelag1={}
#for i in fileload['savelag'].keys():
    #if (i=='x' or i=='y' or i=='z' or i=='time'):
        ##savelag1[i]=fileload['savelag'][i].value  
        #savelag1[i]=fileload['savelag'][i].value[0,:]
        
        
#x=savelag1['x']
#y=savelag1['y']
#z=savelag1['z']

#myz=z*0
#myz[:10000]=-3

mat=sio.loadmat('/home/moe46/workspace_matlab/lagtracker/element_starts/vh_high_test_large.mat')

locations=np.vstack([mat['x'],mat['y'],mat['z']]).T

print(locations.shape)

options={}
options['starttime']=0
options['endtime']=5+options['starttime']
options['interpolationratio']=15
options['outputratio']=1
options['ncformat']='NETCDF3_64BIT'
options['zlib']=False
options['useLL']=False
options['layer']=0
options['gridDim']='2D'


start = time.clock()
mypy=pyticle('/media/moe46/runs/vh_high/2012-02-01_2012-02-05/output/vh_high_0001.nc',locations,'test.nc',options=options)



mypy.run()
print('run in: %f' % (time.clock() - start))


#data=loadnc('',singlename='test.nc')


#s=100
#e=150
#plt.plot(data['x'][:300,s:e],data['y'][:300,s:e],'r')   
#plt.plot(data['x'][:300,s:e],data['y'][:300,s:e],'r.')   
#plt.plot(data['x'][0,s:e],data['y'][0,s:e],'k*')
#plt.plot(x[:,s:e],y[:,s:e],'b')
#plt.plot(x[:,s:e],y[:,s:e],'b.')
#plt.plot(x[0,s:e],y[0,s:e],'k*')
#plt.show()
 
        
        

