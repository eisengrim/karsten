import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import netCDF4 as net
import seaborn

# load in variables
path2pkl="/EcoII/acadia_uni/projects/drifters/std_bathy_3D_dngridCSR/"
path2nc="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/sample_grid/"
ncfile="dngridCSR_sample.nc"
std_bathy = pickle.load(open(path2pkl + 'PP_stdbathy_triArea_dngridCSR.p', 'rb'))

nc = net.Dataset(path2nc+ncfile)
lonc = nc.variables['lonc'][:]
latc = nc.variables['latc'][:]
x = nc.variables['x'][:]
y = nc.variables['y'][:]
trinodes = nc.variables['nv'][:, :].T - 1
nc.close()

print 'trinodes: ', trinodes.shape
print 'std_bathy: ', std_bathy.shape

'''
# define GP bounding box and find points within box
westb = [-66.36, 44.24]
northb = [-66.31, 44.3]
'''
'''
# MP box
westb = [-65.5, 45.0]
northb = [-63.3, 46.0]
'''

# define PP bounding box and find points within box
westb = [-66.23, 44.37]
northb = [-66.19, 44.41]

'''
# define DG bounding box and find points within box
westb = [-65.84, 44.64]
northb = [-65.73, 44.72]
'''

print 'generating plot...'
log_bathy = np.log(std_bathy)
log_bathy[np.where(log_bathy < -3)] = -3

log_bathy[np.isnan(log_bathy)] = -3

plt.figure()
plt.gca().set_aspect('equal')
colors = cm.jet

plt.tripcolor(x, y, trinodes, facecolors=log_bathy, cmap=colors)
plt.colorbar()
# plt.colorbar(sp)
plt.show()
