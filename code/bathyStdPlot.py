import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import netCDF4 as net
import seaborn

# load in variables
std_bathy_file = '/EcoII/acadia_uni/projects/ecoEII/MP_std_bathy_triArea.p'
std_bathy = pickle.load(open(std_bathy_file, 'rb'))
nc_path = '/EcoII/acadia_uni/projects/acadia_bay_of_fundy_numerical_model_r50/work/simulations/acadia_bay_of_fundy_numerical_model_r50_v01/calibration/Cd/2011-08-12_2011-11-03_0.00225/output/'
nc = net.Dataset(nc_path + 'acadia_BoF_0001.nc')
lonc = nc.variables['lonc'][:]
latc = nc.variables['latc'][:]
x = nc.variables['x'][:]
y = nc.variables['y'][:]
trinodes = nc.variables['nv'][:, :].T - 1
nc.close()

'''
# define GP bounding box and find points within box
westb = [-66.36, 44.24]
northb = [-66.31, 44.3]
'''

# MP box
westb = [-64.40359722778196, 45.362024744162525]
northb = [-64.43292-.01, 45.37389]

'''
# define PP bounding box and find points within box
westb = [-66.23, 44.37]
northb = [-66.19, 44.41]
'''

'''
# define DG bounding box and find points within box
westb = [-65.84, 44.64]
northb = [-65.73, 44.72]
'''

GP_ind = np.where((lonc >= westb[0]) & (lonc <= northb[0]) &
                  (latc >= westb[1]) & (latc <= northb[1]))[0].tolist()

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
