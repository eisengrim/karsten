# stdBathy
# Written by: Joel Culina

# This code calculates bathymetric standard deviation field sigma(x,y)
# Varying bottom roughness is then found as:
# z0(x,y) = alpha*sigma(x,y), where alpha is tuning (scalar) coefficient
# Based on Anderson and Meneveau (2011), doi:10.1017/jfm.2011.137
# Paper located in /EcoII/force/force_acadia_project/references/

# There are some parameters and functions that were chosen
# somewhat arbitrarily in calculating sigma(x,y):
# -Std of Gaussian filter/kernel is defined below as
# as the square root of element area (at each element)
# -Gaussian filter is used
# -


def triangleArea(x, y, trinodes):
    '''
    Calculates the area of the triangles in the grid.
    '''

    x1 = x[trinodes[:, 0] - 1]
    x2 = x[trinodes[:, 1] - 1]
    x3 = x[trinodes[:, 2] - 1]
    y1 = y[trinodes[:, 0] - 1]
    y2 = y[trinodes[:, 1] - 1]
    y3 = y[trinodes[:, 2] - 1]

    a = np.sqrt((x1 - x2)**2. + (y1 - y2)**2.)
    b = np.sqrt((x2 - x3)**2. + (y2 - y3)**2.)
    c = np.sqrt((x3 - x1)**2. + (y3 - y1)**2.)
    s = (a + b + c) / 2.

    tri_area = np.sqrt(s * (s - a) * (s - b) * (s - c))

    return tri_area


# Import modules

import numpy as np
import netCDF4 as net
import pyproj
import matplotlib.pyplot as plt
import cPickle as pickle

####################
#### User Input ####
####################

#### Multibeam Depth File ####
path_MB='/array/home/119865c/karsten/bathym/'
file_MB='GP_2m_MSL_Data.csv'

#### Model Files ####
path_model='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/sample_grid/'
file_model_nc='dngridCSR_sample.nc'

#### Output File ####
path_output='/array/home/119865c/karsten/std_bathy/'
output_file='GP_stdbathy_triArea_dngridCSR.p'

#### Bounding Box ####
# Use DG bounding box to define region
# westb = [-65.84, 44.64]
# northb = [-65.73, 44.72]

# Use PP bounding box to define region
# westb = [-66.23, 44.37]
# northb = [-66.19, 44.41]

# Use GP bounding box to define region
westb = [-66.38, 44.21]
northb = [-66.29, 44.32]

#### Projections ####
# dngridCSR projection
proj_lcc = pyproj.Proj(proj=
	'lcc +lon_0=-64.55880 +lat_0=41.78504 +lat_1=39.69152 +lat_2=43.87856')
proj_utm20 = pyproj.Proj(proj=
	'utm +zone=20T, +north +ellps=WGS84 +datum=WGS84 +unites=m +no_defs')

# acadia_BoF projection
# proj_lcc = pyproj.Proj(proj=
#     'lcc +lon_0=-64.55880 +lat_0=41.84493 +lat_1=39.72147 +lat_2=43.96838')
# proj_utm20 = pyproj.Proj(proj=
#     'utm +zone=20T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

###################
#### Load Data ####
###################

# Initialize
print 'initializing...'

xobs = []
yobs = []
latobs = []
lonobs = []
hobs = []
xmod = []
ymod = []
hmod = []

# MB dep
print 'opening file...\ngathering data...'
f = open(path_MB + file_MB, 'r')
for i, line in enumerate(f):
    if i == 0:
        continue
    line = line.strip()
    columns = line.split(',')
    lonobs.append(float(columns[0]))
    latobs.append(float(columns[1]))
    hobs.append(float(columns[2]))
f.close()
latobs = np.asarray(latobs)
lonobs = np.asarray(lonobs)
hobs = np.asarray(hobs)
h2obs = hobs ** 2.

# Convert lon/lat into x/y values
print 'converting lon/lat to x/y...'
xobs, yobs = proj_utm20(lonobs, latobs)

# Model grid and depth
print 'extracting netcdf data...'
nc = net.Dataset(path_model + file_model_nc)
trinodes = nc.variables['nv'][:, :].T
art1 = nc.variables['art1'][:]
hmod = nc.variables['h'][:]
lonc = nc.variables['lonc'][:]
latc = nc.variables['latc'][:]
x_vert = nc.variables['x'][:]
y_vert = nc.variables['y'][:]
nc.close()

# Find values at element centroids:
print 'finding values at element centroids...'
hmod = (hmod[trinodes[:, 0] - 1] + hmod[trinodes[:, 1] - 1]
        + hmod[trinodes[:, 2] - 1]) / 3.
art1 = (art1[trinodes[:, 0] - 1] + art1[trinodes[:, 1] - 1]
        + art1[trinodes[:, 2] - 1]) / 3.

# Change model proj to utm20 to match obs and calculate triangle area
print 'calculating triangle area...'
xmod, ymod = proj_utm20(lonc, latc)
tri_area = triangleArea(x_vert, y_vert, trinodes)

# Find find subset of elements (ROImod) over which find std
# Use rectangle bounding box instead of ellipse
print 'finding elements in bounding box...'
ROImod = np.where((lonc >= westb[0]) & (lonc <= northb[0]) &
                  (latc >= westb[1]) & (latc <= northb[1]))[0].tolist()

# Find local bathy standard deviation from obs at each model element:
print 'calculating sig...'
# Std of Gaussian filter/kernel, use triangle area
# sig = art1 ** (1/2.)
sig = np.sqrt(tri_area / np.pi)

# Initialize for filter loop
hobs_filt = np.zeros(hmod.shape[0])
std_bathy = np.zeros(hmod.shape[0])

# Filter at each point
print 'filtering at each point...'
for i, roi in enumerate(ROImod):
    print i, np.size(ROImod)

    diffx = xobs - xmod[roi]
    diffy = yobs - ymod[roi]
    dist = (diffx**2 + diffy**2)**(1/2.)
    near_node = np.where(dist <= (3. * sig[roi]))
    if near_node[0].size == 0:
        std_bathy[roi] = np.nan
        continue
    print 'Near node points: {}'.format(near_node)

    # Define pdf to be 0 beyond 3 std
    gauss = (np.sqrt(2*np.pi) * sig[roi]) ** (-2.) * \
            np.exp(-((xobs[near_node] - xmod[roi]) ** 2.
                   + (yobs[near_node] - ymod[roi]) ** 2.)
                   / (2 * (sig[roi]) ** 2.))
    gauss = np.asarray(gauss)
    gauss = gauss / np.sum(gauss)

    hobs_filt[roi] = np.sum(hobs[near_node]*gauss)

    h2obs_filt = np.sum(h2obs[near_node]*gauss)
    hobs_filt_2 = (np.sum(hobs[near_node]*gauss)) ** 2.
    std_bathy[roi] = (h2obs_filt - hobs_filt_2) ** (1/2.)

#### Save ####
print 'saving data to pickle...'
pickle.dump(std_bathy, open(path_output + output_file, 'wb'))
