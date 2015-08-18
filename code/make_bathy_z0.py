import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle

# parameters
kappa = 0.4
min_depth = 3

# passage bounding boxes
GP = [-66.36, -66.31, 44.24, 44.3]
PP = [-66.23, -66.19, 44.37, 44.41]
DG = [-65.84, -65.73, 44.64, 44.72]
MP = [-65.5, -63.3, 45.0, 46.0]

# for later, save calculated initial guesses for geographic z0
MP_z = 0.005
GP_z = 0.0098
PP_z = 0.0323
DG_z = 0.0101

# initial guesses for bottom roughness
GP_cd = 0.0025
PP_cd = 0.00275
DG_cd = 0.00175


def bathyBR(alpha_set, cd, grid, outfile):
    '''
    Create a variable bottom roughness data file given user inputs of the
    cd and H parameters, and the grid on which we're building the file.
    Creates an nc file readable by FVCOM.
    '''

    cd = float(cd)

    # read in the grid nc file variables
    grid_file = nc.Dataset(grid)
    lon = grid_file.variables['lonc'][:]
    lat = grid_file.variables['latc'][:]
    trinodes = grid_file.variables['nv']
    h = grid_file.variables['h']
    h = np.array(h)
    htri = (h[trinodes[0, :] - 1] + h[trinodes[1, :] - 1] + h[trinodes[2, :] - 1]) / 3

    # calculate indices for each passage and mean heights
    mp_region = np.argwhere((lon >= MP[0]) & (lon <= MP[1]) &
                            (lat >= MP[2]) & (lat <= MP[3]))
    mp_region = mp_region.flatten()
    # mp_H = np.mean(htri[dg_region])
    # print 'MP mean height: {}'.format(mp_H)

    dg_region = np.argwhere((lon >= DG[0]) & (lon <= DG[1]) &
                            (lat >= DG[2]) & (lat <= DG[3]))
    dg_region = dg_region.flatten()
    # dg_H = np.mean(htri[dg_region])
    # print 'DG mean height: {}'.format(dg_H)

    gp_region = np.argwhere((lon >= GP[0]) & (lon <= GP[1]) &
                            (lat >= GP[2]) & (lat <= GP[3]))
    gp_region = gp_region.flatten()
    # gp_H = np.mean(htri[gp_region])
    # print 'GP mean height: {}'.format(gp_H)

    pp_region = np.argwhere((lon >= PP[0]) & (lon <= PP[1]) &
                            (lat >= PP[2]) & (lat <= PP[3]))
    pp_region = pp_region.flatten()
    # pp_H = np.mean(htri[pp_region])
    # print 'PP mean height: {}'.format(pp_H)

    # load in bathymetry z0 files
    GP_file = 'std_bathy_files/GP_std_bathy_triArea.p'
    DG_file = 'std_bathy_files/DG_std_bathy_triArea.p'
    PP_file = 'std_bathy_files/PP_std_bathy_triArea.p'
    MP_file = 'std_bathy_files/MP_std_bathy_triArea_combined.p'
    GP_bathy = pickle.load(open(GP_file, 'rb'))
    DG_bathy = pickle.load(open(DG_file, 'rb'))
    PP_bathy = pickle.load(open(PP_file, 'rb'))
    MP_bathy = pickle.load(open(MP_file, 'rb'))    

    # set z0 for entire grid, begin with default
    # take region z0 values from bathymetry std, scale by passage alpha 
    z0 = htri * np.exp(-(kappa * cd**(-0.5) + 1))
    z0[mp_region] = MP_bathy[mp_region] * alpha_set[0]
    z0[gp_region] = GP_bathy[gp_region] * alpha_set[1]
    z0[pp_region] = PP_bathy[pp_region] * alpha_set[2]
    z0[dg_region] = DG_bathy[dg_region] * alpha_set[3]
    z0[np.where(htri <= min_depth)] = \
        min_depth * np.exp(-(kappa * cd**(-0.5) + 1))
    z0[np.where(np.isnan(z0))] = \
        htri * np.exp(-(kappa * cd**(-0.5) + 1))
    z0[z0 == 0] = np.min(z0[z0 > 0])

    print np.sort(z0)

    '''
    # calculate z0, begin with default across entire grid
    z0 = htri * np.exp(-(kappa * cd**(-0.5) + 1))
    z0[dg_region] = dg_H * np.exp(-(kappa * dg_cd**(-0.5) + 1))
    z0[gp_region] = gp_H * np.exp(-(kappa * gp_cd**(-0.5) + 1))
    z0[pp_region] = pp_H * np.exp(-(kappa * pp_cd**(-0.5) + 1))
    z0[np.where(htri <= min_depth)] = \
        min_depth * np.exp(-(kappa * cd**(-0.5) + 1))
    z0 = z0.astype(np.float64)
    '''

    # pass the data onto an nc file
    br2nc(z0, outfile)


def br2nc(z0, outfile):
    '''
    Converts a bottom roughness array to a netcdf3 file, readable
    by FVCOM.
    '''

    # create netcdf4 file
    n = z0.size
    nc_id = nc.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')

    # set up the nc file and put all the data into it
    nc_id.createDimension('nele', n)
    z0b = nc_id.createVariable('z0b', 'f8', ('nele'))
    z0b.name = 'bottom roughness lengthscale'
    z0b.units = 'metres'
    z0b[:] = z0

    nc_id.close()
