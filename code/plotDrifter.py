#!usr/bin/python2.7
# encoding: utf-8

import sys, os.path
from pyseidon import *
import numpy as np
import matplotlib.pyplot as plt
from interpolation_utils import *

PATH_TO_SIM_FILE="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/DG/2013_Nov_05_3D/output/subdomain_DG1_0001.nc"
PATH_TO_OBS_FILE="/EcoII/acadia_uni/workspace/observed/DG/Drifter/DG_E_20131105_78_1_001_S5.mat"
LOCATION = 'DG'

"""
The program plots a drifter trajectory taken from a drifter class onto a
spatial varying map of the flow speed (specifically, the flood or ebb
velocity norm). Assumes the FVCOM and Drifter objects already loaded.
The program takes a location tag and works around that domain.
"""


if __name__ == '__main__':

    if not os.path.isfile(PATH_TO_SIM_FILE):
	sys.exit('FVCOM file invalid  / not found')
    if not os.path.isfile(PATH_TO_OBS_FILE):
	sys.exit('Drifter file invalid / not found')

    model = FVCOM(PATH_TO_SIM_FILE, debug=False)
    drift = Drifter(PATH_TO_OBS_FILE, debug=False)
    location = LOCATION

    # location = raw_input('Please enter the location tag (GP, DG, PP, MP): ')
    if location == 'GP':
        loc = [-66.33906, 44.26898]
    elif location == 'DG':
	loc = [-65.76000, 44.67751]
    elif location == 'PP':
        loc = [-66.20692, 44.38937]
    elif location == 'MP':
	loc = [-64.40725, 45.34758]
    else: sys.exit("Not a valid location tag.")

    fI, eI, pa, pav = model.Util2D.ebb_flood_split_at_point(loc[0], loc[1])
    model.Util3D.velo_norm()
    tide = str(drift.Data['water_level'].tide)

    tModel = model.Variables.matlabTime
    tDrift = drift.Variables.matlabTime
    win1 = (np.abs(tModel-tDrift.min())).argmin()
    win2 = (np.abs(tModel-tDrift.max())).argmin()

    #if tide == 'flood':
    #	tideNorm = np.mean(model.Variables.velo_norm[fI,:,:],0)
    #elif tide == 'ebb':
    #	tideNorm = np.mean(model.Variables.velo_norm[eI,:,:],0)

    if tide == 'flood':
        tideNorm = np.mean(model.Variables.velo_norm[win1:win2,:,:],0)
    elif tide == 'ebb':
        tideNorm = np.mean(model.Variables.velo_norm[win1:win2,:,:],0)

    model.Plots.colormap_var(tideNorm[0,:], mesh=False)
    plt.hold('on')

    pts = closest_points(drift.Variables.lon[:], drift.Variables.lat[:], \
                model.Grid.lon, model.Grid.lat)
    xObs = model.Grid.x[pts]
    yObs = model.Grid.y[pts]

    x = drift.Variables.lon
    y = drift.Variables.lat
    u = drift.Variables.u
    v = drift.Variables.v

    plt.quiver(x,y,u,v)
    # plt.scatter(x,y)
    plt.show()



