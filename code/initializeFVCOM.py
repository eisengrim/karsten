from pyseidon import *
from interpolation_utils import *
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta

# ----------Grand Passage:
# PATH_TO_SIM_FILE="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/GP/2013_Aug_01_3D/output/subdomain_GP1_0001.nc"
# PATH_TO_OBS_FILE="/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_F_20130801_78_2_001_SE15.mat"
# PATH_TO_ADCP='/EcoII/acadia_uni/workspace/observed/GP/ADCP/Flow_GP-130730-TA_avg15.mat'

# ----------Digby Gut:

LOC = 'GP'
DATE = '2013_Aug_01_3D'
BF = ['0.015', '0.012', '0.009', '0.020']


PATH2SIM="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/" \
        + "drifter_runs/BFRIC_" + BF[0] + "/" + LOC + "/" + DATE + "/output/" \
        + "subdomain_" + LOC + "1_0001.nc"
PATH2SIM2='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/' \
        + 'drifter_runs/BFRIC_' + BF[1] + '/' + LOC + "/" + DATE + '/output/' \
        + 'subdomain_' + LOC + '1_0001.nc'
PATH2SIM3='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/' \
        + 'drifter_runs/BFRIC_' + BF[2] + '/' + LOC + "/" + DATE + '/output/' \
        + 'subdomain_' + LOC + '1_0001.nc'

PATH2OBS="/EcoII/acadia_uni/workspace/observed/" + LOC + "/Drifter/" \
        + "GP_F_20130801_78_1_001_SE15.mat"


if __name__ == '__main__':
    """
    The program initializes the FVCOM, Drifter, and Validation classes for the
    given FVCOM grid models and drifter file(s).
    """

    model = FVCOM(PATH2SIM, debug=True)
    model2= FVCOM(PATH2SIM2, debug=True)
    model3= FVCOM(PATH2SIM3, debug=True)

    drift = Drifter(PATH2OBS, debug=True)
    # adcp = ADCP(PATH_TO_ADCP, debug=True)

    # create validation objects
    valid = Validation(drift, model,debug=True)
    valid2 = Validation(drift, model2, debug=True)
    valid3 = Validation(drift, model3, debug=True)

