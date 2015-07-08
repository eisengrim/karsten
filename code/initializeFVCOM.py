from pyseidon import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from datetime import datetime, timedelta

# ----------Grand Passage:
# PATH_TO_SIM_FILE="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/GP/2013_Aug_01_3D/output/subdomain_GP1_0001.nc"
# PATH_TO_OBS_FILE="/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_F_20130801_78_2_001_SE15.mat"
# PATH_TO_ADCP='/EcoII/acadia_uni/workspace/observed/GP/ADCP/Flow_GP-130730-TA_avg15.mat'

# ----------Digby Gut:
PATH_TO_SIM_FILE="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/DG/2013_Oct_10_3D/output/subdomain_DG1_0001.nc"
PATH_TO_OBS_FILE="/EcoII/acadia_uni/workspace/observed/DG/Drifter/DG_F_20131010_78_1_002_NNE05.mat"

PATH_TO_SIM_FILE2='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/DG/2013_Oct_10_3D/output/subdomain_DG1_0001.nc'

ax =''
tx1=''
tx2=''

if __name__ == '__main__':
    """
    The program initializes the FVCOM, Drifter, and Validation classes for the
    given FVCOM grid models and drifter file(s).
    """

    print 'using {} as model path and {} as observed path...'.format(PATH_TO_SIM_FILE, PATH_TO_OBS_FILE)
    model = FVCOM(PATH_TO_SIM_FILE, debug=True)
    model2= FVCOM(PATH_TO_SIM_FILE2, debug=True)
    drift = Drifter(PATH_TO_OBS_FILE, debug=True)
    # adcp = ADCP(PATH_TO_ADCP, debug=True)
    valid = Validation(drift, model,debug=True)
    valid2 = Validation(drift, model2, debug=True)

    mTimes = valid.Variables.struct['mod_time']
    oU = valid.Variables.struct['obs_timeseries']['u']
    oV = valid.Variables.struct['obs_timeseries']['v']
    mU = valid.Variables.struct['mod_timeseries']['u']
    mV = valid.Variables.struct['mod_timeseries']['v']
    speedS = np.asarray(np.sqrt(mU**2 + mV**2))
    speedO = np.asarray(np.sqrt(oU**2 + oV**2))

