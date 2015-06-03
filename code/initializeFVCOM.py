from pyseidon import *
import matplotlib.pyplot as plt
import numpy as np

# ----------Grand Passage:
PATH_TO_SIM_FILE="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/GP/2013_Aug_01_3D/output/subdomain_GP1_0001.nc"
PATH_TO_OBS_FILE="/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_F_20130801_78_2_001_SE15.mat"

# ----------Digby Gut:
#PATH_TO_SIM_FILE="/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/DG/2013_Nov_05_3D/output/subdomain_DG1_0001.nc"
#PATH_TO_OBS_FILE="/EcoII/acadia_uni/workspace/observed/DG/Drifter/DG_E_20131105_78_1_001_S5.mat"


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
    drift = Drifter(PATH_TO_OBS_FILE, debug=True)
    valid = Validation(drift,model,debug=True)
