import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from datetime import datetime, timedelta
import netCDF4 as nc

# local imports - might interfere if using a different pyseidon
from drifter_statistics import *
from plot_utils import *
from utils import *
from color_map import createColorMap

# load pyseidon
from pyseidon_dvt import *
from interpolation_utils import *

model = FVCOM("/array/home/119865c/workspace/acadia_bof_v2_3d_20170106.nc", debug=False)
modelb = FVCOM("/array/home/119865c/workspace/acadia_bof_v2_3d_20170106_18.nc", debug=False)
drift = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_1.mat", debug=False)

driftmp = Drifter("/array/home/119865c/workspace/MP_BMP_20170106.mat", debug=False)
driftgp = Drifter("/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_all.mat", debug=False)

ddir = "/EcoII/acadia_uni/workspace/observed/GP/Drifter/"
ddirmp = "/array/home/119865c/workspace/mp_bmp/"
model2 = FVCOM("/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/" + \
        "drifter_runs/BFRIC_0.015/GP/2013_Aug_01_3D/output/subdomain_GP1_0001.nc", debug=False)
model1 = FVCOM("/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/" + \
        "drifter_runs/BFRIC_0.015/GP/2013_Aug_02_3D/output/subdomain_GP1_0001.nc", debug=False)
model8 = FVCOM("/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/" + \
        "drifter_runs/BFRIC_0.015/GP/2013_Aug_08_3D/output/subdomain_GP1_0001.nc", debug=False)
