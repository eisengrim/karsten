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

model=FVCOM("/array/home/119865c/workspace/acadia_bof_v2_3d_20170106.nc", debug=False)
drift1 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_1.mat", debug=False)
drift2 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_2.mat", debug=False)
drift3 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_3.mat", debug=False)
drift4 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_4.mat", debug=False)
drift5 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_5.mat", debug=False)
drift6 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_6.mat", debug=False)
drift7 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_7.mat", debug=False)
drift8 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_8.mat", debug=False)
drift9 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_9.mat", debug=False)
drift10 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_10.mat", debug=False)
drift11 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_11.mat", debug=False)
drift12 = Drifter("/array/home/119865c/workspace/mp_bmp/mp_bmp_20170106_12.mat", debug=False)
