import numpy as np
import os
import cPickle as pkl
from pyseidon_dvt import FVCOM
from scipy.interpolate import interp1d
from matplotlib import cm
import matplotlib.pyplot as plt
import netCDF4 as nc
import gc

limits = [[7.0,21.0], [8.0, 20.0], [9.0,19.0], [7.0, 17.0]]
clearance = 6.0
header = "PP_dngrid"
pklfile = "nautricity_curves.p"
powcurves = ["p14m", "p12m", "p10m", "p10m"]
path2files = "/home/thomas/Desktop/Acadia/Nautricity/"

# import data
data = nc.Dataset(path2files+'data_'+header+'.nc')
speed = data.variables['speed'][:]
depth = data.variables['depth'][:]
data.close()

# reload
fvcom = FVCOM(path2files+header+'CSR_2014-07-08_2014-07-15.nc')

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
feat = pkl.load(open(path2files+pklfile,"rb"))
for powcurv, bb in zip(powcurves, limits):
    downlimit = bb[0]
    uplimit = bb[1]

    pc = interp1d(feat["speed"], feat[powcurv], bounds_error=False, fill_value=0.0)
    power = pc(np.asarray(speed))
    power = np.ma.masked_where(depth<downlimit, power)
    power = np.ma.masked_where(depth>uplimit, power)
    speedma = np.ma.masked_where(depth<downlimit, speed)
    speedma = np.ma.masked_where(depth>uplimit, speed)
    meanspeed = np.squeeze(speed.mean(axis=1))
    meanpower = np.squeeze(np.apply_over_axes(np.nanmean, power, (0, 1)))
    # filtering location where turbine stick out of the water
    mindepth = np.nanmin(np.squeeze(depth[:,0,:]), axis=0)
    meanpower[np.where(mindepth < uplimit + clearance)] = np.nan

    fvcom.Plots.colormap_var(0.5*1.025*(meanspeed.mean(axis=0))**3.0,
                             title='Power density (kW/m2) between '+
                             str(int(downlimit))+'-'+str(int(uplimit))+' meters'+'-'+header,
                             #png=True, kmz=True)
                             shapefile=True, png=True, kmz=True)
    fvcom.Plots.colormap_var(0.5*1.025*(meanspeed.mean(axis=0))**3.0,
                             title='Power density (kW/m2) between '+
                              str(int(downlimit))+'-'+str(int(uplimit))+
                             ' meters - distance in meters',
                             #png=True, degree=False, kmz=True)
                             png=True, degree=False)

    #ax.plot(feat["speed"], feat[powcurv], label=powcurv)
    fvcom.Plots.colormap_var(meanpower,
                             title="Mean power (kW) between "+
                             str(int(downlimit)) + '-' + str(int(uplimit)) +
                             " - "+header+" - "+powcurv,
                             shapefile=True, png=True, kmz=True,
                             #png=True, kmz=True)#,
                             cmin=150.0, cmax=250.0, cmap=cm.jet)
    #plt.savefig("mean_power_"+str(int(downlimit))+'-'+str(int(uplimit))+
    #            "_"+header+"_"+powcurv+".png", bbox_inches='tight')
    fvcom.Plots.colormap_var(meanpower,
                             title="Mean power (kW) between "+
                             str(int(downlimit))+'-'+str(int(uplimit))+
                             " - "+header+" - "+powcurv+' - distance in meters',
                             #shapefile=True, png=True, kmz=True,
                             png=True,
                             cmin=150.0, cmax=250.0, cmap=cm.jet, degree=False)

    del meanpower, meanspeed, power, speedma
    gc.collect()

#ax.legend()

#min = np.min(fvcom.Grid.depth, axis=1)
#meanmin = min.mean(axis=0)
#for powcurv in powcurves:
#    pc = interp1d(feat["speed"], feat[powcurv])
#    power = pc(meanspeed)
#    meanpower = np.squeeze(np.mean(power, axis=0))
#    meanpower[np.where(meanmin>-25.0)] = np.nan
#    fvcom.Plots.colormap_var(meanpower,
#                             title="Mean power (kW) - "+header+" - "+powcurv+" and depth filter",
#                             shapefile=True,
#                             cmin=20.0, cmax=120.0, cmap=cm.jet)
#    plt.savefig("mean_power_"+header+"_"+powcurv+"_depthfiltered"+".png", bbox_inches='tight')
#    fvcom.Plots.colormap_var(meanpower,
#                             title="Mean power (kW) - "+header+" - "+powcurv+" and depth filter",
#                             cmin=20.0, cmax=120.0, cmap=cm.jet, degree=False)
#

#plt.show()
raw_input("Press Enter to exit...")
