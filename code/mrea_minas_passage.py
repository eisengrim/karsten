#!/usr/bin/python2.7
# encoding: utf-8

from __future__ import division
import numpy as np
from math import atan2
from pyseidon_dvt import FVCOM
import matplotlib.pyplot as plt
import multiprocessing as mp
from progressbar import ProgressBar
import cPickle as pkl
import os
import gc
import xarray as xr

### HEADER
# Global variables #
#  Number of cpus
cpus = int(mp.cpu_count() * 0.8)  # using 70% ish of the cpu resource
#  Paths
#path2files = '/EcoII/acadia_uni/projects/OERA_MREA/work/simulations/2011_runs/fvcom_files/'
#path2files = '/EcoII/force/force_acadia_project/work/simulations/acadia_bof_v2/2d/2011_monthly_runs/fvcom_files/'
#path2files = '/EcoII/force/force_acadia_project/work/simulations/acadia_bof_v2/2d/BF_0.00275/2011_monthly_runs/fvcom_files/'
path2files = '/EcoII/force/force_acadia_project/work/simulations/acadia_bof_v2/2d/BF_0.00275/2011-10-01_2011-10-31/output/'
#filename = 'acadia_BoF_0001.nc'
savename = 'fvcom_data_2011'
#  Threshold for slack water (m/s)
slack = 0.1  # 0.5
#time range of 1 lunar month = 57 tidal cycles =  29 days, 11h, 56m, 24s
time_range=['2011-10-01 00:00:00','2011-10-29 11:56:24']
#area = [-64.886019,-64.605919,45.243679,45.35815]
#area = [-64.82,-64.7,45.25,45.3]
#area = [-64.6,-64.3,45.3,45.4]
area = [-65.0,-63.4,45.0,46]



#  Colormaps
#cmap_spec = plt.get_cmap('viridis')  # color map for error
cmap_spec = plt.get_cmap('Spectral')  # color map for error
cmap_std = plt.get_cmap('rainbow')  # color map for vorticity
cmap_jet = plt.get_cmap('jet')  # directional spread
### END HEADER

### PROCESS
# define principal axis function
def principal_axis(u, v):

    # create velocity matrix
    U = np.vstack((u,v)).T
    # eliminate NaN values
    U = U[~np.isnan(U[:, 0]), :]
    # convert matrix to deviate form
    rep = np.tile(np.mean(U, axis=0), [len(U), 1])
    U -= rep
    # compute covariance matrix
    R = np.dot(U.T, U) / (len(U) - 1)

    # calculate eigenvalues and eigenvectors for covariance matrix
    R[np.where(np.isnan(R))] = 0.0
    R[np.where(np.isinf(R))] = 0.0
    R[np.where(np.isneginf(R))] = 0.0
    lamb, V = np.linalg.eig(R)
    # sort eigenvalues in descending order so that major axis is given by first eigenvector
    # sort in descending order with indices
    ilamb = sorted(range(len(lamb)), key=lambda k: lamb[k], reverse=True)
    lamb = sorted(lamb, reverse=True)
    # reconstruct the eigenvalue matrix
    lamb = np.diag(lamb)
    #reorder the eigenvectors
    V = V[:, ilamb]

    # rotation angle of major axis in radians relative to cartesian coordiantes
    #ra = np.arctan2(V[0,1], V[1,1])
    #ra = atan2(V[0,1], V[1,1])
    # express principal axis in compass coordinates
    #PA = -ra * 180 / np.pi + 90
    #TR: still not sure here
    PA = np.rad2deg(np.arctan2(V[0,0], V[1,0]))
    # variance captured by principal
    varxp_PA = np.diag(lamb[0]) / np.trace(lamb)

    return PA, varxp_PA

# define main loop
def main_loop(sesteps, ua, va, elc, powden, results):
    # Initialiasize array
    length = sesteps[1] - sesteps[0]
    pa = np.zeros(length)
    paf = np.zeros(length)
    pae = np.zeros(length)
    stda = np.zeros(length)
    stdf = np.zeros(length)
    stde = np.zeros(length)
    dirasy = np.zeros(length)
    speedFmean = np.zeros(length)
    speedEmean = np.zeros(length)
    speedFmedian = np.zeros(length)
    speedEmedian = np.zeros(length)
    speedFpercentile = np.zeros(length)
    speedEpercentile = np.zeros(length)
    speedFmax = np.zeros(length)
    speedEmax = np.zeros(length)
    powdenFmean = np.zeros(length)
    powdenEmean = np.zeros(length)
    powdenFmax = np.zeros(length)
    powdenEmax = np.zeros(length)

    # for ii in range(length):
    # Progress bar
    pbar = ProgressBar()
    for ii in pbar(range(length)):
        U = ua[:, ii].copy()
        V = va[:, ii].copy()
        EL = elc[:, ii].copy()
        powd = powden[:, ii]
        # masking out slack water speeds
        speed = np.sqrt(U**2.0 + V**2.0)
        U[np.where(speed < slack)] = np.nan
        V[np.where(speed < slack)] = np.nan
        EL[np.where(speed < slack)] = np.nan
        powd[np.where(speed < slack)] = np.nan
        speed[np.where(speed < slack)] = np.nan

        # compute flow's principal axis and its variance
        pr_axis, pr_ax_var = principal_axis(U, V)
        #note we use u,V to get directions in compass coords.
        dir_all = np.rad2deg(np.arctan2(U, V))
#        ind = np.where(dir_all < 0.0)
#        dir_all[ind] = dir_all[ind] + 360.0
        # sign speed - eliminating wrap-around
        dir_PA = dir_all - pr_axis
        dir_PA[dir_PA < -180.0] += 360.0
        dir_PA[dir_PA > 180.0] -= 360.0
        #if 0.0 < pr_axis < 90.0:
        #    pr_axis = 90.0 - pr_axis
        #elif 90.0 < pr_axis <180.0:
        #    pr_axis = 360.0 - (pr_axis - 90.0)
        #elif pr_axis < 0.0:
        #    pr_axis = 90.0 + np.abs(pr_axis)
        #general direction of flood passed as input argument
        fI = np.where((dir_PA >= -90.0) & (dir_PA < 90.0))[0]
        eI = np.arange(dir_PA.shape[0])
        eI = np.delete(eI, fI[:])
        #check if tide is not rising on flood, swith fI and eI in this case
        diff_EL=EL;
        diff_EL[0:-1]=np.diff(EL);
        diff_EL[-1]=0
        if np.nanmean(diff_EL[fI])<0:
           temp=fI
           fI=eI
           eI=temp		
           if pr_axis>=0:
              pr_axis=pr_axis-180
           else:
              pr_axis=pr_axis+180
           dir_PA[dir_PA>=0] -=180
           dir_PA[dir_PA<0] +=180
        pa[ii]=pr_axis
        # principal axis and stuff
        #pa, pav = principal_axis(U[fI], V[fI])
        #paf[ii] = np.average(dir_all[fI], weights=speed[fI]) # pa
        pf = np.rad2deg(np.arctan2(np.nanmean(U[fI]), np.nanmean(V[fI]) ))
        # True north conversion
        #if 0.0 < pf < 90.0:
        #    pf = 90.0 - pf
        #elif 90.0 < pf <180.0:
        #    pf = 360.0 - (pf - 90.0)
        #elif pf < 0.0:
        #    pf = 90.0 + np.abs(pf)
        paf[ii] = pf
        dirF = dir_PA[fI]
        stdf[ii] = np.std(np.ma.masked_where(np.isnan(dirF), dirF))
        #pa, pav = principal_axis(U[eI], V[eI])
        #pae[ii] = np.average(dir_all[eI], weights=speed[eI]) # pa
        pe = np.rad2deg(np.arctan2(np.nanmean(U[eI]), np.nanmean(V[eI])))
        # True north conversion
        #if 0.0 < pe < 90.0:
        #    pe = 90.0 - pe
        #elif 90.0 < pe < 180.0:
        #    pe = 360.0 - (pe - 90.0)
        #elif pe < 0.0:
        #    pe = 90.0 + np.abs(pe)
        pae[ii] = pe
        #if pr_axis < 0:
        #    pae[ii] = pr_axis + 180.0
        #else:
        #    pae[ii] = pr_axis - 180.0
        #stde[i] = np.sqrt(np.abs(pav[0,0]))
        dirE = dir_PA[eI]
        dirE[dirE>0]-=180
        dirE[dirE<0]+=180
        stde[ii] = np.std(np.ma.masked_where(np.isnan(dirE), dirE))
        stda[ii] = np.std(np.ma.masked_where(np.isnan(dir_PA), dir_PA))
        dirasy[ii] = np.abs(np.abs(paf[ii]-pae[ii])-180.0)

        # ebb flood velo
        speedFmean[ii] = np.nanmean(speed[fI])
        speedEmean[ii] = np.nanmean(speed[eI])
        speedFmedian[ii] = np.nanmedian(speed[fI])
        speedEmedian[ii] = np.nanmedian(speed[eI])
        speedFpercentile[ii] = np.nanpercentile(speed[fI], 70.0)
        speedEpercentile[ii] = np.nanpercentile(speed[eI], 70.0)
        try:
            speedFmax[ii] = np.nanmax(speed[fI])
        except ValueError:
            #speedFmax[ii] = 0.0
            speedFmax[ii] = np.nan
        try:
            speedEmax[ii] = np.nanmax(speed[eI])
        except ValueError:
            #speedEmax[ii] = 0.0
            speedEmax[ii] = np.nan
        powdenFmean[ii] = np.nanmean(powd[fI])
        powdenEmean[ii] = np.nanmean(powd[eI])
        try:
            powdenFmax[ii] = np.nanmax(powd[fI])
        except ValueError:
            #powdenFmax[ii] = 0.0
            powdenFmax[ii] = np.nan
        try:
            powdenEmax[ii] = np.nanmax(powd[eI])
        except ValueError:
            #powdenEmax[ii] = 0.0
            powdenEmax[ii] = np.nan

    output = {'pa': pa, 'paf': paf, 'pae': pae, 'stda': stda, 'stdf': stdf, 'stde': stde, 'dirasy': dirasy,
              'speedFmean': speedFmean, 'speedEmean': speedEmean,
              'speedFmedian': speedFmedian, 'speedEmedian': speedEmedian,
              'speedFpercentile': speedFpercentile, 'speedEpercentile': speedEpercentile,
              'speedFmax': speedFmax, 'speedEmax': speedEmax,
              'powdenFmean': powdenFmean, 'powdenEmean': powdenEmean,
              'powdenFmax': powdenFmax, 'powdenEmax': powdenEmax}

    print "main loop key: "+str(sesteps)
    results[str(sesteps)] = output
    # output.put((sesteps, results))

if __name__ == "__main__":
    # listing
    list_files = os.listdir(path2files)
    trim_list = []
    for name in list_files:
        if ".nc" in name:
            if not "station" in name.lower():
                if not "restart" in name.lower():
                    trim_list.append(name)
    trim_list.sort()

    # Serial extraction
    if not len(trim_list) == 1:
        for ii, file in enumerate(trim_list):
            fvcom=FVCOM(path2files+file)
#            fvcom.Util2D.depth_averaged_power_density()
            #fvcom.Util2D.vorticity()
            if ii == 0:
                flowspeed = fvcom.Variables.hori_velo_norm[:]
                powden = fvcom.Variables.depth_av_power_density[:]
                ua = fvcom.Variables.ua[:]
                va = fvcom.Variables.va[:]
                el = fvcom.Variables.el[:]
                #vort = fvcom.Variables.depth_av_vorticity[:]
            else:
                flowspeed = np.concatenate((flowspeed, fvcom.Variables.hori_velo_norm[:]), axis=0)
                powden = np.concatenate((powden, fvcom.Variables.depth_av_power_density[:]),
                                         axis=0)
                ua = np.concatenate((ua, fvcom.Variables.ua[:]), axis=0)
                va = np.concatenate((va, fvcom.Variables.va[:]), axis=0)
                el = np.concatenate((el, fvcom.Variables.el[:]), axis=0)
                #vort = np.concatenate((vort, fvcom.Variables.depth_av_vorticity[:]), axis=0)
            del fvcom
            gc.collect()
    else:
        filename=trim_list[0]
        fvcom = FVCOM(path2files+filename,tx=time_range,ax=area)
#        fvcom = FVCOM(path2files+filename)
        print "...computing vorticity"
        fvcom.Util2D.vorticity(debug=True)
        vorticity = np.nanmean(np.abs(fvcom.Variables.depth_av_vorticity), axis=0)
    	f = open("mrea_vorticity.p", "wb")
        pkl.dump({'vorticity': vorticity,
              'fvcom file': '/media/thomas/Data/Acadia/R50_v02/upper_BoF/acadia_upper_bof_v2_2d_0001.nc'}, f)
        f.close()
        fvcom.Variables.depth_av_vorticity=[]
        print "...computing power density"
        fvcom.Util2D.depth_averaged_power_density()
        powden = fvcom.Variables.depth_av_power_density[:]
        flowspeed = fvcom.Variables.hori_velo_norm[:]
        fvcom.Variables.depth_av_power_density=[]
        fvcom.Variables.hori_velo_norm=[]
        print "...saving variables"
        ua = fvcom.Variables.ua[:]
        va = fvcom.Variables.va[:]
        el = fvcom.Variables.el[:]
        gc.collect()

    elc=el[:,fvcom.Grid.trinodes[:,:]].mean(axis=2)
    # chunking spatial domain
    nsteps = ua.shape[1]
    dstep = int(nsteps/(cpus-1.0))
    sesteps = []
    uas = []
    vas =[]
    elcs =[]
    powdens = []
    estep = 0
    sstep = dstep
    for ii in range(cpus):
        if sstep > nsteps:
            sstep = nsteps
        if estep >= nsteps:
            break
        sesteps.append([estep, sstep])
        uas.append(ua[:, estep:sstep])
        vas.append(va[:, estep:sstep])
        elcs.append(elc[:, estep:sstep])
        powdens.append(powden[:,estep:sstep])
        # print "Dictionnary key: " + str([estep, sstep])
        estep += dstep
        sstep += dstep

    print "start parallel processing..."
    # Standard parallel block
    mgr = mp.Manager()
    d = mgr.dict()
    processes = [mp.Process(target=main_loop,
                 args=(sestep, u, v, elc, powd, d))
                 for sestep, u, v, elc, powd in zip(sesteps, uas, vas, elcs, powdens)]

    for ii, p in enumerate(processes):
        print "starting processes nb. " + str(ii) + "..."
        p.start()
    # wait until all processes started
    waiting = 1
    while waiting:
        for p in processes:
            print p.pid
            if p.pid is None:
                continue
        waiting = 0
        print "...all processes have started..."

    # joining processes
    for ii, p in enumerate(processes):
        print "joining processes nb. " + str(ii) + "..."
        p.join()

    print "...rebuilding matrices..."
    # rebuild matrix
    #  initialize matrices
    pa = np.zeros(nsteps)
    paf = np.zeros(nsteps)
    pae = np.zeros(nsteps)
    stda = np.zeros(nsteps)
    stdf = np.zeros(nsteps)
    stde = np.zeros(nsteps)
    dirasy = np.zeros(nsteps)
    speedFmean = np.zeros(nsteps)
    speedEmean = np.zeros(nsteps)
    speedFmedian = np.zeros(nsteps)
    speedEmedian = np.zeros(nsteps)
    speedFpercentile = np.zeros(nsteps)
    speedEpercentile = np.zeros(nsteps)
    speedFmax = np.zeros(nsteps)
    speedEmax = np.zeros(nsteps)
    powdenFmean = np.zeros(nsteps)
    powdenEmean = np.zeros(nsteps)
    powdenFmax = np.zeros(nsteps)
    powdenEmax = np.zeros(nsteps)

    # store in dict
    output_list = {'pa': pa, 'paf': paf, 'pae': pae, 'stda': stda, 'stdf': stdf, 'stde': stde, 'dirasy': dirasy,
                   'speedFmean': speedFmean, 'speedEmean': speedEmean,
                   'speedFmedian': speedFmedian, 'speedEmedian': speedEmedian,
                   'speedFpercentile': speedFpercentile, 'speedEpercentile': speedEpercentile,
                   'speedFmax': speedFmax, 'speedEmax': speedEmax,
                   'powdenFmean': powdenFmean, 'powdenEmean': powdenEmean,
                   'powdenFmax': powdenFmax, 'powdenEmax': powdenEmax}
    for ind in sesteps:
        for key in output_list.keys():
            if key in ['powdenFmean','powdenEmean',  'powdenFmax','powdenEmax']:
                output_list[key][ind[0]:ind[1]] = d[str(ind)][key][:]/1000.0
            else:
                output_list[key][ind[0]:ind[1]] = d[str(ind)][key][:]
    del mgr, d, processes
    print "...done parallel processing"

    # Dump in pickle file
    try:
        f = open(savename + ".p", "wb")
        pkl.dump(output_list, f)
        f.close()
    except:
        pass

#    # Load Dataset
#    fvcomData = xr.open_dataset(path2files+filename, decode_times=False, drop_variables=['siglev', 'siglay'])
#    # Vorticity
#    #t = np.arange(fvcomData.time.shape[0])
#    t = np.arange(fvcomData.time.shape[0]/2.)
#    print t
#    # Surrounding elements
#    n1 = fvcom.Grid.triele[:, 0]
#    n2 = fvcom.Grid.triele[:, 1]
#    n3 = fvcom.Grid.triele[:, 2]
#    ##change end bound indices
#    test = -1
#    n1[np.where(n1 == test)[0]] = 0
#    n2[np.where(n2 == test)[0]] = 0
#    n3[np.where(n3 == test)[0]] = 0
#    test = fvcom.Grid.nele - 1
#    n1[np.where(n1 > test)[0]] = -1
#    n2[np.where(n2 > test)[0]] = -1
#    n3[np.where(n3 > test)[0]] = -1
#    N1 = n1[:].astype(int)
#    N2 = n2[:].astype(int)
#    N3 = n3[:].astype(int)
#    #dvdx = xr.DataArray(np.zeros(fvcomData.ua.shape)) #, coords=flowSpeed.coords, dims=flowSpeed.dims)
#    #dudy = xr.DataArray(np.zeros(fvcomData.ua.shape)) #, coords=flowSpeed.coords, dims=flowSpeed.dims)
#    vort = xr.DataArray(np.zeros(fvcomData.ua.shape)) #, coords=flowSpeed.coords, dims=flowSpeed.dims)
#    j = 0
#    pbar = ProgressBar()
#    print "...computing vorticity"
#    for i in pbar(t):
#        vort[j, :] = (fvcomData.Grid.a1u[0, :] * fvcomData.va.values[i, :] \
#                     + fvcomData.Grid.a1u[1, :] * fvcomData.va.values[i, N1].T \
#                     + fvcomData.Grid.a1u[2, :] * fvcomData.va.values[i, N2].T \
#                     + fvcomData.Grid.a1u[3, :] * fvcomData.va.values[i, N3].T) \
#                   - (fvcomData.Grid.a2u[0, :] * fvcomData.ua.values[i, :] \
#                     + fvcomData.Grid.a2u[1, :] * fvcomData.ua.values[i, N1].T \
#                     + fvcomData.Grid.a2u[2, :] * fvcomData.ua.values[i, N2].T \
#                     + fvcomData.Grid.a2u[3, :] * fvcomData.ua.values[i, N3].T)
#        j += 1
#
    # Plot vorticity

    os.chdir('/EcoII/acadia_uni/projects/OERA_MREA/scripts/mrea_output')
    fvcom.Plots.colormap_var(vorticity, title='Horizontal vorticity', cmin=0.0,  # cmax=3.7, # isoline='var',
                             units='1/s', mesh=False, cmap=cmap_std, kmz=True, png=True, shapefile=True)
    
#    del vorticity, fvcomData
#    gc.collect()
### END PROCESS

### PLOT
    # Bathymetry
#    fvcom = FVCOM(path2files + filename)  # reload
    fvcom.Plots.colormap_var(fvcom.Grid.h, title='Bathymetry',isostep=10,
                             units='m', mesh=False, kmz=True, png=True, shapefile=True)
    # Slope
#    fvcom.Util2D.slope()
#    fvcom.Plots.colormap_var(fvcom.Grid.slope, title='Bathymetric slope',
#                             units='degrees', mesh=False, kmz=True, png=True, shapefile=True)
    # LLW
    llw = el.min(axis=0)
    fvcom.Plots.colormap_var(fvcom.Grid.h + llw, title='Low Low Water', isoline='var',
                             units='m', mesh=False, kmz=True, png=True, shapefile=True)
    del llw
    gc.collect()
    
    # Plot u mean
#    fvcom.Plots.colormap_var(np.nanmean(ua, axis=0), title='Mean u velocity component', #cmin=0.0, cmax=3.7,  # isoline='var',
#                             units='m/s', mesh=False, cmap=cmap_jet, kmz=True, png=True, shapefile=True)
    del ua
    gc.collect()
    # Plot v mean
 #   fvcom.Plots.colormap_var(np.nanmean(va, axis=0), title='Mean v velocity component',
                             # cmin=0.0, cmax=3.7,  # isoline='var',
 #                            units='m/s', mesh=False, cmap=cmap_jet, kmz=True, png=True, shapefile=True)
    del va
    gc.collect()
    # Plot mean speed
    fvcom.Plots.colormap_var(np.nanmean(flowspeed, axis=0), title='Mean speed', cmin = 0.0, cmax=2.0, # isoline='var',
                             units='m/s', mesh=False, cmap=cmap_jet,  kmz=True, png=True, shapefile=True)
    plt.close('all')
    # Plot cubic root of mean cubic speed
    crmcs = np.nanmean(flowspeed**3.0, axis=0)**(1.0/3.0)
    fvcom.Plots.colormap_var(crmcs, title='cubic root of mean cubic speed', cmin = 0.0,  # cmax=3.7, # isoline='var',
                             units='m/s', mesh=False, cmap=cmap_jet,  kmz=True, png=True, shapefile=True)
    f = open("mrea_cubic_root_of_mean_cubic_speed.p", "wb")
    pkl.dump({'cubic root of mean cubic speed': np.nanmean(flowspeed**3.0, axis=0)**(1.0/3.0),
              'fvcom file': '/media/thomas/Data/Acadia/R50_v02/upper_BoF/acadia_upper_bof_v2_2d_0001.nc'}, f)
    f.close()
    del crmcs
    gc.collect()
    # Plot 95th percentile of speed
#    speed95percentile = np.nanpercentile(flowspeed, 95.0, axis=0)
#    fvcom.Plots.colormap_var(speed95percentile, title='95th percentile flow speed', cmin = 0.0,  # cmax=3.7, # isoline='var',
#                             units='m/s', mesh=False, cmap=cmap_jet,  kmz=True, png=True, shapefile=True)
#    f = open("mrea_95th_percentile_flow_speed.p", "wb")
#    pkl.dump({'95th percentile flow speed': speed95percentile,
#              'fvcom file': '/media/thomas/Data/Acadia/R50_v02/upper_BoF/acadia_upper_bof_v2_2d_0001.nc'}, f)
#     f.close()
#     del speed95percentile
#     gc.collect()
    # Plot max speed
#     fvcom.Plots.colormap_var(flowspeed.max(axis=0), title='Max speed', cmin = 0.0, cmax=8.0,
#                              units='m/s', mesh=False, kmz=True, cmap=cmap_jet, png=True, shapefile=True)
#     del flowspeed
#     gc.collect()
    plt.close('all')
    # Plot mean power density
    fvcom.Plots.colormap_var(np.nanmean(powden/ 1000.0 , axis=0), title='Mean power density', cmin = 0.0, cmax=20.0, # conversion to kW.m-2
                             units='kW/m2', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    # Plot max power density
#    fvcom.Plots.colormap_var(np.nanmax(powden/ 1000.0 , axis=0), title='Max power density', #cmin = 3.0, cmax=9.0, # conversion to kW.m-2
#                             units='kW/m2', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    # Plot mean speed flood and ebb
    fvcom.Plots.colormap_var(output_list['speedFmean'][:], title='Mean speed flood', cmin = 0,  cmax=2.0,
                             units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
#    fvcom.Plots.colormap_var(output_list['speedFmax'][:], title='Max speed flood', cmin = 0.0, cmax=5.0,
#                             units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    fvcom.Plots.colormap_var(output_list['speedEmean'][:], title='Mean speed ebb', cmin = 0, cmax=2.0,
                             units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
#    fvcom.Plots.colormap_var(output_list['speedEmax'][:], title='Max speed ebb', cmin = 0.0, cmax=6.4,
#                             units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    plt.close('all')
    # Plot median speed flood and ebb
#     fvcom.Plots.colormap_var(output_list['speedFmedian'][:], title='Median speed flood', cmin = 0,  cmax=2.5,
#                              units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
#     fvcom.Plots.colormap_var(output_list['speedEmedian'][:], title='Median speed ebb', cmin = 0, cmax=4.2,
#                              units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    # Plot median speed flood and ebb
#     fvcom.Plots.colormap_var(output_list['speedFpercentile'][:], title='70th percentile speed flood', cmin = 0.0,  cmax=5.0,
#                              units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
#     fvcom.Plots.colormap_var(output_list['speedEpercentile'][:], title='70th percentile speed ebb', cmin = 0.0,  cmax=5.7,
#                              units='m/s', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    # Plot mean power density flood and ebb
    fvcom.Plots.colormap_var(output_list['powdenFmean'][:],
                              title='Mean power density flood', cmin = 0.0, cmax=10.0,
                              units='kW/m2', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    fvcom.Plots.colormap_var(output_list['powdenFmax'][:],
                              title='Max power density flood', cmin = 0.0, cmax=10.0,
                              units='kW/m2', mesh=False, kmz=True,  cmap=cmap_jet, png=True, shapefile=True)
    fvcom.Plots.colormap_var(output_list['powdenEmean'][:],
                              title='Mean power density ebb', cmin = 0.0, cmax=10.0,
                              units='kW/m2', mesh=False, kmz=True, cmap=cmap_jet, png=True, shapefile=True)
    fvcom.Plots.colormap_var(output_list['powdenEmax'][:],
                             title='Max power density ebb', cmin = 0.0,  cmax=10.0,
                              units='kW/m2', mesh=False, kmz=True, cmap=cmap_jet, png=True, shapefile=True)
#     plt.close('all')
    # Plot power asymetry
    powasy = output_list['powdenEmean'][:]/output_list['powdenFmean'][:]
    powasy[np.where(np.isnan(powasy))] = 50.0
    fvcom.Plots.colormap_var(powasy, title='Mean power density asymmetry ebb vs flood', cmin = 0.0, cmax=10.0,
                             mesh=False, kmz=True, cmap=cmap_jet, png=True, shapefile=True)
    # Plot directional plots
    fvcom.Plots.colormap_var(output_list['pa'][:], title='Principal axis', cmin = -180.0, cmax= 180.0,
                              units='deg.', mesh=False, kmz=True, cmap=cmap_spec, png=True, shapefile=True)
    fvcom.Plots.colormap_var(output_list['paf'][:], title='Principal axis flood', cmin = -180.0, cmax= 180.0,
                              units='deg.', mesh=False, kmz=True, cmap=cmap_spec, png=True, shapefile=True)
    # TR: weird things happening here with ebb direction
    fvcom.Plots.colormap_var(output_list['pae'][:], title='Principal axis ebb', cmin=-180.0,  cmax=180.0,
                              units='deg.', mesh=False, kmz=True, cmap=cmap_spec, png=True, shapefile=True)
#    fvcom.Plots.colormap_var(output_list['paf'][:] + output_list['dirasy'][:], title='Principal axis ebb', cmin = 130.0, # cmax=180.0,
#                             units='deg.', mesh=False, kmz=True, cmap=cmap_spec, png=True, shapefile=True)
    stdf = np.ma.masked_where(np.isnan(stdf), stdf)
    output_list['stdf'][np.isnan(output_list['stdf'][:])] = 90.0
    fvcom.Plots.colormap_var(output_list['stdf'][:], title='Standard deviation flood', cmin=0.0, #cmax=65.0,
                              mesh=False, kmz=True, cmap=cmap_std, png=True, shapefile=True)
    stde = np.ma.masked_where(np.isnan(stde), stde)
    output_list['stde'][np.isnan(output_list['stde'][:])] = 90.0
    fvcom.Plots.colormap_var(output_list['stde'][:], title='Standard deviation ebb', cmin=0.0, #cmax=65.0,
                              mesh=False, kmz=True, cmap=cmap_std, png=True, shapefile=True)
    output_list['stda'][np.isnan(output_list['stda'][:])] = 90.0
    fvcom.Plots.colormap_var(output_list['stda'][:], title='Standard deviation', cmin=0.0, #cmax=65.0,
                              mesh=False, kmz=True, cmap=cmap_std, png=True, shapefile=True)
    fvcom.Plots.colormap_var(output_list['dirasy'][:], title='Directional asymetry', cmin = 0.0, cmax=180.0,
                             units='deg.', mesh=False, kmz=True, cmap=cmap_std, png=True, shapefile=True)

    
    raw_input("---Press enter to exit---")
