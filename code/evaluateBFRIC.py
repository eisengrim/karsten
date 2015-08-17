from pyseidon import *
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta


# ----------File Parameters:

LOC = 'DG'
DATE = '2013_Oct_10_3D'
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
        + 'DG_F_20131010_78_2_002_NNE05.mat'



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

    # The speed calculations we want!
    mTimes = valid.Variables.struct['mod_time']
    oU = valid.Variables.struct['obs_timeseries']['u']
    oV = valid.Variables.struct['obs_timeseries']['v']
    speedO = np.asarray(np.sqrt(oU**2 + oV**2))

    datetimes = np.asarray([datetime.fromordinal(int(time)) + \
                          timedelta(days=time % 1) - \
                          timedelta(days=366)for time in mTimes])

    mU = valid.Variables.struct['mod_timeseries']['u']
    mV = valid.Variables.struct['mod_timeseries']['v']
    speedS = np.asarray(np.sqrt(mU**2 + mV**2))

    mU2 = valid2.Variables.struct['mod_timeseries']['u']
    mV2 = valid2.Variables.struct['mod_timeseries']['v']
    speedS2 = np.asarray(np.sqrt(mU2**2 + mV2**2))

    mU3 = valid3.Variables.struct['mod_timeseries']['u']
    mV3 = valid3.Variables.struct['mod_timeseries']['v']
    speedS3 = np.asarray(np.sqrt(mU3**2 + mV3**2))

    # create figure
    ax = plt.figure()
    fig = ax.add_subplot(1,1,1,axisbg='#F5DEB3')
    fig.plot(datetimes, speedO, '#8B0000', label='Drifter', linewidth=2)
    fig.hold('on')
    fig.plot(datetimes, speedS, 'b--',label='Model (BFRIC=0.015)', linewidth=2)
    fig.hold('on')
    fig.plot(datetimes, speedS2, 'g:', label='Model (BFRIC=0.012)',linewidth=2)
    fig.hold('on')
    fig.plot(datetimes, speedS3, 'r-.', label='Model (BFRIC=0.009)', linewidth=2)
    fig.set_ylabel('Speed (m/s)')
    fig.set_xlabel('Time (HH:MM:SS)')
    fig.set_title('Observed and Simulated Speed vs. Time')
    fig.legend(loc='upper left')
    plt.gcf().autofmt_xdate()
    fig.grid(True)

    print 'using: \n' + DATE + '\n' + PATH2OBS
