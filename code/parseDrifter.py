import scipy.io as sio
import pandas as pd
import numpy as np

from utils import dn2dt

if "__name__" == "__main__":

    filename = "/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_all.mat"
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)

    dates = ["20130801", "20130802", "20130808", "20120605", "20120705"]
    dtt = np.asarray([d.seconds for d in np.diff(dt)])
    breaks = np.squeeze(np.where(dtt > 45))
    num = len(breaks) + 1
    breaks = np.insert(breaks, 0, -1)

    for k in dates:
        idx = [k in x for x in data['fn']]
        bks = breaks[np.where(idx)]
        u = []
        v = []
        angle = []
        speed = []
        beta = []
        alpha = []
        time = []
        lon = []
        lat = []
        for i in np.squeeze(np.where(idx)):
            u = np.append(u, data['velocity'].u[breaks[i]+1:breaks[i+1]])
            v = np.append(v, data['velocity'].v[breaks[i]+1:breaks[i+1]])
            angle = np.append(angle, data['velocity'].direction_deg[breaks[i]+1:breaks[i+1]])
            speed = np.append(speed, data['velocity'].speed[breaks[i]+1:breaks[i+1]])
            time = np.append(time, data['velocity'].vel_time[breaks[i]+1:breaks[i+1]])
            lon = np.append(lon, data['velocity'].vel_lon[breaks[i]+1:breaks[i+1]])
            lat = np.append(lat, data['velocity'].vel_lat[breaks[i]+1:breaks[i+1]])
            beta = np.append(beta, data['water_level'].beta[breaks[i]+1:breaks[i+1]])

    for k, fn in enumerate(data['fn'], start=0):
        pass

    drifters = {}
    drifters['u'] = data['velocity'].u
    drifters['v'] = data['velocity'].v
    drifters['speed'] = data['velocity'].speed
    drifters['time'] = data['velocity'].vel_time
    drifters['alpha'] = data['velocity'].alpha
    drifters['beta'] = data['water_level'].beta
    drifters['tide'] = data['water_level'].tide
    drifters['orientation'] = data['velocity'].direction_deg
    drifters['lat'] = data['velocity'].vel_lat
    drifters['lon'] = data['velocity'].vel_lon
    drifters['range'] = data['water_level'].tide_range
    drifters['ID'] = data['fn']

    for id, d in drifters.iteritems():
        ids.append(id)
        frames.append(pd.DataFrame(d).loc[~pd.DataFrame(d).index.duplicated(keep='first')])
    pd.concat(frames, keys=ids)
    drift = pd.concat(frames, keys=ids)
    drift.to_csv("~/drift_gp_all.csv")

