import scipy.io as sio
import pandas as pd

if "__name__" == "__main__":

    # filename = "/EcoII/Luna/Drifter_data/MinasPassage/MP_BMP_20170106_processed.mat"
    filename = "/array/home/119865c/MP_BMP_20170106.mat"
    data = sio.loadmat(filename, struct_as_record=False, sqeeze_me=True)

    drifters = {}
    drifters['speed'] = data['speed']
    drifters['v'] = data['v']
    drifters['u'] = data['u']
    drifters['time'] = data['time']
    drifters['alpha'] = data['tide_time']
    drifters['beta'] = data['Tr']
    drifters['ID'] = data['ID']
    drifters['orientation'] = data['angle_deg_N0']
    drifters['lat'] = data['lat']
    drifters['lon'] = data['lon']
    drifters['distance'] = data['dist']
    drifters['range'] = data['range']

    for id, d in drifters.iteritems():
        ids.append(id)
        frames.append(pd.DataFrame(d).loc[~pd.DataFrame(d).index.duplicated(keep='first')])
    pd.concat(frames, keys=ids)
    drift = pd.concat(frames, keys=ids)
    drift.to_csv("~/drift.csv")

