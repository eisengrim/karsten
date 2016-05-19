#!/usr/bin/env python2.7

import pandas as pd

track_cols=["ID", "Timestamp", "Pos_X", "Pos_Y", "Pos_Z", "Vel_X", "Vel_Y", "Vel_Z", "GVel_X", "GVel_Y", "GVel_Z", "Yaw", "Pitch", "Roll", "Force_X", "Force_Y", "Force_Z", "LocalVF_X", "LocalVF_Y", "LocalVF_Z", "Target_X", "Target_Y", "Target_Z", "O_Rand", "V_Rand", "Mass", "Last Element", "Status", "Exited", "Exit Status"]
def read_track_file(filename, nrows=None, particles=None, chunksize=5000):
    if particles:
        parts = pd.read_csv(filename, delimiter="\t", skipinitialspace=True, names=track_cols, skiprows=1, iterator=True, chunksize=chunksize, nrows=nrows)
        return pd.concat([part[part["ID"].isin(particles)] for part in parts])

    tracks=pd.read_csv(filename, delimiter="\t", skipinitialspace=True, names=track_cols, nrows=nrows, skiprows=1)
    return tracks
