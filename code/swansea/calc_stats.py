# coding: utf-8
import matplotlib
import numpy as np
import seaborn as sb
import pandas as pd

AA = []
AE = []
CC = []
EE = []

def do_calcs():
    data = pd.read_csv("gathered_thesis.txt")
    data0 = data[data["ts"]==0].copy()
    dataF = data.drop(data0.index).copy()

    # Have to do the reset *after* using data0 to prune dataF!
    data0 = data0.reset_index()
    dataF = dataF.reset_index()

    AA = dataF.Simulation=="A-A-HP5MM"
    AE = dataF.Simulation=="A-E-HP5MM"
    EE = dataF.Simulation=="E-E-HP5MM"
    CC = dataF.Simulation=="C-C-HP5MM"

    maxvals = dataF[dataF["Batch"]==2000]

    dataF.loc[AE, "mpx_divisor"] = maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "mean_pos_x"].values
    dataF.loc[CC, "mpx_divisor"] = maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "mean_pos_x"].values
    dataF.loc[EE, "mpx_divisor"] = maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "mean_pos_x"].values
    dataF.loc[AA, "mpx_divisor"] = maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "mean_pos_x"].values

    dataF["rel_mpx_err"] = abs(dataF["mean_pos_x"] - dataF["mpx_divisor"]) / dataF["mpx_divisor"] # Relative Error between mean_pos_x and 2000 batch case

    dataF.loc[AE, "mpy_divisor"] = maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "mean_pos_y"].values
    dataF.loc[CC, "mpy_divisor"] = maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "mean_pos_y"].values
    dataF.loc[EE, "mpy_divisor"] = maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "mean_pos_y"].values
    dataF.loc[AA, "mpy_divisor"] = maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "mean_pos_y"].values
    dataF["rel_mpy_err"] = abs(dataF["mean_pos_y"] - dataF["mpy_divisor"]) / dataF["mpy_divisor"] # Relative Error between mean_pos_y and 2000 batch case

    dataF["delta_mpx"] = abs(dataF.mean_pos_x - data0.mean_pos_x)
    dataF["delta_mpy"] = abs(dataF.mean_pos_y - data0.mean_pos_y)

    maxvals = dataF[dataF["Batch"]==2000]

    ae_b2k_x = maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "delta_mpx"].values
    cc_b2k_x = maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "delta_mpx"].values
    ee_b2k_x = maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "delta_mpx"].values
    aa_b2k_x = maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "delta_mpx"].values

    ae_b2k_y = maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "delta_mpy"].values
    cc_b2k_y = maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "delta_mpy"].values
    ee_b2k_y = maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "delta_mpy"].values
    aa_b2k_y = maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "delta_mpy"].values

    dataF.loc[AE, "error_mpx"] = abs(dataF.loc[AE, "delta_mpx"] - ae_b2k_x) / ae_b2k_x
    dataF.loc[CC, "error_mpx"] = abs(dataF.loc[CC, "delta_mpx"] - cc_b2k_x) / cc_b2k_x
    dataF.loc[EE, "error_mpx"] = abs(dataF.loc[EE, "delta_mpx"] - ee_b2k_x) / ee_b2k_x
    dataF.loc[AA, "error_mpx"] = abs(dataF.loc[AA, "delta_mpx"] - aa_b2k_x) / aa_b2k_x

    dataF.loc[AE, "error_mpy"] = abs(dataF.loc[AE, "delta_mpy"] - ae_b2k_y) / ae_b2k_y
    dataF.loc[CC, "error_mpy"] = abs(dataF.loc[CC, "delta_mpy"] - cc_b2k_y) / cc_b2k_y
    dataF.loc[EE, "error_mpy"] = abs(dataF.loc[EE, "delta_mpy"] - ee_b2k_y) / ee_b2k_y
    dataF.loc[AA, "error_mpy"] = abs(dataF.loc[AA, "delta_mpy"] - aa_b2k_y) / aa_b2k_y

    maxvals = dataF[dataF["Batch"]==2000]

    dataF.loc[AE, "error_stddev_x"] = abs(dataF.loc[AE, "stddev_pos_x"] - maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "stddev_pos_x"].values) / maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "stddev_pos_x"].values
    dataF.loc[CC, "error_stddev_x"] = abs(dataF.loc[CC, "stddev_pos_x"] - maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "stddev_pos_x"].values) / maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "stddev_pos_x"].values
    dataF.loc[EE, "error_stddev_x"] = abs(dataF.loc[EE, "stddev_pos_x"] - maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "stddev_pos_x"].values) / maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "stddev_pos_x"].values
    dataF.loc[AA, "error_stddev_x"] = abs(dataF.loc[AA, "stddev_pos_x"] - maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "stddev_pos_x"].values) / maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "stddev_pos_x"].values

    dataF.loc[AE, "error_stddev_y"] = abs(dataF.loc[AE, "stddev_pos_y"] - maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "stddev_pos_y"].values) / maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "stddev_pos_y"].values
    dataF.loc[CC, "error_stddev_y"] = abs(dataF.loc[CC, "stddev_pos_y"] - maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "stddev_pos_y"].values) / maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "stddev_pos_y"].values
    dataF.loc[EE, "error_stddev_y"] = abs(dataF.loc[EE, "stddev_pos_y"] - maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "stddev_pos_y"].values) / maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "stddev_pos_y"].values
    dataF.loc[AA, "error_stddev_y"] = abs(dataF.loc[AA, "stddev_pos_y"] - maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "stddev_pos_y"].values) / maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "stddev_pos_y"].values

    dataF["delta_stddev_x"] =  dataF["stddev_pos_x"] - data0["stddev_pos_x"]
    dataF["delta_stddev_y"] =  dataF["stddev_pos_y"] - data0["stddev_pos_y"]

    maxvals = dataF[dataF["Batch"]==2000]

    dataF.loc[AE, "error_dstd_x"] = abs(dataF.loc[AE, "delta_stddev_x"] - maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "delta_stddev_x"].values) / maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "delta_stddev_x"].values
    dataF.loc[CC, "error_dstd_x"] = abs(dataF.loc[CC, "delta_stddev_x"] - maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "delta_stddev_x"].values) / maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "delta_stddev_x"].values
    dataF.loc[EE, "error_dstd_x"] = abs(dataF.loc[EE, "delta_stddev_x"] - maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "delta_stddev_x"].values) / maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "delta_stddev_x"].values
    dataF.loc[AA, "error_dstd_x"] = abs(dataF.loc[AA, "delta_stddev_x"] - maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "delta_stddev_x"].values) / maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "delta_stddev_x"].values

    dataF.loc[AE, "error_dstd_y"] = abs(dataF.loc[AE, "delta_stddev_y"] - maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "delta_stddev_y"].values) / maxvals.loc[maxvals.Simulation == "A-E-HP5MM", "delta_stddev_y"].values
    dataF.loc[CC, "error_dstd_y"] = abs(dataF.loc[CC, "delta_stddev_y"] - maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "delta_stddev_y"].values) / maxvals.loc[maxvals.Simulation == "C-C-HP5MM", "delta_stddev_y"].values
    dataF.loc[EE, "error_dstd_y"] = abs(dataF.loc[EE, "delta_stddev_y"] - maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "delta_stddev_y"].values) / maxvals.loc[maxvals.Simulation == "E-E-HP5MM", "delta_stddev_y"].values
    dataF.loc[AA, "error_dstd_y"] = abs(dataF.loc[AA, "delta_stddev_y"] - maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "delta_stddev_y"].values) / maxvals.loc[maxvals.Simulation == "A-A-HP5MM", "delta_stddev_y"].values

    maxvals = dataF[dataF["Batch"]==2000]

    return dataF, AA, AE, CC, EE
