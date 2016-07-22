# coding: utf-8
import matplotlib
import numpy as np
import seaborn as sb
import pandas as pd

from calc_stats import do_calcs

pyplot = sb.plt # Hacky, but gets the nice seaborn settings for things
params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [8.5, 4.5]
   }
pyplot.rcParams.update(params)

y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
mean_of_means_style = 'k^'
print "Running calculations..."

dataF, AA, AE, CC, EE = do_calcs()

def plot_group(data, filebase, function):
    locations = [AA, AE, CC, EE]
    colour_indices = [5, 0, 1, 2]
    names = ["A-A-HP5MM", "A-E-HP5MM", "C-C-HP5MM", "E-E-HP5MM"]
    for (l,c,n) in zip(locations, colour_indices, names):
        function(data, l, c, "%s_%s.eps" % (filebase, n))

def mean_x_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Mean X Position")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "mean_pos_x"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["mean_pos_x"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    pyplot.legend(['Sample means', 'Mean of Means'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def mean_y_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Mean Y Position")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "mean_pos_y"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["mean_pos_y"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    pyplot.legend(['Sample means', 'Mean of Means'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def error_x_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Error in X Mean")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "rel_mpx_err"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["rel_mpx_err"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    a,b = pyplot.ylim()
    if b < 0.003:
        b = 0.003
    pyplot.ylim([-0.05 * b, b])
    pyplot.legend(['Sample error', 'Mean error'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def error_y_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Error in Y Mean")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "rel_mpy_err"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["rel_mpy_err"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    a,b = pyplot.ylim()
    if b < 0.00035:
        b = 0.00035
    pyplot.ylim([-0.05 * b, b])
    pyplot.legend(['Sample error', 'Mean error'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def delta_mpx_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Change in X Mean")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "delta_mpx"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["delta_mpx"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    pyplot.legend(['Sample', 'Mean'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def delta_mpy_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Change in Y Mean")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "delta_mpy"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["delta_mpy"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    pyplot.legend(['Sample', 'Mean'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def error_dmpx_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Relative Error")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "error_mpx"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["error_mpx"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    pyplot.ylim([-0.01, 0.25])
    pyplot.legend(['Sample', 'Mean'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def error_dmpy_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Relative Error")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "error_mpy"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    gm=dat[selection].groupby("Batch").mean()["error_mpy"]
    gm=gm[gm.index.values <= 1000]
    pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    pyplot.ylim([-0.01, 0.35])
    pyplot.legend(['Sample', 'Mean'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def stddev_x_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Standard Deviation")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "stddev_pos_x"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    #gm=dat[selection].groupby("Batch").mean()["stddev_pos_x"]
    #gm=gm[gm.index.values <= 1000]
    #pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    #pyplot.legend(['Sample means', 'Mean of Means'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def stddev_y_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Standard Deviation")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "stddev_pos_y"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    #gm=dat[selection].groupby("Batch").mean()["stddev_pos_y"]
    #gm=gm[gm.index.values <= 1000]
    #pyplot.plot(gm.index, gm, mean_of_means_style)
    pyplot.xlim([0,2100])
    #pyplot.legend(['Sample means', 'Mean of Means'])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def stddev_err_x_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Relative Error")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "error_stddev_x"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    pyplot.ylim([-0.001, 0.08])
    pyplot.xlim([0,2100])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def stddev_err_y_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Relative Error")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "error_stddev_y"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    pyplot.xlim([0,2100])
    pyplot.ylim([-0.003, 0.15])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def d_stddev_x_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Change in Standard Deviation")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "delta_stddev_x"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    pyplot.xlim([0,2100])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def d_stddev_y_plot(dat, selection, colour_index, filename):
    pyplot.figure()
    pyplot.gca().yaxis.set_major_formatter(y_formatter)
    pyplot.xlabel("Population Size")
    pyplot.ylabel("Change in Standard Deviation")
    pyplot.plot(dat.loc[selection, "Batch"], dat.loc[selection, "delta_stddev_y"], 'o', color=pyplot.rcParams['axes.color_cycle'][colour_index])
    pyplot.xlim([0,2100])
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def d_stddev_xy_plot(data, selections, colour_indices, legends, filename):
    pyplot.figure()
    for (l,c) in zip(selections, colour_indices):
        pyplot.plot(data.loc[l, "delta_stddev_x"], data.loc[l, "delta_stddev_y"], ".", color=pyplot.rcParams['axes.color_cycle'][c])
    pyplot.xlabel("Change in Standard Deviation (in X)")
    pyplot.ylabel("Change in Standard Deviation (in Y)")
    pyplot.legend(legends, loc='best', fancybox=True)
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def d_mean_xy_plot(data, selections, colour_indices, legends, filename):
    pyplot.figure()
    for (l,c) in zip(selections, colour_indices):
        pyplot.plot(data.loc[l, "delta_mpx"], data.loc[l, "delta_mpy"], ".", color=pyplot.rcParams['axes.color_cycle'][c])
    pyplot.xlabel("Change in Mean Position (in X)")
    pyplot.ylabel("Change in Mean Position (in Y)")
    pyplot.legend(legends, loc='best', fancybox=True)
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

print "Calculations complete. Plotting..."

plot_group(dataF, "mean_x", mean_x_plot)
plot_group(dataF, "mean_y", mean_y_plot)

plot_group(dataF, "error_x", error_x_plot)
plot_group(dataF, "error_y", error_y_plot)

plot_group(dataF, "mdisp_x", delta_mpx_plot)
plot_group(dataF, "mdisp_y", delta_mpy_plot)

plot_group(dataF, "error_dmpx", error_dmpx_plot)
plot_group(dataF, "error_dmpy", error_dmpy_plot)

plot_group(dataF, "stddev_x", stddev_x_plot)
plot_group(dataF, "stddev_y", stddev_y_plot)

plot_group(dataF, "stddev_err_x", stddev_err_x_plot)
plot_group(dataF, "stddev_err_y", stddev_err_y_plot)

plot_group(dataF, "delta_stddev_x", d_stddev_x_plot)
plot_group(dataF, "delta_stddev_y", d_stddev_y_plot)

d_stddev_xy_plot(dataF, [AA, CC, EE], [5, 1, 2],["A-A-HP5MM", "C-C-HP5MM", "E-E-HP5MM"], "delta_stddev_xy.eps")
d_stddev_xy_plot(dataF, [CC, EE], [1, 2], ["C-C-HP5MM", "E-E-HP5MM"], "delta_stddev_xy2.eps")


d_mean_xy_plot(dataF, [AA, CC, EE], [5, 1, 2],["A-A-HP5MM", "C-C-HP5MM", "E-E-HP5MM"], "delta_mean_xy.eps")
d_mean_xy_plot(dataF, [CC, EE], [1, 2], ["C-C-HP5MM", "E-E-HP5MM"], "delta_mean_xy2.eps")

print "Done."
print

