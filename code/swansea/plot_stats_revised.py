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

def d_stddev_xy_plot(data, selections, colour_indices, legends, filename):
    pyplot.figure()
    for (l,c) in zip(selections, colour_indices):
        pyplot.scatter(data.loc[l, "delta_stddev_x"], data.loc[l, "delta_stddev_y"], c=data.loc[l, "Batch"], color=pyplot.rcParams['axes.color_cycle'][c], cmap=sb.dark_palette(pyplot.rcParams['axes.color_cycle'][c], as_cmap=True, reverse=True))
    pyplot.xlabel("Change in Standard Deviation (in X)")
    pyplot.ylabel("Change in Standard Deviation (in Y)")
    pyplot.legend(legends, loc='upper left', fancybox=True)
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

def d_mean_xy_plot(data, selections, colour_indices, legends, filename):
    pyplot.figure()
    for (l,c) in zip(selections, colour_indices):
        pyplot.scatter(data.loc[l, "delta_mpx"], data.loc[l, "delta_mpy"], c=data.loc[l, "Batch"], color=pyplot.rcParams['axes.color_cycle'][c], cmap=sb.dark_palette(pyplot.rcParams['axes.color_cycle'][c], as_cmap=True, reverse=True))
    pyplot.xlabel("Change in Mean Position (in X)")
    pyplot.ylabel("Change in Mean Position (in Y)")
    pyplot.legend(legends, loc='upper left', fancybox=True)
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close()

print "Calculations complete. Plotting..."

d_stddev_xy_plot(dataF, [AA, CC, EE], [5, 1, 2],["A-A-HP5MM", "C-C-HP5MM", "E-E-HP5MM"], "delta_stddev_xy.new.eps")
d_stddev_xy_plot(dataF, [CC, EE], [1, 2], ["C-C-HP5MM", "E-E-HP5MM"], "delta_stddev_xy2.new.eps")


d_mean_xy_plot(dataF, [AA, CC, EE], [5, 1, 2],["A-A-HP5MM", "C-C-HP5MM", "E-E-HP5MM"], "delta_mean_xy.new.eps")
d_mean_xy_plot(dataF, [CC, EE], [1, 2], ["C-C-HP5MM", "E-E-HP5MM"], "delta_mean_xy2.new.eps")

print "Done."
print

