__author__ = 'boris'
from scipy.stats.mstats import winsorize
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy
import math
import mod_lib
import mod_utils
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines


def plot_mutated_nts_pie(libraries, out_prefix, subtract_background = False):
    #Makes an array of pie charts, 1 per library
    if subtract_background:
        #if subtracting background, need to only look at those which have a defined control
        libraries = [library for library in libraries if library.lib_settings.sample_name in
                     library.experiment_settings.get_property('experimentals')]
    num_subplots = len(libraries)
    num_plots_wide = math.ceil(math.sqrt(num_subplots))
    num_plots_high = num_plots_wide
    fig = plt.figure(figsize=(4*num_plots_wide, 4*num_plots_high))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plot_index =1
    for library in libraries:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
        mutated_nts_count = library.count_mutations_by_nucleotide(subtract_background = subtract_background)
        labels = sorted(mutated_nts_count.keys())
        sizes = [mutated_nts_count[nt] for nt in labels]
        total = float(sum(sizes))
        merged_labels = ['%s %.3f' % (labels[i], sizes[i]/total) for i in range(len(sizes))]
        plot.pie(sizes, labels = merged_labels, colors = mod_utils.rainbow)
        plot.set_title(library.lib_settings.sample_name)
        plot_index += 1
    if subtract_background:
        plt.suptitle('background-subtracted fractions of mutated nts')
    else:
        plt.suptitle('fractions of mutated nts')
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()
