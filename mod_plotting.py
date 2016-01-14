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
        mutated_nts_count = library.count_mutation_rates_by_nucleotide(subtract_background = subtract_background)
        labels = sorted(mutated_nts_count.keys())
        sizes = numpy.array([mutated_nts_count[nt] for nt in labels])
        total = float(sum(sizes))
        sizes = sizes/total
        merged_labels = ['%s %.3f' % (labels[i], sizes[i]) for i in range(len(sizes))]
        plot.pie(sizes, labels = merged_labels, colors = mod_utils.rainbow)
        plot.set_title(library.lib_settings.sample_name)
        plot_index += 1
    if subtract_background:
        plt.suptitle('background-subtracted mutation rate fractions')
    else:
        plt.suptitle('mutation rate fractions')
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def plot_mutation_rate_cdfs(libraries, out_prefix):
    #Makes 2 CDF plots. One of all libraries, showing the coverage-normalized mutation rates
    # and one showing background-subtracted mutation rates

    fig = plt.figure(figsize=(16,8))
    plot = fig.add_subplot(121)
    colormap = plt.get_cmap('jet')
    colorindex = 0
    for library in libraries:
        all_mutation_rates = [math.log(val, 10) for val in  library.list_mutation_rates(subtract_background=False)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("log10 mutation rate")
    plot.set_title('raw mutation rates')
    plot = fig.add_subplot(122)
    colorindex = 0
    libraries = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries:
        all_mutation_rates = [math.log(val, 10) for val in  library.list_mutation_rates(subtract_background=True)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("background-subtracted log10 mutation rate")
    plot.set_title('normalized mutation rates')

    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()
