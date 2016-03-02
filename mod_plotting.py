__author__ = 'boris'
from scipy.stats.mstats import winsorize
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy
import math
import mod_lib
import mod_utils
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines


def plot_mutated_nts_pie(libraries, out_prefix, subtract_background=False, subtract_control=False, exclude_constitutive=False):
    #Makes an array of pie charts, 1 per library
    if subtract_background or subtract_control:
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
        mutated_nts_count = library.count_mutation_rates_by_nucleotide(subtract_background=subtract_background, subtract_control=subtract_control,
                                                                       exclude_constitutive=exclude_constitutive)
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
    if subtract_control:
        plt.suptitle('control-subtracted mutation rate fractions')
    else:
        plt.suptitle('mutation rate fractions')
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def plot_mutation_rate_cdfs(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False):
    #Makes 2 CDF plots. One of all libraries, showing the coverage-normalized mutation rates
    # and one showing background-subtracted mutation rates

    fig = plt.figure(figsize=(24, 16))
    plots = []
    plot = fig.add_subplot(231)
    plots.append(plot)
    colormap = plt.get_cmap('spectral')
    colorindex = 0
    for library in libraries:
        all_mutation_rates = [val for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("mutation rate")
    plot.set_title('raw mutation rates')
    lg=plt.legend(loc=4,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-0.001, 0.02)

    plot = fig.add_subplot(232)
    plots.append(plot)
    colorindex = 0
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [val for val in
                              library.list_mutation_rates(subtract_background=True, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("background-subtracted mutation rate")
    plot.set_title('normalized mutation rates')
    lg=plt.legend(loc=4,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-0.001, 0.02)

    plot = fig.add_subplot(233)
    plots.append(plot)
    colorindex = 0
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [val for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=True,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive)]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("control-subtracted mutation rate")
    plot.set_title('control normalized mutation rates')
    lg=plt.legend(loc=4,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-0.001, 0.02)

    plot = fig.add_subplot(234)
    plots.append(plot)
    colormap = plt.get_cmap('spectral')
    colorindex = 0
    for library in libraries:
        all_mutation_rates = [math.log(val, 10) for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive) if val>0]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("log10 mutation rate")
    plot.set_title('raw mutation rates')
    lg=plt.legend(loc=2,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-5, -1)

    plot = fig.add_subplot(235)
    plots.append(plot)
    colorindex = 0
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [math.log(val, 10) for val in
                              library.list_mutation_rates(subtract_background=True, subtract_control=False,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive) if val>0]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("background-subtracted log10 mutation rate")
    plot.set_title('normalized mutation rates')
    lg=plt.legend(loc=2,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-5, -1)

    plot = fig.add_subplot(236)
    plots.append(plot)
    colorindex = 0
    libraries_to_plot = [library for library in libraries if library.lib_settings.sample_name in
             library.experiment_settings.get_property('experimentals')]
    for library in libraries_to_plot:
        all_mutation_rates = [math.log(val, 10) for val in
                              library.list_mutation_rates(subtract_background=False, subtract_control=True,
                                                          nucleotides_to_count=nucleotides_to_count, exclude_constitutive=exclude_constitutive) if val>0]
        plot.hist(all_mutation_rates, 10000, normed=1, cumulative=True, histtype='step', color=colormap(colorindex/float(len(libraries))),
                  label=library.lib_settings.sample_name, lw=2)
        colorindex += 1
    plot.set_xlabel("control-subtracted log10 mutation rate")
    plot.set_title('control normalized mutation rates')
    lg=plt.legend(loc=2,prop={'size':6}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xlim(-5, -1)

    for plot in plots:
        plot.set_ylabel("cumulative fraction of %s nucleotides" % (nucleotides_to_count))
        plot.set_ylim(0, 1)
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def plot_changes_vs_control_interactive(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False,
                                        max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library use bokeh to plot an interactive plot of magnitude of change (experimental-control)
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
            All nucleotides will be labelled on mouseover.
    """
    from bokeh.plotting import figure, output_file, show, ColumnDataSource, gridplot
    from bokeh.models import Range1d
    from bokeh.models import HoverTool
    from collections import OrderedDict

    # output to static HTML file
    output_file("%s.html" % (out_prefix))
    plot_figs=[]

    for library in libraries:
        mag_change, fold_change, annotation = [], [], []
        prot_mag_change, prot_fold_change, prot_annotation = [], [], []
        deprot_mag_change, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        print 'deprotected ', nucleotide
                        deprot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        print 'protected ', nucleotide
                        prot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        source = ColumnDataSource(data=dict(x = mag_change, y = fold_change, label = annotation))
        prot_source = ColumnDataSource(data=dict(x = prot_mag_change, y = prot_fold_change, label = prot_annotation))
        deprot_source = ColumnDataSource(data=dict(x = deprot_mag_change, y = deprot_fold_change,
                                                   label = deprot_annotation))
        TOOLS = "pan,wheel_zoom,reset,save,hover"
        PlotFig = figure(x_axis_label = "mutation rate [%s] - [%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_label = "fold change [%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_type="log", tools=TOOLS, toolbar_location="right")
        PlotFig.circle("x", "y", size = 5, source=source, color=mod_utils.bokeh_black)
        PlotFig.circle("x", "y", size = 5, source=prot_source, color=mod_utils.bokeh_vermillion)
        PlotFig.circle("x", "y", size = 5, source=deprot_source, color=mod_utils.bokeh_bluishGreen)
        PlotFig.x_range = Range1d(start=-0.2, end=0.2)
        PlotFig.y_range = Range1d(start=.001, end=100)

        #adjust what information you get when you hover over it
        Hover = PlotFig.select(dict(type=HoverTool))
        Hover.tooltips = OrderedDict([("nuc", "@label")])
        plot_figs.append([PlotFig])
    p = gridplot(plot_figs)
    show(p)

def ma_plots_interactive(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False,
                         max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library use bokeh to plot an interactive plot of average magnitude of signal (experimental+control)/2
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
            All nucleotides will be labelled on mouseover.
    """
    from bokeh.plotting import figure, output_file, show, ColumnDataSource, gridplot
    from bokeh.models import Range1d
    from bokeh.models import HoverTool
    from collections import OrderedDict

    # output to static HTML file
    output_file("%s.html" % (out_prefix))
    plot_figs=[]

    for library in libraries:
        mag, fold_change, annotation = [], [], []
        prot_mag, prot_fold_change, prot_annotation = [], [], []
        deprot_mag, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    avg_mutation_rate = (nucleotide.mutation_rate+nucleotide.get_control_nucleotide().mutation_rate)/2.0
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag.append(avg_mutation_rate)
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        print 'deprotected ', nucleotide
                        deprot_mag.append(avg_mutation_rate)
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        print 'protected ', nucleotide
                        prot_mag.append(avg_mutation_rate)
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        source = ColumnDataSource(data=dict(x = mag, y = fold_change, label = annotation))
        prot_source = ColumnDataSource(data=dict(x = prot_mag, y = prot_fold_change, label = prot_annotation))
        deprot_source = ColumnDataSource(data=dict(x = deprot_mag, y = deprot_fold_change,
                                                   label = deprot_annotation))
        TOOLS = "pan,wheel_zoom,reset,save,hover"
        PlotFig = figure(x_axis_label = "avg signal ([%s] + [%s])/2" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_label = "fold change [%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),
                         y_axis_type="log", x_axis_type="log", tools=TOOLS, toolbar_location="right")
        PlotFig.circle("x", "y", size = 5, source=source, color=mod_utils.bokeh_black)
        PlotFig.circle("x", "y", size = 5, source=prot_source, color=mod_utils.bokeh_vermillion)
        PlotFig.circle("x", "y", size = 5, source=deprot_source, color=mod_utils.bokeh_bluishGreen)
        PlotFig.x_range = Range1d(start=0.00001, end=1)
        PlotFig.y_range = Range1d(start=.001, end=100)

        #adjust what information you get when you hover over it
        Hover = PlotFig.select(dict(type=HoverTool))
        Hover.tooltips = OrderedDict([("nuc", "@label")])
        plot_figs.append([PlotFig])
    p = gridplot(plot_figs)
    show(p)

def plot_changes_vs_control(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library make a plot of magnitude of change (experimental-control)
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
    """
    output_file = "%s.pdf" % (out_prefix)
    plot_figs=[]

    num_subplots = len(libraries)
    num_plots_wide = math.ceil(math.sqrt(num_subplots))
    num_plots_high = num_plots_wide
    fig = plt.figure(figsize=(4*num_plots_wide, 4*num_plots_high))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plot_index =1
    for library in libraries:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
        mag_change, fold_change, annotation = [], [], []
        prot_mag_change, prot_fold_change, prot_annotation = [], [], []
        deprot_mag_change, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        print 'deprotected ', nucleotide
                        deprot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        print 'protected ', nucleotide
                        prot_mag_change.append(nucleotide.get_control_sub_mutation_rate())
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        plot.set_xlabel("[%s] - [%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name), fontsize = 8)
        plot.set_ylabel("[%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name),  fontsize = 8)
        plot.set_yscale('log')
        plot.scatter(mag_change, fold_change, color=mod_utils.black, s=2)
        plot.scatter(prot_mag_change, prot_fold_change, color=mod_utils.vermillion, s=2)
        plot.scatter(deprot_mag_change, deprot_fold_change, color=mod_utils.bluishGreen, s=2)
        plot.set_xlim(-0.2,0.2)
        plot.set_ylim(.001,100)

        plot_figs.append(plot)
        plot_index+=1
    plt.savefig(output_file, transparent='True', format='pdf')

def ma_plots(libraries, out_prefix, nucleotides_to_count='ATCG', exclude_constitutive=False,
             max_fold_reduction=0.001, max_fold_increase=100):
    """

    :param libraries:
    :param out_prefix:
    :param nucleotides_to_count:
    :param exclude_constitutive:
    :return: for each library use bokeh to plot an interactive plot of magnitude of signal (experimental+control)/2
            vs log10 fold change (experimental/control).
            Protected and de-protected calls will be colored, based on a fold change cutoff and confidence interval.
            All nucleotides will be labelled on mouseover.
    """
    output_file = "%s.pdf" % (out_prefix)
    plot_figs=[]

    num_subplots = len(libraries)
    num_plots_wide = math.ceil(math.sqrt(num_subplots))
    num_plots_high = num_plots_wide
    fig = plt.figure(figsize=(4*num_plots_wide, 4*num_plots_high))
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plot_index =1
    for library in libraries:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
        mag, fold_change, annotation = [], [], []
        prot_mag, prot_fold_change, prot_annotation = [], [], []
        deprot_mag, deprot_fold_change, deprot_annotation = [], [], []
        for rRNA_name in library.rRNA_mutation_data:
            for position in library.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = library.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if (exclude_constitutive and nucleotide.exclude_constitutive)or nucleotide.identity not in nucleotides_to_count:
                    pass
                else:
                    protection_call = nucleotide.determine_protection_status(confidence_interval=library.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                   fold_change_cutoff=library.experiment_settings.get_property('fold_change_cutoff'))
                    control_fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                    avg_mutation_rate = (nucleotide.mutation_rate+nucleotide.get_control_nucleotide().mutation_rate)/2.0
                    if control_fold_change == 0:
                        control_fold_change = max_fold_reduction
                    elif control_fold_change == float('inf'):
                        control_fold_change = max_fold_increase
                    if protection_call == 'no_change':
                        mag.append(avg_mutation_rate)
                        fold_change.append(control_fold_change)
                        annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'deprotected':
                        print 'deprotected ', nucleotide
                        deprot_mag.append(avg_mutation_rate)
                        deprot_fold_change.append(control_fold_change)
                        deprot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
                    elif protection_call == 'protected':
                        print 'protected ', nucleotide
                        prot_mag.append(avg_mutation_rate)
                        prot_fold_change.append(control_fold_change)
                        prot_annotation.append('%s_%s%d' %(rRNA_name,nucleotide.identity,position))
        plot.set_xlabel("([%s] + [%s])/2" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name), fontsize = 8)
        plot.set_ylabel("[%s]/[%s]" % (library.lib_settings.sample_name, library.get_normalizing_lib_with_mod().lib_settings.sample_name), fontsize = 8)
        plot.set_yscale('log')
        plot.set_xscale('log')
        plot.scatter(mag, fold_change, color=mod_utils.black, s=2)
        plot.scatter(prot_mag, prot_fold_change, color=mod_utils.vermillion, s=2)
        plot.scatter(deprot_mag, deprot_fold_change, color=mod_utils.bluishGreen, s=2)
        plot.set_xlim(0.00001,1)
        plot.set_ylim(.001,100)
        plot_figs.append(plot)
        plot_index+=1
    plt.savefig(output_file, transparent='True', format='pdf')