__author__ = 'boris'
"""
inputs:
    outfolder - where to put all the results
    control_file_name - a pickled dict of [strand][chromosome][position] = background-subtracted mutations/coverage
        output from normalize_to_control_make_wig.py - this is already a comparison of modifier to no modifier.
    experimental_file_names - any number of files, of same format as control file, to be compared to the control.
outputs:
    for each file, a WIG file with the max rescaled to 1
    ?add pseudocounts?
    for each experimental file, a WIG of the log2 ratio of experimental/control
    for each experimental file, a WIG of the log2 ratio of experimental/control, which were rescaled before comparison.
    for each experimental file, a pickled dictionary of the log2 ratio, in the same format as above
    for each experimental file, a pickled dictionary of the log2 ratio, previously rescaled as above
"""

import sys
import mod_utils
import os
import gzip
from scipy.stats.mstats import winsorize
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy
import math
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines

def write_wig(mutation_dict, sample_name, output_prefix):
    """

    :param mutation_dict: a pickled dict of [strand][chromosome][position] = mutations/coverage
    :param output_prefix:
    :return:
    """
    score_table = open(output_prefix+'_scores.txt', 'w')
    plusWig = gzip.open(output_prefix+'_plus.wig.gz', 'w')
    #minusWig = gzip.open(output_prefix+'_minus.wig.gz', 'w')
    plusWig.write('track type=wiggle_0 name=%s\n' % (sample_name+'_plus'))
    #minusWig.write('track type=wiggle_0 name=%s\n' % (sample_name+'_minus'))
    allChrs=set(mutation_dict['+'].keys()).union(set(mutation_dict['-'].keys()))
    for chr in allChrs:
        #minPos = min([min(weightedReads['+'][chr].keys()), min(weightedReads['-'][chr].keys())])
        #maxPos = max([max(weightedReads['+'][chr].keys()), max(weightedReads['-'][chr].keys())])
        plusWig.write('variableStep chrom=%s\n' % (chr))
        #minusWig.write('variableStep chrom=%s\n' % (chr))
        if chr in mutation_dict['+']:
            for i in sorted(mutation_dict['+'][chr].keys()):
                plusWig.write('%d\t%f\n' % (i, mutation_dict['+'][chr][i]))
                score_table.write('%s_%d\t%f\n' % (chr, i, mutation_dict['+'][chr][i]))
        #if chr in mutation_dict['-']:
            #for i in sorted(mutation_dict['-'][chr].keys()):
                #minusWig.write('%d\t%f\n' % (i, mutation_dict['-'][chr][i]/mutation_dict))
    plusWig.close()
    score_table.close()
    #minusWig.close()

def normalize_dict_to_max(mutation_dict, winsorize_data = False, winsorization_limits = (0, 0.95)):
    all_values = []
    normed_dict = {}
    for strand in mutation_dict:
        normed_dict[strand] = {}
        for chromosome in mutation_dict[strand]:
            normed_dict[strand][chromosome] = {}
            #print mutation_dict[strand][chromosome].values()
            all_values += mutation_dict[strand][chromosome].values()
            #print all_values
    if winsorize_data:
        winsorize(all_values, limits = (winsorization_limits[0], 1-winsorization_limits[1]), inplace = True)

    max_value = float(max(all_values))
    for strand in mutation_dict:
        for chromosome in mutation_dict[strand]:
            for position in mutation_dict[strand][chromosome]:
                val = mutation_dict[strand][chromosome][position]
                if val < min(all_values):
                    val = min(all_values)
                if val > max(all_values):
                    val = max(all_values)
                normed_dict[strand][chromosome][position] = val/max_value
    return normed_dict

def plot_weighted_nts_pie(background_subtracted, fasta_genome, title, out_prefix):
    genome = mod_utils.convertFastaToDict(fasta_genome)
    fig = plt.figure(figsize=(8,8))
    plot = fig.add_subplot(111)#a pie chart of mutated nts weighted by background-subtracted counts
    labels = "ATCG"
    nt_counts = defaultdict(float)
    for strand in background_subtracted:
        for chromosome in background_subtracted[strand]:
            for position in background_subtracted[strand][chromosome]:
                nt = genome[chromosome][position-1]
                nt_counts[nt] += background_subtracted[strand][chromosome][position]
    sizes = numpy.array([nt_counts[nt] for nt in labels])
    total = float(sum(sizes))
    sizes = sizes/total
    merged_labels = ['%s %.3f' % (labels[i], sizes[i]) for i in range(len(sizes))]
    plot.pie(sizes, labels = merged_labels, colors = mod_utils.rainbow)
    plot.set_title(title)

    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def compare_to_control(experiment_dict, control_dict):
    """

    :param experiment_dict:
    :param normalization_dict:
    :return: subtracted_dict - normalization value subtracted from experimental at every position
    """
    ratio_dict = {}
    for strand in experiment_dict:
        ratio_dict[strand] = {}
        for chromosome in experiment_dict[strand]:
            ratio_dict[strand][chromosome] = {}
            for position in experiment_dict[strand][chromosome]:
                try:
                    #dividing by zero will cause issues here
                    ratio_dict[strand][chromosome][position] = math.log(float(experiment_dict[strand][chromosome][position])/float(control_dict[strand][chromosome][position]), 2)
                except:
                    pass
    return ratio_dict

def normed_mutation_rate_histogram(normalized_mutations, dataset_names, output_prefix, title = '', xlim= (0,1), min = 0, max = 1, step = 0.01):
    fig = plt.figure(figsize=(16,16))
    plot = fig.add_subplot(111)
    bins = numpy.arange(xlim[0],xlim[1],step)
    bins = numpy.append(bins, max+step)
    bins = numpy.append([min-step], bins)
    for i in range(len(dataset_names)):
        mutation_densities = []
        for strand in normalized_mutations[i]:
            for chromosome in normalized_mutations[i][strand]:
                mutation_densities = mutation_densities + [val for val in normalized_mutations[i][strand][chromosome].values()]
        counts, edges = numpy.histogram(mutation_densities, bins = bins)
        #plot.hist(mutation_densities, color = mod_utils.colors  [i], bins = bins, label=dataset_names[i])
        counts = [0]+list(counts)+[0]
        edges = [0]+list(edges)+[edges[-1]]
        plot.fill(edges[:-1], numpy.array(counts), alpha = 0.3, color = mod_utils.colors[i], lw=0)
        plot.plot(edges[:-1], numpy.array(counts), color = mod_utils.colors[i], label=dataset_names[i], lw=2)
    plot.set_xlim(xlim[0]-step,xlim[1]+step)
    #plot.set_xticks(numpy.arange(0,10)+0.5)
    #plot.set_xticklabels(numpy.arange(0,10))
    plot.set_xlabel('mutations/coverage')
    plot.set_ylabel("# positions")
    plot.set_title(title)
    lg=plt.legend(loc=2,prop={'size':10}, labelspacing=0.2)
    lg.draw_frame(False)
    #plot.set_yscale('log')
    #plot.set_title(title)
    plt.savefig(output_prefix + '_mut_density.pdf', transparent='True', format='pdf')
    plt.clf()

def main():
    outfolder, genome_fasta, normalization_file_name = sys.argv[1:4]
    experimental_file_names = sys.argv[4:]

    control_dict = mod_utils.unPickle(normalization_file_name)
    rescaled_control_dict = normalize_dict_to_max(control_dict)
    norm_name = '.'.join(os.path.basename(normalization_file_name).split('.')[:-1])
    experimental_dict_names = ['.'.join(os.path.basename(file_name).split('.')[:-1]) for file_name in experimental_file_names]
    experimental_dicts = [mod_utils.unPickle(file_name) for file_name in experimental_file_names]
    rescaled_experimental_dicts = [normalize_dict_to_max(exp_dict) for exp_dict in experimental_dicts]

    print experimental_dict_names, norm_name

    normed_mutation_rate_histogram(rescaled_experimental_dicts, experimental_dict_names, os.path.join(outfolder, '%s_rescaled_mutation_rate_histogram' % norm_name), title='mutation rate, rescaled to max', xlim = (0, 0.1), min = 0, max =1, step = 0.001)
    comparisons = []
    rescaled_comparisons = []
    write_wig(control_dict, norm_name, os.path.join(outfolder, norm_name))
    for i in range(len(experimental_dict_names)):
        write_wig(rescaled_experimental_dicts[i], experimental_dict_names[i], os.path.join(outfolder, experimental_dict_names[i]))
        comparison_log2_ratios = compare_to_control(experimental_dicts[i], control_dict)
        rescaled_comparison_log2_ratios = compare_to_control(rescaled_experimental_dicts[i], rescaled_control_dict)
        comparisons.append(comparison_log2_ratios)
        rescaled_comparisons.append(rescaled_comparison_log2_ratios)
        mod_utils.makePickle(comparison_log2_ratios, os.path.join(outfolder, experimental_dict_names[i]+'_comparison_log2.pkl'))
        #mod_utils.makePickle(rescaled_comparison_log2_ratios, os.path.join(outfolder, experimental_dict_names[i]+'_rescaled_comparison_log2.pkl'))
        write_wig(comparison_log2_ratios, experimental_dict_names[i]+'_comparison_log2', os.path.join(outfolder, experimental_dict_names[i]+'_comparison_log2'))
        #write_wig(rescaled_comparison_log2_ratios, experimental_dict_names[i]+'_rescaled)comparison_log2', os.path.join(outfolder, experimental_dict_names[i]+'_rescaled_comparison_log2'))
        #try:
        #    plot_weighted_nts_pie(background_subtracted, genome_fasta, '%s backround-subtracted fractions' % experimental_dict_names[i], os.path.join(outfolder, experimental_dict_names[i]+'_sub_pie'))
        #except:
        #    pass
    #print comparisons
    #print rescaled_comparisons
    normed_mutation_rate_histogram(comparisons, experimental_dict_names, os.path.join(outfolder, '%s_comparison_histogram' % norm_name), title='log2 experiment/control', xlim = (-10, 10), min = -100, max =100, step = 0.1)
    #normed_mutation_rate_histogram(rescaled_comparisons, experimental_dict_names, os.path.join(outfolder, '%s_rescaled_comparison_histogram' % norm_name), title='log2 rescaled experiment/control', xlim = (-10, 10), min = -100, max =100, step = 0.1)



main()