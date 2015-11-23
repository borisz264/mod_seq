__author__ = 'boris'
"""
inputs:
    outfolder - where to put all the results
    normalization_file_name - a pickled dict of [strand][chromosome][position] = mutations/coverage
        output from count_reads_and_mismatches.py - this is from a sample where no modifying reagent was added.
    experimental_file_names - any number of files, of same format as normalization file, to be normalized by the normalization file
outputs:
    for each file, a WIG file of the mutation rate, normalized by the coverage within the same sample
        set maximum to 1
        winsorize 0-95%?
    for each experimental file, a WIG of the coverage-normalized mutation rate, normalized again by the control.
    for each experimental file, a pickled dictionary of the coverage-normalized mutation rate, normalized again by the control, in the same format as above
"""

import sys
import mod_utils
import os
import gzip
from scipy.stats.mstats import winsorize
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines

def write_wig(mutation_dict, sample_name, output_prefix):
    """

    :param mutation_dict: a pickled dict of [strand][chromosome][position] = mutations/coverage
    :param output_prefix:
    :return:
    """

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
        #if chr in mutation_dict['-']:
            #for i in sorted(mutation_dict['-'][chr].keys()):
                #minusWig.write('%d\t%f\n' % (i, mutation_dict['-'][chr][i]/mutation_dict))
    plusWig.close()
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

def subtract_background(experiment_dict, normalization_dict):
    """

    :param experiment_dict:
    :param normalization_dict:
    :return: subtracted_dict - normalization value subtracted from experimental at every position
    """
    subtracted_dict = {}
    for strand in experiment_dict:
        subtracted_dict[strand] = {}
        for chromosome in experiment_dict[strand]:
            subtracted_dict[strand][chromosome] = {}
            for position in experiment_dict[strand][chromosome]:
                subtracted_dict[strand][chromosome][position] = max(experiment_dict[strand][chromosome][position]-normalization_dict[strand][chromosome][position], 0)
    return subtracted_dict

def main():
    outfolder, genome_fasta, normalization_file_name = sys.argv[1:4]
    experimental_file_names = sys.argv[4:]

    normalization_dict = mod_utils.unPickle(normalization_file_name)
    norm_name = '.'.join(os.path.basename(normalization_file_name).split('.')[:-2])
    experimental_dict_names = ['.'.join(os.path.basename(file_name).split('.')[:-2]) for file_name in experimental_file_names]
    experimental_dicts = [mod_utils.unPickle(file_name) for file_name in experimental_file_names]

    write_wig(normalization_dict, norm_name, os.path.join(outfolder, norm_name))
    for i in range(len(experimental_dict_names)):
        write_wig(experimental_dicts[i], experimental_dict_names[i], os.path.join(outfolder, experimental_dict_names[i]))
        background_subtracted = subtract_background(experimental_dicts[i], normalization_dict)
        mod_utils.makePickle(background_subtracted, os.path.join(outfolder, experimental_dict_names[i]+'_subtracted.pkl'))
        write_wig(background_subtracted, experimental_dict_names[i]+'_subtracted', os.path.join(outfolder, experimental_dict_names[i]+'_subtracted'))
        try:
            plot_weighted_nts_pie(background_subtracted, genome_fasta, '%s backround-subtracted fractions' % experimental_dict_names[i], os.path.join(outfolder, experimental_dict_names[i]+'_sub_pie'))
        except:
            pass
main()