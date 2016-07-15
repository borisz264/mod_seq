__author__ = 'boris'
"""
THIS IS FOR TROUBLESHOOTING AND COMPARING DIFFERENT TRUE POSITIVE AND TRUE NEGATIVE ANNOTATIONS


5'e end data is a pickled dict of form srt_dict[strand][chrom][position] = counts at position
take the 5' end data from count_reads_and_mismatches.py, as well as any number of files output by
    compute_true_positive_negative.py

and compute:
    1) 90% windorize the input data (All data above 95th percentile set to 95th percentile)
    2) normalize to highest position in rRNA (should just stick to the rRNA  of interest, which will be 18S for my initial test)
    3) slide a cutoff from 0 to 1 in ~10,000 steps, computing % of true positives and true negatives called positive at each step
    4) plot these two percentages against each other for each step, also output these values as a spreadsheet
        also plot y=x for reference
"""

import sys, mod_utils, os
import numpy
from scipy.stats.mstats import winsorize
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
from collections import defaultdict

def winsorize_norm_chromosome_data(read_5p_ends, chromosome, strand, genome_dict, nucs_to_count, to_winsorize = True, low = 0, high = 0.95):
    """


    :param read_5p_ends:
    :param chromosome:
    :param strand:
    :param genome_dict:
    :param nucs_to_count:
    :param low:
    :param high:
    :return: an array (now zero-indexed from 1-indexed) of densities for the given chromosome on the given strand, winsorized, and only for the given nucleotides
    """
    max_position = max(read_5p_ends[strand][chromosome].keys())
    density_array =numpy.array([0] * max_position)
    for position in read_5p_ends[strand][chromosome].keys():
        if genome_dict[chromosome][position-1] in nucs_to_count:
            density_array[position-1] = read_5p_ends[strand][chromosome][position]
    if to_winsorize:
        winsorize(density_array, limits = (low, 1-high), inplace = True)
    normed_array = density_array/float(max(density_array))
    return  normed_array

def get_tp_tn(tp_tn_file):
    TP = set()
    TN = set()
    f = open(tp_tn_file)
    for line in f:
        ll= line.strip('\n').split('\t')
        if ll[2] == 'TP':
            TP.add(int(ll[0]))
        if ll[2] =='TN':
            TN.add(int(ll[0]))
    f.close()
    return TP, TN

def call_positives(density_array, chromosome, strand, genome_dict, nucs_to_count, cutoff):
    """

    :param density_array:
    :return:a set of called positive positions
            I've reverted these to 1-indexed to match the TP and TN calls from the structures
    """
    positives = set()

    for i in range(len(density_array)):
        if genome_dict[chromosome][i-1] in nucs_to_count:
            if density_array[i] >= cutoff:
                positives.add(i)#adding 1 not necessary, since the modified nucleotide is the one 1 upstream of the RT stop!!!

    return positives

def plot_ROC_curves(roc_curves, out_prefix):
    fig = plt.figure(figsize=(8,8))
    plot = fig.add_subplot(111)#first a pie chart of mutated nts
    color_index = 0
    for name in roc_curves:
        x, y = roc_curves[name]
        plot.plot(x, y, lw =2, label = name, color = mod_utils.rainbow[color_index])
        color_index +=1
    plot.plot(x, x, lw =1, ls = 'dashed', color = mod_utils.rainbow[color_index], label = 'y=x')
    plot.set_xlabel('False positive rate (%) (100-specificity)')
    plot.set_ylabel('True positive rate (%) (sensitivity)')
    lg=plt.legend(loc=2,prop={'size':10}, labelspacing=0.2)
    lg.draw_frame(False)
    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def pie_read_5p_ends(read_5p_ends, genome_dict, out_prefix):
    nuc_counts = defaultdict(int)
    for chromosome in read_5p_ends['+']:
        for position in read_5p_ends['+'][chromosome]:
            if position-2 > 0 :
                nuc = genome_dict[chromosome][position-1]
                nuc_counts[nuc] += read_5p_ends['+'][chromosome][position]
    fig = plt.figure(figsize=(8,8))

    plot = fig.add_subplot(111)#first a pie chart of mutated nts
    labels = sorted(nuc_counts.keys())
    sizes = [nuc_counts[nt] for nt in labels]
    plot.pie(sizes, labels = labels, colors = mod_utils.rainbow)
    plot.set_title('nt exactly at read 5p ends across rRNA')

    plt.savefig(out_prefix + '_nt_5p_ends.pdf', transparent='True', format='pdf')
    plt.clf()
def main():
    read_5p_ends_file, genome_fasta, outprefix = sys.argv[1:4]
    tp_tn_annotations = sys.argv[4:]#true positive and true negative annotations
    genome_dict = mod_utils.convertFastaToDict(genome_fasta)
    read_5p_ends = mod_utils.unPickle(read_5p_ends_file)
    normed_density_array = winsorize_norm_chromosome_data(read_5p_ends, 'S.c.18S_rRNA', '+', genome_dict, 'ACTG')
    real_tp_tn_data = []
    for filename in tp_tn_annotations:
        real_tp, real_tn = get_tp_tn(filename)
        real_tp_tn_data.append((os.path.basename(filename), real_tp, real_tn))

    roc_curves = {}
    for entry in real_tp_tn_data:
        roc_curves[entry[0]] = [[],[]]#x and y value arrays for each

    stepsize = 0.0001
    for cutoff in numpy.arange(0,1.+5*stepsize, stepsize):
        called_p = call_positives(normed_density_array, 'S.c.18S_rRNA', '+', genome_dict, 'AC', cutoff)
        for entry in real_tp_tn_data:
            #print called_p.intersection(entry[1])

            num_tp_called = len(called_p.intersection(entry[1]))#how many true positives called at this cutoff
            num_fp_called = len(called_p.intersection(entry[2]))#how many fp positives called at this cutoff
            roc_curves[entry[0]][0].append(100.*num_fp_called/float(len(entry[2])))#FP rate on x axis
            roc_curves[entry[0]][1].append(100.*num_tp_called/float(len(entry[1])))#TP rate on y axis

    plot_ROC_curves(roc_curves, outprefix)
    #pie_read_5p_ends(read_5p_ends, genome_dict, outprefix)
if __name__ == '__main__':
    main()