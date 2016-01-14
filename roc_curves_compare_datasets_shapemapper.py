__author__ = 'boris'
"""

5'e end data is a pickled dict of form srt_dict[chrom][position] = counts at position
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

def winsorize_norm_chromosome_data(mut_density, chromosome, genome_dict, nucs_to_count, to_winsorize = False, low = 0, high = 0.95):
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
    max_position = max(mut_density[chromosome].keys())
    density_array =numpy.array([0.0] * max_position)
    for position in mut_density[chromosome].keys():
        if genome_dict[chromosome][position-1] in nucs_to_count:
            density_array[position-1] = mut_density[chromosome][position]
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

def call_positives(density_array, chromosome, genome_dict, nucs_to_count, cutoff):
    """

    :param density_array:
    :return:a set of called positive positions
            I've reverted these to 1-indexed to match the TP and TN calls from the structures
    """
    positives = set()

    for i in range(len(density_array)):
        if genome_dict[chromosome][i] in nucs_to_count:
            if density_array[i] >= cutoff:
                positives.add(i+1)#adding 1 not necessary for RT stops, since the modified nucleotide is the one 1 upstream of the RT stop!!!

    return positives

def plot_ROC_curves(roc_curves, out_prefix):
    fig = plt.figure(figsize=(8,8))
    plot = fig.add_subplot(111)#first a pie chart of mutated nts
    colormap = plt.get_cmap('spectral')
    color_index = 0
    for name in sorted(roc_curves.keys()):
        x, y = roc_curves[name]
        area_under_curve = numpy.trapz(numpy.array(y[::-1])/100., x=numpy.array(x[::-1])/100.)
        plot.plot(x, y, lw =2, label = '%s   %.3f' % (name, area_under_curve), color = colormap(color_index/float(len(roc_curves))))
        color_index +=1
    plot.plot(numpy.arange(0,100,0.1), numpy.arange(0,100,0.1), lw =1, ls = 'dashed', color = mod_utils.black, label = 'y=x')
    plot.set_xlabel('False positive rate (%) (100-specificity)')
    plot.set_ylabel('True positive rate (%) (sensitivity)')
    lg=plt.legend(loc=4,prop={'size':10}, labelspacing=0.2)
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
    tp_tn_annotations, genome_fasta, outprefix = sys.argv[1:4]
    density_files = sys.argv[4:]
    sample_names = [os.path.basename(filename) for filename in density_files]
    mutation_densities = [mod_utils.unPickle(pickled_density) for pickled_density in density_files]

    genome_dict = mod_utils.convertFastaToDict(genome_fasta)
    normed_density_arrays = [winsorize_norm_chromosome_data(mutation_density, 'S.c.18S_rRNA', genome_dict, 'AC') for mutation_density in mutation_densities]
    real_tp, real_tn = get_tp_tn(tp_tn_annotations)
    roc_curves = {}
    for sample_name in sample_names:
        roc_curves[sample_name] = [[],[]]#x and y value arrays for each

    stepsize = 0.0001
    for cutoff in numpy.arange(0,1.+5*stepsize, stepsize):
        for i in range(len(sample_names)):
            called_p = call_positives(normed_density_arrays[i], 'S.c.18S_rRNA', genome_dict, 'AC', cutoff)
            num_tp_called = len(called_p.intersection(real_tp))#how many true positives called at this cutoff
            num_fp_called = len(called_p.intersection(real_tn))#how many fp positives called at this cutoff
            roc_curves[sample_names[i]][1].append(100.*num_tp_called/float(len(real_tp)))#TP rate on y axis
            roc_curves[sample_names[i]][0].append(100.*num_fp_called/float(len(real_tn)))#FP rate on x axis

    plot_ROC_curves(roc_curves, outprefix)
    #pie_read_5p_ends(read_5p_ends, genome_dict, outprefix)
if __name__ == '__main__':
    main()