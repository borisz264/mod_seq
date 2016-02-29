__author__ = 'boris'
"""
inputs:
    outprefix
    bundle 1
    bundle 2
    bundle 3
    bundle 4
    bundle 5 - the 5 pdb files from the 4v88 bundle
    reactivity_values - a pickled dict of [chromosome][position] = reactivity_value or change, such as from compare_samples.py

outputs:
    the 5 PDB files in the bundle, with b factors replaced with reactivity scores for the corresponding rRNA residues.
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
    allChrs=set(mutation_dict.keys()).union(set(mutation_dict['-'].keys()))
    for chr in allChrs:
        #minPos = min([min(weightedReads['+'][chr].keys()), min(weightedReads['-'][chr].keys())])
        #maxPos = max([max(weightedReads['+'][chr].keys()), max(weightedReads['-'][chr].keys())])
        plusWig.write('variableStep chrom=%s\n' % (chr))
        #minusWig.write('variableStep chrom=%s\n' % (chr))
        if chr in mutation_dict:
            for i in sorted(mutation_dict[chr].keys()):
                plusWig.write('%d\t%f\n' % (i, mutation_dict[chr][i]))
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
    for chromosome in mutation_dict:
        normed_dict[chromosome] = {}
        #print mutation_dict[strand][chromosome].values()
        all_values += mutation_dict[chromosome].values()
    if winsorize_data:
        winsorize(all_values, limits = (winsorization_limits[0], 1-winsorization_limits[1]), inplace = True)

    max_value = float(max(all_values))
    min_value = float(min(all_values))
    norm_value = max(max_value, -1*min_value)
    print max_value
    for chromosome in mutation_dict:
        for position in mutation_dict[chromosome]:
            val = mutation_dict[chromosome][position]
            if val < min(all_values):
                val = min(all_values)
            if val > max(all_values):
                val = max(all_values)
            normed_dict[chromosome][position] = val/norm_value
    return normed_dict

def split_by_n(line, n=6):
    """
    trims endline, and returns set of chunks of length 6, each of which has whitespace stripped
    :param line:
    :param n:
    :return:
    """
    return [line[i:i+n].strip() for i in range(0, len(line), 6)]

#rRNA_assignments = {3:{'d':"S.c.18S_rRNA"}, 1:{'A':"S.c.18S_rRNA"}, 2:{'A':"S.c.25S__rRNA", 'B':"S.c.5S___rRNA", 'C':"S.c.5.8S_rRNA"} , 4:{'K':"S.c.25S__rRNA", 'L':"S.c.5S___rRNA", 'M':"S.c.5.8S_rRNA"}}
#for 4V88:
rRNA_assignments = {3:{'d':"S.c.18S_rRNA"}, 1:{'A':"S.c.18S_rRNA"}, 2:{'A':"S.c.25S__rRNA", 'B':"S.c.5S___rRNA", 'C':"S.c.5.8S_rRNA"} , 4:{'I':"S.c.25S__rRNA", 'L':"S.c.5S___rRNA", 'M':"S.c.5.8S_rRNA"}}
#for 4V6I:
#rRNA_assignments = {2:{'S':"S.c.18S_rRNA"}, 3:{'A':"S.c.25S__rRNA", 'C':"S.c.5S___rRNA", 'B':"S.c.5.8S_rRNA"}}
#for 4U3U:
#rRNA_assignments = {1:{'A':"S.c.18S_rRNA"}, 2:{'A':"S.c.25S__rRNA", 'C':"S.c.5S___rRNA", 'B':"S.c.5.8S_rRNA"}}

def main():
    outprefix, bundle1, bundle2, bundle3, bundle4, bundle5, datafile_name  = sys.argv[1:8]

    bundles = [bundle1, bundle2, bundle3, bundle4, bundle5]
    reactivities = mod_utils.unPickle(datafile_name)

    scaling_factor = 1000. #all values will be multiplied byt his factor. This is necessary because the values produced are often very small, and the B factors in the PDB file are limited to 6 characters

    for i in range(1,6):
        infile = open(bundles[i-1])
        outfile = open(outprefix+'_bundle'+str(i)+'.pdb' ,'w')

        for line in infile:
            if line.startswith('ATOM'):
                chain = line[21]
                resi = int(line[22:28].strip())
                #print i, chain, resi
                #print rRNA_assignments[i][chain]
                #print reactivities[rRNA_assignments[i][chain]]
                if i in rRNA_assignments and chain in rRNA_assignments[i] and resi in reactivities[rRNA_assignments[i][chain]]:
                    insert = '%6.3f' % (scaling_factor*reactivities[rRNA_assignments[i][chain]][resi])
                    new_line = '%s%s%s' % (line[:60], insert[:6], line[66:])
                    assert len(line) == len(new_line)
                else:
                    new_line = '%s%6.4f%s' % (line[:60], 0.0, line[66:])
                    assert len(line) == len(new_line)

            elif line.startswith("ANISOU") or line.startswith("HETATM"):
                new_line = '' #remove the anisotropic b factors, I don't need them
            else:
                new_line = line

            outfile.write(new_line)
        infile.close()
        outfile.close()





main()