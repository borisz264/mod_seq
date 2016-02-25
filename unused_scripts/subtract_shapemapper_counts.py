__author__ = 'boris'
"""
takes:
    all_counts - pickled dict of mutation counts,  all_counts[rRNA_name][sample_name] = counts_table
    all_Depths - pickled dict of coverage counts , all_depths[rRNA_name][sample_name] = depth_table
    min_mutations: if a position has less coverage than less mutations than this across both datasets, the subtracted total will be set to zero, and the ratio to 1
    output_prefix
    comparisons: any number of sample,control name pairs

makes for each comparison_pair:
    a dict of subtraction-normed mutation rates
    a dict of variance for subtracted rates
    a dict of division-normed mutation rates
    a dict of variance for division-normed rates


"""

import sys, math, mod_utils
from collections import defaultdict


def write_out_counts(subtracted_rates, subtraction_errors, divided_rates, division_errors, rRNA, output_filename):
    f = open(output_filename, 'w')
    f.write('position\tsubtracted\tsub error\tdivided\tdiv error\n')
    for position in sorted(subtracted_rates[rRNA]):

        line = '%d\t%f\t%f\t%f\t%f\n' % (position, subtracted_rates[rRNA][position], subtraction_errors[rRNA][position], divided_rates[rRNA][position],division_errors[rRNA][position])
        f.write(line)

    f.close()

def standard_error(mutations, coverage):
    return math.sqrt(float(mutations))/float(coverage)

def subtraction_norm(all_counts, all_depths, min_mutations, comparison):
    sample = comparison[0]
    control = comparison[1]

    normalized_rates = defaultdict(dict)
    normalized_errors = defaultdict(dict)
    for rRNA_name in all_counts:
        for position in all_counts[rRNA_name][sample]:
            if  all_counts[rRNA_name][sample][position] >= min_mutations and all_counts[rRNA_name][control][position] >= min_mutations:
                sample_ratio = float(all_counts[rRNA_name][sample][position])/float(all_depths[rRNA_name][sample][position])
                sample_error = standard_error(all_counts[rRNA_name][sample][position], all_depths[rRNA_name][sample][position])
                control_ratio = float(all_counts[rRNA_name][control][position])/float(all_depths[rRNA_name][control][position])
                control_error = standard_error(all_counts[rRNA_name][control][position], all_depths[rRNA_name][control][position])

                normalized_rates[rRNA_name][position] = sample_ratio-control_ratio
                normalized_errors[rRNA_name][position] = math.sqrt((sample_error**2)+(control_error**2))
            else:
                normalized_rates[rRNA_name][position] = 0
                normalized_errors[rRNA_name][position] = 0
    return normalized_rates, normalized_errors

def division_norm(all_counts, all_depths, min_mutations, comparison):
    sample = comparison[0]
    control = comparison[1]

    normalized_rates = defaultdict(dict)
    normalized_errors = defaultdict(dict)
    for rRNA_name in all_counts:
        for position in all_counts[rRNA_name][sample]:
            if  all_counts[rRNA_name][sample][position] >= min_mutations and all_counts[rRNA_name][control][position] >= min_mutations:
                sample_ratio = float(all_counts[rRNA_name][sample][position])/float(all_depths[rRNA_name][sample][position])
                sample_error = standard_error(all_counts[rRNA_name][sample][position], all_depths[rRNA_name][sample][position])
                control_ratio = float(all_counts[rRNA_name][control][position])/float(all_depths[rRNA_name][control][position])
                control_error = standard_error(all_counts[rRNA_name][control][position], all_depths[rRNA_name][control][position])

                normalized_rates[rRNA_name][position] = math.log(sample_ratio/control_ratio, 2)
                normalized_errors[rRNA_name][position] = ((sample_ratio/control_ratio)*math.sqrt((sample_error/sample_ratio)**2+(control_error/control_ratio)**2))/((sample_ratio/control_ratio)*math.log(2))
                #this is a standard error propogation formula
            else:
                normalized_rates[rRNA_name][position] = 0
                normalized_errors[rRNA_name][position] = 0
    return normalized_rates, normalized_errors

def main():
    all_counts_file, all_depths_file, min_mutations, output_prefix = sys.argv[1:5]
    min_mutations = int(min_mutations)
    all_counts = mod_utils.unPickle(all_counts_file)
    all_depths = mod_utils.unPickle(all_depths_file)

    comparisons = (pair.split(',') for pair in sys.argv[5:])

    for comparison in comparisons:
        subtracted_rates, subtraction_errors = subtraction_norm(all_counts, all_depths, min_mutations, comparison)
        divided_rates, division_errors = division_norm(all_counts, all_depths, min_mutations, comparison)
        mod_utils.makePickle(subtracted_rates, '%s_%s_%s_sub_norm.pkl' % (output_prefix, comparison[0], comparison[1]))
        mod_utils.makePickle(subtraction_errors, '%s_%s_%s_sub_err.pkl' % (output_prefix, comparison[0], comparison[1]))
        mod_utils.makePickle(divided_rates, '%s_%s_%s_div_norm.pkl' % (output_prefix, comparison[0], comparison[1]))
        mod_utils.makePickle(division_errors, '%s_%s_%s_div_err.pkl' % (output_prefix, comparison[0], comparison[1]))
        for rRNA in subtracted_rates:
            write_out_counts(subtracted_rates, subtraction_errors, divided_rates, division_errors, rRNA, '%s_%s_%s_%s.txt' % (output_prefix, comparison[0], comparison[1], rRNA))
main()