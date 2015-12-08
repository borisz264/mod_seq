__author__ = 'boris'

"""
Takes a directory of mutation counts from the smapemapper pipeline, then for each rRNA output spreadsheets of mutation
rates and coverage depth for each dataset. (2 spreadsheets per rRNA)
"""

import sys, os, mod_utils
from collections import defaultdict

def parse_mutations_columns(filename):

    counts_table = {}
    depth_table = {}

    f= open(filename, 'rU')
    lines =  f.readlines()
    sample_name = lines[0].split(',')[0]
    rRNA_name = lines[1].split(',')[0]
    for line in lines[3:]:
        ll = line.strip().split(',')
        try:
            nt_pos = int(ll[0])
            counts = sum([int(ll[i]) for i in range(2, 18)])
            depth = int(ll[19])
            counts_table[nt_pos] = counts
            depth_table[nt_pos] = depth
        except:
            continue
    f.close()

    return sample_name, rRNA_name, counts_table, depth_table

def normalize(counts_table, depth_table):
    rate_table = {}
    for position in counts_table:
        rate_table[position] = float(counts_table[position])/float(depth_table[position])


def write_out_counts(counts_dict, output_filename):
    f = open(output_filename, 'w')
    samples = sorted(counts_dict.keys())
    header = 'nt\t'+'\t'.join(samples)+'\n'
    f.write(header)
    for position in sorted(counts_dict[samples[0]]):
        values = '\t'.join(str(counts_dict[sample][position]) for sample in samples)
        line = '%d\t%s\n' % (position, values)
        f.write(line)

    f.close()

def main():
    infolder, outprefix = sys.argv[1:3]
    all_counts = defaultdict(dict)
    all_depths = defaultdict(dict)
    all_rates = defaultdict(dict)
    for filename in os.listdir(infolder):
        if filename.endswith('.csv'):
            sample_name, rRNA_name, counts_table, depth_table = parse_mutations_columns(os.path.join(infolder, filename))
            all_counts[rRNA_name][sample_name] = counts_table
            all_depths[rRNA_name][sample_name] = depth_table
            #all_rates[sample_name][rRNA_name] = normalize(counts_table, depth_table)
    for rRNA_name in all_counts:
        write_out_counts(all_counts[rRNA_name], '%s_%s_mutations.txt' % (outprefix, rRNA_name))
        write_out_counts(all_depths[rRNA_name], '%s_%s_coverage.txt' % (outprefix, rRNA_name))

    mod_utils.makePickle(all_counts, outprefix+'_mutations.pkl')
    mod_utils.makePickle(all_depths, outprefix+'_depths.pkl')
    #mod_utils.makePickle(all_depths, outprefix+'_mutation_rates.pkl')
main()