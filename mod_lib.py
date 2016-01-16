from collections import defaultdict
import matplotlib.pyplot as plt
import re
import scipy.stats
import subprocess
import os
import cPickle
import mod_utils
import numpy as np
from collections import Counter

#TODO: consider adding options for ignoring nucleotides that are modified in vivo
#TODO: add methods to take in shapemapper processed data with errors, and use them for comparing libraries.
class ModLib:
    def __init__(self, experiment, experiment_settings, lib_settings):
        """
        Constructor for Library class
        """
        self.experiment = experiment
        self.experiment_settings = experiment_settings
        self.lib_settings = lib_settings
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir
        self.rRNA_mutation_data = {}    #maps rRNA names to rRNA_mutations objects, which are containers for nucleotide
                                        # objects for that rRNA
        self.parse_shapemapper_output_files()


    def parse_shapemapper_output_files(self):
        shapemapper_output_dir = os.path.join(os.path.dirname(self.experiment_settings.get_shapemapper_config_file()),
                                                  'output', 'counted_mutations_columns')
        sample_name = self.lib_settings.sample_name
        for rRNA_name in self.experiment_settings.rRNA_seqs:
            shapemapper_output_file = os.path.join(shapemapper_output_dir, sample_name+'_'+rRNA_name+'.csv')
            assert mod_utils.file_exists(shapemapper_output_file)
            self.rRNA_mutation_data[rRNA_name] = rRNA_mutations(self, self.lib_settings, self.experiment_settings,
                                                                shapemapper_output_file)

    def count_mutation_rates_by_nucleotide(self, subtract_background = False):
        """
        counts, over all RNAs, the total number of mutation rates at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.
        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        total_counts = defaultdict(int)

        for rRNA_name in self.rRNA_mutation_data:
            rRNA_counts = self.rRNA_mutation_data[rRNA_name].count_mutation_rates_by_nucleotide(subtract_background=subtract_background)
            for nucleotide_type in rRNA_counts:
                total_counts[nucleotide_type] += rRNA_counts[nucleotide_type]
        return total_counts

    def list_mutation_rates(self, subtract_background = False, nucleotides_to_count = 'ATCG'):
        all_mutation_rates = []
        for rRNA_name in self.rRNA_mutation_data:
            all_mutation_rates.extend(self.rRNA_mutation_data[rRNA_name].
                                      list_mutation_rates(subtract_background = subtract_background,
                                                          nucleotides_to_count = nucleotides_to_count))
        return all_mutation_rates

    def get_normalizing_lib(self):
        """
        #returns the library that is the normalization for this one (no-modification control)
        """
        if self.lib_settings.sample_name in self.experiment_settings.get_property('experimentals'):
            lib_index = self.experiment_settings.get_property('experimentals').index(self.lib_settings.sample_name)
            normalizing_lib_name = self.experiment_settings.get_property('no_mod_controls')[lib_index]
            return self.experiment.get_lib_from_name(normalizing_lib_name)
        else:
            return None

    def get_mutation_count_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].total_mutation_counts

    def get_coverage_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].sequencing_depth

    def get_mutation_rate_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].mutation_rate

    def pickle_mutation_rates(self, output_name, subtract_background = False):
        """
        stores mutation rates as a simple pickle, of {rRNA_name:{position:mutation rate}}
        :param subtract_background:
        :return:
        """
        output_dict = {}
        for rRNA in self.rRNA_mutation_data:
            output_dict[rRNA] = {}
            for position in self.rRNA_mutation_data[rRNA].nucleotides:
                nucleotide =  self.rRNA_mutation_data[rRNA].nucleotides[position]
                if subtract_background:
                    output_dict[rRNA][position] = max((nucleotide.mutation_rate - self.get_normalizing_lib().
                                                get_mutation_rate_at_position(rRNA, nucleotide.position)), 0.)
                else:
                    output_dict[rRNA][position] = nucleotide.mutation_rate
        mod_utils.makePickle(output_dict, output_name)

    def write_mutation_rates_to_wig(self, output_prefix, subtract_background = False):
        """
        write out mutation rates to a wig file that can be opened with a program like IGV or mochiview,
        given the corresponding rRNA fasta as a genome, of course
        :param output_prefix:
        :param subtract_background:
        :return:
        """
        wig = gzip.open(output_prefix+'.wig.gz', 'w')
        wig.write('track type=wiggle_0 name=%s\n' % (self.lib_settings.sample_name))
        for rRNA_name in self.rRNA_mutation_data:
            wig.write('variableStep chrom=%s\n' % (chr))
            for position in sorted(self.rRNA_mutation_data[rRNA_name].nucleotides.keys()):
                if subtract_background:
                    wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                            nucleotides[position].mutation_rate))
                else:
                    wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                            nucleotides[position].mutation_rate))
        wig.close()





class rRNA_mutations:
    def __init__(self, lib, lib_settings, experiment_settings, mutation_filename):
        self.lib = lib
        self.lib_settings = lib_settings
        self.experiment_settings = experiment_settings
        self.nucleotides = {}
        self.parse_mutations_columns(mutation_filename)


    def parse_mutations_columns(self, filename):
        f= open(filename, 'rU')
        lines =  f.readlines()
        sample_name = lines[0].split(',')[0]
        assert sample_name == self.lib_settings.sample_name
        self.rRNA_name = lines[1].split(',')[0]
        self.sequence = self.experiment_settings.rRNA_seqs[self.rRNA_name]
        headers = lines[2].strip().split(',')
        for line in lines[3:]:
            if line.strip().strip(',') != '':
                nucleotide_data = Nucleotide(self, headers, line)
                self.nucleotides[nucleotide_data.position] = nucleotide_data
        f.close()

    def count_mutation_rates_by_nucleotide(self, subtract_background = False):
        """
        counts, over this RNA, the total number of mutations at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.

        NOTE that this will set any background-subtracted rate of less than zero to zero

        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        counts = defaultdict(int)
        for nucleotide in self.nucleotides.values():
            if subtract_background:
                counts[nucleotide.identity] += max((nucleotide.mutation_rate - self.lib.get_normalizing_lib().
                                                get_mutation_rate_at_position(self.rRNA_name, nucleotide.position)), 0.)
            else:
                counts[nucleotide.identity] += nucleotide.mutation_rate
        return counts

    def list_mutation_rates(self, subtract_background=False, nucleotides_to_count='ATCG'):
        """
        #note that these values may be less than zero when background is subtracted
        :param subtract_background:
        :return:
        """
        rates = []
        for nucleotide in self.nucleotides.values():
            if nucleotide.identity in nucleotides_to_count:
                if subtract_background:

                    rates.append((nucleotide.mutation_rate - self.lib.get_normalizing_lib().
                                                    get_mutation_rate_at_position(self.rRNA_name, nucleotide.position)))
                else:
                    rates.append(nucleotide.mutation_rate)
        return rates

class Nucleotide:
    def __init__(self, rRNA, headers, mutation_data_line):
        self.rRNA = rRNA
        self.mutations_by_type = {} #will map each type of mutation to the number of such mutations detected
        self.parse_mutation_data_line(headers, mutation_data_line)

    def parse_mutation_data_line(self, headers, mutation_data_line):
        ll = mutation_data_line.strip().split(',')
        self.position = int(ll[0])
        self.identity = ll[1]
        assert self.rRNA.sequence[self.position-1] == self.identity #the rRNA is 1-indexed, but python strings 0-indexed
        self.total_mutation_counts = sum([float(ll[i]) for i in range(2, 18)])
        self.sequencing_depth = float(ll[19])
        try:
            self.mutation_rate = self.total_mutation_counts/self.sequencing_depth
        except:
            self.mutation_rate = 0
        for i in range(2, 18):
            self.mutations_by_type[headers[i]] = float(ll[i])
        self.back_sub_mutation_rate = self.mutation_rate - \
                                      self.rRNA.lib.get_mutation_rate_at_position(self.rRNA.rRNA_name, self.position)



