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

    def count_mutations_by_nucleotide(self, subtract_background = False):
        """
        counts, over all RNAs, the total number of mutations at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.
        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        total_counts = defaultdict(int)

        for rRNA_name in self.rRNA_mutation_data:
            rRNA_counts = self.rRNA_mutation_data[rRNA_name].count_mutations_by_nucleotide(subtract_background=subtract_background)
            for nucleotide_type in rRNA_counts:
                total_counts[nucleotide_type] += rRNA_counts[nucleotide_type]
        return total_counts

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

    def count_mutations_by_nucleotide(self, subtract_background = False):
        """
        counts, over this RNA, the total number of mutations at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.
        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        counts = defaultdict(int)
        for nucleotide in self.nucleotides.values():
            if subtract_background:
                counts[nucleotide.identity] += (nucleotide.total_mutation_counts - self.lib.get_normalizing_lib().
                                                get_mutation_count_at_position(self.rRNA_name, nucleotide.position))
            else:
                counts[nucleotide.identity] += nucleotide.total_mutation_counts
        return counts

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
        self.total_mutation_counts = sum([int(ll[i]) for i in range(2, 18)])
        self.sequencing_depth = int(ll[19])
        for i in range(2, 18):
            self.mutations_by_type[headers[i]] = int(ll[i])



