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
import bzUtils

class ModLib:
    def __init__(self, experiment_settings, lib_settings):
        """
        Constructor for Library class
        """
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
            self.rRNA_mutation_data[rRNA_name] = rRNA_mutations(self.lib_settings, self.experiment_settings,
                                                                shapemapper_output_file)

class rRNA_mutations:
    def __init__(self, lib_settings, experiment_settings, mutation_filename):
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

