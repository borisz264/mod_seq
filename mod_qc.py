from collections import defaultdict
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import scipy.stats as stats
import subprocess
import os
import cPickle
import mod_utils
import numpy as np
import itertools
import pysam
import math
import gzip
import bzUtils

class TPS_qc:
    def __init__(self, tpse, experiment_settings, threads):
        """
        Constructor for Library class
        """
        self.threads = threads
        self.tpse = tpse
        self.experiment_settings = experiment_settings
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir
        mod_utils.make_dir(self.tpse.rdir_path('QC'))

    def identify_contaminating_sequences(self):
        for lib_settings in self.experiment_settings.iter_lib_settings():
            self.map_for_contaminating_sequences_one_lib(lib_settings)
        for lib_settings in self.experiment_settings.iter_lib_settings():
            self.write_mapping_summary(lib_settings.get_rRNA_mapping_stats(), lib_settings.get_pool_mapping_stats(), lib_settings.get_genome_mapping_stats(), lib_settings.get_overall_contamination_summary())


    def map_for_contaminating_sequences_one_lib(self, lib_settings):
        #first, take unmapped sequences and map them to yeast rRNA, counting mapping stats
        if not mod_utils.file_exists(lib_settings.get_rRNA_unmapped_reads()):
            subprocess.Popen('bowtie2 -f -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 -x %s -p %d -U %s --un-gz %s 2>>%s | samtools view -bS - > %s 2>>%s ' % (self.experiment_settings.get_rRNA_bowtie_index(), self.threads,
                                                                                                       lib_settings.get_unmappable_reads(), lib_settings.get_rRNA_unmapped_reads(), lib_settings.get_rRNA_mapping_stats(),
                                                                                                       lib_settings.get_rRNA_mapped_reads(), lib_settings.get_log(),
                                                                                                       ), shell=True).wait()
        if not mod_utils.file_exists(lib_settings.get_genome_unmapped_reads()):
            #take still unmapped sequences and map them to the rest of the yeast genome, counting mapping stats
            subprocess.Popen('bowtie2 -f -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 -x %s -p %d -U %s --un-gz %s 2>>%s | samtools view -bS - > %s 2>>%s ' % (self.experiment_settings.get_genome_bowtie_index(), self.threads,
                                                                                               lib_settings.get_rRNA_unmapped_reads(), lib_settings.get_genome_unmapped_reads(), lib_settings.get_genome_mapping_stats(),
                                                                                               lib_settings.get_genome_mapped_reads(), lib_settings.get_log(),
                                                                                               ), shell=True).wait()

    def write_mapping_summary(self, rRNA_file, pool_file, genome_file, output_file):
        pool_stats = self.parse_mapping_stats(pool_file)
        rRNA_stats = self.parse_mapping_stats(rRNA_file)
        genome_stats = self.parse_mapping_stats(genome_file)

        f = open(output_file, 'w')
        f.write('\tunique_pool\tmultiple_pool\tunique_rRNA\tmultiple_rRNA\tunique_genome\tmultiple_genome\n')
        f.write('total\t%d\t%d\t%d\t%d\t%d\t%d\n' % (pool_stats[2], pool_stats[3], rRNA_stats[2], rRNA_stats[3],
                                                     genome_stats[2], genome_stats[3]))
        f.close()


    def parse_mapping_stats(self, alignment_summary_file):
        '''
        example alignment summary:
        8333978 reads; of these:
          8333978 (100.00%) were unpaired; of these:
            7905371 (94.86%) aligned 0 times
            276859 (3.32%) aligned exactly 1 time
            151748 (1.82%) aligned >1 times
        5.14% overall alignment rate
        '''
        f = open(alignment_summary_file)
        lines = f.readlines()
        total_reads = int(lines[0].strip().split()[0])
        unaligned_reads = int(lines[2].strip().split()[0])
        uniquely_aligned_reads = int(lines[3].strip().split()[0])
        multiply_aligned_reads = int(lines[4].strip().split()[0])
        overall_alignment_percent = float(lines[5].strip().split()[0][:-1])
        f.close()
        return total_reads, unaligned_reads, uniquely_aligned_reads, multiply_aligned_reads, overall_alignment_percent

    def print_library_count_concordances(self):
        out_name =  os.path.join(self.experiment_settings.get_rdir(), 'QC',
          'count_concordances.txt')
        f = open(out_name, 'w')
        header = 'sample1\tsample2\tpearson r\t pearson p\t spearman r\t spearman p\n'
        f.write(header)
        for libi, libj in itertools.combinations(self.tpse.libs, 2):
            pearsonR, spearmanR, pearsonP, spearmanP = self.get_library_count_correlation(libi, libj)
            line = '%s\t%s\t%f\t%f\t%f\t%f\n' % (libi.get_sample_name(), libj.get_sample_name(),
                                                         pearsonR, pearsonP, spearmanR, spearmanP)
            f.write(line)
        f.close()


    def plot_average_read_positions(self):
        for lib in self.tpse.libs:
            self.plot_average_read_positions_one_lib(lib)

    def plot_average_read_positions_one_lib(self, lib, min_x = 0, max_x = 150):
        positions = np.array(range(min_x, max_x+1))
        averages = [np.average([pool_sequence_mapping.fraction_at_position(position) for pool_sequence_mapping in lib.pool_sequence_mappings.values() if pool_sequence_mapping.total_passing_reads>0]) for position in positions]

        fig = plt.figure(figsize=(8,8))
        plot = fig.add_subplot(111)
        plot.bar(positions , averages,color=bzUtils.rainbow[0], lw=0)
        plot.set_xticks(positions[::10]+0.5)
        plot.set_xticklabels(positions[::10])
        plot.set_xlabel("position of read 5' end from RNA end")
        plot.set_ylabel("average read fraction")
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          '%(sample_name)s.read_positions.pdf' % {'sample_name': lib.get_sample_name ()})
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()


    def plot_count_distributions(self):
        num_libs = len(self.tpse.libs)
        fig = plt.figure(figsize=(16,16))
        plot_index = 1
        cutoff = 100
        hbins = np.arange(0, 400, 10)
        hbins = np.append(hbins, 10000000)
        for lib in self.tpse.libs:
            plot = fig.add_subplot(math.sqrt(bzUtils.next_square_number(num_libs)), math.sqrt(bzUtils.next_square_number(num_libs)), plot_index)
            sample_name = lib.lib_settings.sample_name
            dist = self.get_library_count_distribution(lib)
            plot.hist(dist, bins = hbins, color=bzUtils.skyBlue, histtype='stepfilled', edgecolor = None, lw = 0)
            plot.set_xlabel("# reads", fontsize = 10)
            plot.set_ylabel("# genes (%d have >= %d reads)" % (bzUtils.number_passing_cutoff(dist, cutoff), cutoff), fontsize = 10)
            plot.set_xlim(0, 400)
            #plot.set_ylim(0,1)
            plot.axvline(cutoff, ls = 'dashed')
            plot.set_title(sample_name, fontsize = 8)
            plot_index += 1
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15, wspace=0.4, hspace=0.6)
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          'count_distributions.pdf')
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()

        """
        def plot_insert_size_distributions(self):
            #plot distribution of insert sizes from cutadapt output
            TODO - need to parse log file to get this info
            num_libs = len(self.tpse.libs)
            fig = plt.figure(figsize=(16,16))
            plot_index = 1
            cutoff = 100
            hbins = np.arange(0, 51, 1)
            for lib in self.tpse.libs:
                plot = fig.add_subplot(math.sqrt(bzUtils.next_square_number(num_libs)), math.sqrt(bzUtils.next_square_number(num_libs)), plot_index)
                sample_name = lib.lib_settings.sample_name
                dist = self.get_insert_sizes(lib)
                plot.hist(dist, bins = hbins, color=bzUtils.skyBlue, histtype='stepfilled', edgecolor = None, lw = 0)
                plot.set_xlabel("insert size", fontsize = 10)
                plot.set_ylabel("fraction of reads" % (bzUtils.number_passing_cutoff(dist, cutoff), cutoff), fontsize = 10)
                plot.set_xlim(0, 400)
                #plot.set_ylim(0,1)
                plot.axvline(cutoff, ls = 'dashed')
                plot.set_title(sample_name, fontsize = 8)
                plot_index += 1
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15, wspace=0.4, hspace=0.6)
            out_name =  os.path.join(
              self.experiment_settings.get_rdir(),
              'QC',
              'insetrt_size_distributions.pdf')
            plt.savefig(out_name, transparent='True', format='pdf')
            plt.clf()
        """