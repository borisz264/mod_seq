import operator
import aColors

__author__ = 'boris'
"""
Intended for processing of toeprint-seq data from defined RNA pools
Based on Alex Robertson's original RBNS pipeline, available on github
"""
import sys
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import argparse
import itertools
import collections
from collections import defaultdict
import gzip
import subprocess
import numpy
import scipy.stats as stats

import bzUtils
import tps_settings
import tps_utils
import tps_lib
import tps_qc
import stacked_bar_kmers


class TPSe:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.collapse_identical_reads()
        self.remove_adaptor()
        self.remove_primer()
        self.trim_reads()
        self.trim_reference_pool_fasta()
        self.build_bowtie_index()
        self.map_reads()
        self.initialize_libs()

    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
        tps_utils.make_dir(self.rdir_path('sequence_counts'))
        self.libs = []
        map(lambda lib_settings: self.initialize_lib(lib_settings), self.settings.iter_lib_settings())
        self.settings.write_to_log('initializing libraries, counting reads, done')


    def initialize_lib(self, lib_settings):
        lib = tps_lib.TPS_Lib(self.settings, lib_settings)
        self.libs.append(lib)

    def needs_calculation(self, lib_settings, count_type, k):
        if self.settings.get_force_recount(count_type):
            return True
        return not lib_settings.counts_exist(count_type, k)


    def make_tables(self):
        tps_utils.make_dir(self.rdir_path('tables'))
        self.make_counts_table()

    def make_plots(self):
        tps_utils.make_dir(self.rdir_path('plots'))
        self.plot_AUG_reads()
        self.plot_AUG_reads(unique_only = True,)
        self.plot_last_AUG_reads()
        self.plot_last_AUG_reads(unique_only = True,)
        self.plot_AUG_reads(which_AUG = 2, unique_only = True)
        self.plot_AUG_reads(which_AUG = 2)


    def make_table_header(self, of):
        """
        takes a file handle and writes a good header for it such that
        each lane is a column.
        """
        of.write('#')
        for lib in self.libs:
            of.write('\t' + lib.get_barcode())
        of.write('\n[%s]' % self.settings.get_property('protein_name'))
        for lib in self.libs:
            of.write('\t%s' % lib.get_conc())
        of.write('\nwashes')
        for lib in self.libs:
            of.write('\t%i' % lib.get_washes())
        of.write('\nT (C)')
        for lib in self.libs:
            of.write('\t%s' % lib.get_temperature())
        of.write('\n')

    def collapse_identical_reads(self):
        """
        collapses all identical reads using FASTX toolkit
        :return:
        """
        self.settings.write_to_log('collapsing reads')
        if not self.settings.get_property('force_recollapse'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.collapsed_reads_exist():
                    break
            else:
                return
        tps_utils.make_dir(self.rdir_path('collapsed_reads'))
        if self.settings.get_property('collapse_identical_reads'):
            bzUtils.parmap(lambda lib_setting: self.collapse_one_fastq_file(lib_setting), self.settings.iter_lib_settings(), nprocs = self.threads)
        else:
            bzUtils.parmap(lambda lib_setting: self.fastq_to_fasta(lib_setting), self.settings.iter_lib_settings(), nprocs = self.threads)
        self.settings.write_to_log('collapsing reads complete')

    def collapse_one_fastq_file(self, lib_settings):
        lib_settings.write_to_log('collapsing_reads')
        subprocess.Popen('gunzip -c %s | fastx_collapser -v -Q33 2>>%s | gzip > %s' % (lib_settings.get_fastq_file(),
                                                                                  lib_settings.get_log(),
                                                                                  lib_settings.get_collapsed_reads()
                                                                                  ), shell=True).wait()
        lib_settings.write_to_log('collapsing_reads_done')
    def fastq_to_fasta(self, lib_settings):
        lib_settings.write_to_log('fasta_conversion')
        subprocess.Popen('gunzip -c %s | fastq_to_fasta -v -Q33 2>>%s | gzip > %s' % (lib_settings.get_fastq_file(),
                                                                                  lib_settings.get_log(),
                                                                                  lib_settings.get_collapsed_reads()
                                                                                  ), shell=True).wait()
        lib_settings.write_to_log('fasta_conversion done')
    def remove_adaptor(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.adaptorless_reads_exist():
                    break
            else:
                return

        if self.settings.get_property('trim_adaptor'):
            tps_utils.make_dir(self.rdir_path('adaptor_removed'))
            bzUtils.parmap(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting), self.settings.iter_lib_settings(), nprocs = self.threads)

    def remove_primer(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.primerless_reads_exist():
                    break
            else:
                return

        if self.settings.get_property('trim_adaptor'):
            tps_utils.make_dir(self.rdir_path('primer_removed'))
            bzUtils.parmap(lambda lib_setting: self.remove_primer_one_lib(lib_setting), self.settings.iter_lib_settings(), nprocs = self.threads)

    def remove_adaptor_one_lib(self, lib_settings):
        lib_settings.write_to_log('adaptor trimming')
        command_to_run = 'cutadapt --adapter %s --overlap 3 --minimum-length %d %s --output %s 1>>%s 2>>%s' % (self.settings.get_property('adaptor_sequence'), self.settings.get_property('min_post_adaptor_length'),
                           lib_settings.get_collapsed_reads(), lib_settings.get_adaptor_trimmed_reads(), lib_settings.get_log(),
                           lib_settings.get_log())
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')

    def remove_primer_one_lib(self, lib_settings):
        lib_settings.write_to_log('reverse primer trimming')
        command_to_run = 'cutadapt --adapter %s --overlap 3 --minimum-length %d %s --output %s 1>>%s 2>>%s' % (self.settings.get_property('primer_sequence'), self.settings.get_property('min_post_adaptor_length'),
                           lib_settings.get_adaptor_trimmed_reads(), lib_settings.get_primer_trimmed_reads(), lib_settings.get_log(),
                           lib_settings.get_log())
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('reverse primer trimming done')

    def plot_AUG_reads(self, which_AUG = 1, unique_only = False, min_x = -30, max_x = 30, read_cutoff = 100):
        #1 is for the first AUG, 2 the 2nd and so on. Only TLs with enough AUGs are counted
        assert  which_AUG > 0
        positions = numpy.array(range(min_x, max_x+1))
        mappings_passing_cutoff_in_all_libs = self.libs[0].get_mappings_with_minimum_reads(read_cutoff, names_only = True)
        for lib in self.libs[1:]:
            mappings_passing_cutoff_in_all_libs = \
                mappings_passing_cutoff_in_all_libs.intersection(lib.get_mappings_with_minimum_reads(read_cutoff,
                                                                                                     names_only = True))

        if unique_only:
            out_name =  os.path.join(
              self.settings.get_rdir(),
              'plots',
              'unique_AUG%d_density.pdf' % (which_AUG))
            mapping_names = mappings_passing_cutoff_in_all_libs.intersection(lib.get_single_TL_mappings(names_only = True))
        else:
            out_name =  os.path.join(
              self.settings.get_rdir(),
              'plots',
              'AUG%d_density.pdf' % (which_AUG))
            mapping_names = mappings_passing_cutoff_in_all_libs

        fig = plt.figure(figsize=(8,8))
        plot = fig.add_subplot(111)
        positions = range(min_x, max_x)
        color_index = 0
        genes_plotted = set()
        for lib in self.libs:
            offset_sum = defaultdict(float)
            offset_counts = defaultdict(int)
            num_genes_counted = 0
            for mapping_name in mapping_names:
                mapping = lib.pool_sequence_mappings[mapping_name]
                AUG_positions = mapping.positions_of_subsequence('ATG')
                if len(AUG_positions) >= which_AUG:
                    genes_plotted.add(mapping_name)
                    num_genes_counted += 1
                    alignment_position = AUG_positions[which_AUG-1]
                    for position in positions:
                        AUG_relative_position = alignment_position - position
                        read_fraction_at_position = mapping.fraction_at_position(AUG_relative_position)
                        if read_fraction_at_position != None:
                            offset_sum[position] += read_fraction_at_position
                            offset_counts[position] += 1
            offset_averages = {}
            for position in positions:
                #print position, offset_sum[position], float(offset_counts[position])
                offset_averages[position] = offset_sum[position]/float(offset_counts[position])
            offset_average_array = [offset_averages[position] for position in positions]
            plot.plot(positions, offset_average_array, color=bzUtils.rainbow[color_index], lw=2, label ='%s (%d)' %(lib.get_sample_name(), num_genes_counted))
            color_index += 1
        plot.axvline(16, ls= '--')
        plot.axvline(19, ls= '--')

        lg=plt.legend(loc=2,prop={'size':10}, labelspacing=0.2)
        lg.draw_frame(False)
        plot.set_xticks(positions[::3])
        plot.set_xticklabels(positions[::3])
        plot.set_xlabel("position of read 5' end from AUG %d" %(which_AUG) )
        plot.set_ylabel("average read fraction")
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()
        print genes_plotted
        for mapping_name in genes_plotted:
            self.plot_single_sequence_read_distributions(mapping_name)



    def plot_last_AUG_reads(self, unique_only = False, min_x = -30, max_x = 30, read_cutoff = 100):
        #1 is for the first AUG, 2 the 2nd and so on. Only TLs with enough AUGs are counted
        positions = numpy.array(range(min_x, max_x+1))
        mappings_passing_cutoff_in_all_libs = self.libs[0].get_mappings_with_minimum_reads(read_cutoff, names_only = True)
        for lib in self.libs[1:]:
            mappings_passing_cutoff_in_all_libs = \
                mappings_passing_cutoff_in_all_libs.intersection(lib.get_mappings_with_minimum_reads(read_cutoff,
                                                                                                     names_only = True))

        if unique_only:
            out_name =  os.path.join(
              self.settings.get_rdir(),
              'plots',
              'unique_last_AUG_density.pdf')
            mapping_names = mappings_passing_cutoff_in_all_libs.intersection(lib.get_single_TL_mappings(names_only = True))
        else:
            out_name =  os.path.join(
              self.settings.get_rdir(),
              'plots',
              'last_AUG_density.pdf')
            mapping_names = mappings_passing_cutoff_in_all_libs

        fig = plt.figure(figsize=(8,8))
        plot = fig.add_subplot(111)
        positions = range(min_x, max_x)
        color_index = 0
        for lib in self.libs:
            offset_sum = defaultdict(float)
            offset_counts = defaultdict(int)
            num_genes_counted = 0
            for mapping_name in mapping_names:
                mapping = lib.pool_sequence_mappings[mapping_name]
                AUG_positions = mapping.positions_of_subsequence('ATG')
                if len(AUG_positions) >= 1:
                    num_genes_counted += 1
                    alignment_position = AUG_positions[-1]
                    for position in positions:
                        AUG_relative_position = alignment_position - position
                        read_fraction_at_position = mapping.fraction_at_position(AUG_relative_position)
                        if read_fraction_at_position != None:
                            offset_sum[position] += read_fraction_at_position
                            offset_counts[position] += 1
            offset_averages = {}
            for position in positions:
                #print position, offset_sum[position], float(offset_counts[position])
                offset_averages[position] = offset_sum[position]/float(offset_counts[position])
            offset_average_array = [offset_averages[position] for position in positions]
            plot.plot(positions, offset_average_array, color=bzUtils.rainbow[color_index], lw=2, label ='%s (%d)' %(lib.get_sample_name(), num_genes_counted))
            color_index += 1
        plot.axvline(16, ls='--')
        plot.axvline(19, ls='--')
        lg=plt.legend(loc=2,prop={'size':10}, labelspacing=0.2)
        lg.draw_frame(False)
        plot.set_xticks(positions[::3])
        plot.set_xticklabels(positions[::3])
        plot.set_xlabel("position of read 5' end from last AUG")
        plot.set_ylabel("average read fraction")
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()

    def plot_single_sequence_read_distributions(self, sequence_name):
        fig = plt.figure(figsize=(8,8))
        plot = fig.add_subplot(111)
        colorIndex = 0
        for lib in self.libs:
            mapping = lib.pool_sequence_mappings[sequence_name]
            positions = numpy.array(range(0, len(mapping.full_sequence)))
            fractions = [mapping.fraction_at_position(position) for position in positions]
            plot.plot(positions , fractions,color=bzUtils.rainbow[colorIndex], lw=1, label = lib.lib_settings.sample_name)
            colorIndex+=1
        for AUG_pos in mapping.positions_of_subsequence('ATG'):
            plot.axvline(AUG_pos+16, ls='--')
            plot.axvline(AUG_pos+19, ls='--')

        plot.set_xticks(positions[::10])
        plot.set_xticklabels(positions[::10])
        plot.set_xlim(-1, len(mapping.full_sequence))
        plot.set_xlabel("position of read 5' end from RNA end (--expected AUG toeprints)")
        plot.set_ylabel("read fraction")
        lg=plt.legend(loc=2,prop={'size':10}, labelspacing=0.2)
        lg.draw_frame(False)
        out_name =  os.path.join(
          self.settings.get_rdir(),
          'plots',
          '%(sequence_name)s.read_positions.pdf' % {'sequence_name': sequence_name})
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()

    def trim_reads(self):
        """
        Trim reads by given amount, removing potential random barcoding sequences from 5' end
        Trimming from 3' end can also help if mapping is problematic by reducing chance for indels to prevent mapping
        :return:
        """
        self.settings.write_to_log( 'trimming reads')
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.trimmed_reads_exist():
                    break
            else:
                return
        tps_utils.make_dir(self.rdir_path('trimmed_reads'))
        bzUtils.parmap(lambda lib_setting: self.trim_one_fasta_file(lib_setting), self.settings.iter_lib_settings(), nprocs = self.threads)
        self.settings.write_to_log( 'trimming reads complete')

    def trim_one_fasta_file(self, lib_settings):
        lib_settings.write_to_log('trimming_reads')
        first_base_to_keep = self.settings.get_property('first_base_to_keep') #the trimmer is 1-indexed. 1 means keep every base
        last_base_to_keep = self.settings.get_property('last_base_to_keep')
        if self.settings.get_property('trim_adaptor'):
            subprocess.Popen('gunzip -c %s | fastx_trimmer -f %d -l %d -z -o %s >>%s 2>>%s' % (lib_settings.get_primer_trimmed_reads(),
                                                                                      first_base_to_keep, last_base_to_keep,
                                                                                      lib_settings.get_trimmed_reads(),
                                                                                      lib_settings.get_log(),
                                                                                      lib_settings.get_log()), shell=True).wait()
        else:
            subprocess.Popen('gunzip -c %s | fastx_trimmer -f %d -l %d -z -o %s >>%s 2>>%s' % (lib_settings.get_collapsed_reads(),
                                                                                      first_base_to_keep, last_base_to_keep,
                                                                                      lib_settings.get_trimmed_reads(),
                                                                                      lib_settings.get_log(),
                                                                                      lib_settings.get_log()), shell=True).wait()
        lib_settings.write_to_log('trimming_reads done')
    def get_barcode_match(self, barcode, barcodes):
        """
        takes a barcode and returns the one it matches (hamming <= 1)
        else
        empty string
        """
        if barcode in barcodes:
            return barcode
        for barcode_j in barcodes:
            if tps_utils.hamming_N(barcode, barcode_j) <= self.settings.get_property('mismatches_allowed_in_barcode'):
                return barcode_j
        return ''

    def build_bowtie_index(self):
        """
        builds a bowtie 2 index from the input fasta file
        recommend including barcode+PCR sequences just in case of some no-insert amplicons
        """
        self.settings.write_to_log('building bowtie index')
        if self.settings.get_property('force_index_rebuild') or not self.settings.bowtie_index_exists():
            tps_utils.make_dir(self.rdir_path('bowtie_indices'))
            index = self.settings.get_bowtie_index()
            subprocess.Popen('bowtie2-build -f --offrate 0 %s %s 1>>%s 2>>%s' % (self.settings.get_trimmed_pool_fasta(),
                                                                      self.settings.get_bowtie_index(), self.settings.get_log()+'.bwt',
                                                                      self.settings.get_log()+'.bwt'), shell=True).wait()
        self.settings.write_to_log('building bowtie index complete')

    def trim_reference_pool_fasta(self):
        '''
        Trims the reference sequences to the length of the trimmed reads + a buffer
        '''
        trim_5p = self.settings.get_property('pool_5trim') #nucleotides to cut from 5' end
        trim_3p = self.settings.get_property('pool_3trim') #nucleotides to cut from 3' end
        f = open(self.settings.get_property('pool_fasta'))
        g = open(self.settings.get_trimmed_pool_fasta(), 'w')
        for line in f:
            if not line.strip() == '' and not line.startswith('#'):#ignore empty lines and commented out lines
                if line.startswith('>'):#> marks the start of a new sequence
                    g.write(line)
                else:
                    g.write(self.settings.get_property('pool_prepend')+line.strip()[trim_5p:len(line.strip())-trim_3p]+self.settings.get_property('pool_append')+'\n')
        f.close()
        g.close()

    def map_reads(self):
        """
        map all reads using bowtie
        :return:
        """
        self.settings.write_to_log('mapping reads')
        if not self.settings.get_property('force_remapping'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.mapped_reads_exist():
                    break
            else:
                return
        tps_utils.make_dir(self.rdir_path('mapped_reads'))
        tps_utils.make_dir(self.rdir_path('mapping_stats'))
        tps_utils.make_dir(self.rdir_path('unmapped_reads'))

        bzUtils.parmap(lambda lib_setting: self.map_one_library(lib_setting), self.settings.iter_lib_settings(),
                       nprocs = self.threads)
        self.settings.write_to_log( 'finished mapping reads')

    def map_one_library(self, lib_settings):
        lib_settings.write_to_log('mapping_reads')
        subprocess.Popen('bowtie2 -f -D 20 -R 3 -N 1 -L 15 --norc -i S,1,0.50 -x %s -p %d -U %s --un-gz %s -S %s 1>> %s 2>>%s' % (self.settings.get_bowtie_index(), self.threads,
                                                                                                   lib_settings.get_trimmed_reads(), lib_settings.get_unmappable_reads(), lib_settings.get_mapped_reads_sam(),
                                                                                                                      lib_settings.get_log(), lib_settings.get_pool_mapping_stats()), shell=True).wait()
        #subprocess.Popen('samtools view -b -h -o %s %s 1>> %s 2>> %s' % (lib_settings.get_mapped_reads(), lib_settings.get_mapped_reads_sam(), lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()
        #also, sort bam file, and make an index

        #samtools view -uS myfile.sam | samtools sort - myfile.sorted
        subprocess.Popen('samtools view -uS %s | samtools sort - %s.temp_sorted 1>>%s 2>>%s' % (lib_settings.get_mapped_reads_sam(), lib_settings.get_mapped_reads_sam(),
                                                                          lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()


        #subprocess.Popen('samtools sort %s %s.temp_sorted 1>>%s 2>>%s' % (lib_settings.get_mapped_reads_sam(), lib_settings.get_mapped_reads_sam(),
        #                                                                  lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()
        subprocess.Popen('mv %s.temp_sorted.bam %s' % (lib_settings.get_mapped_reads_sam(),
                                                                          lib_settings.get_mapped_reads()), shell = True).wait()
        subprocess.Popen('samtools index %s' % (lib_settings.get_mapped_reads()), shell = True).wait()
        subprocess.Popen('rm %s' % (lib_settings.get_mapped_reads_sam()), shell = True).wait()
        lib_settings.write_to_log('mapping_reads done')
    def rdir_path(self, *args):
        return os.path.join(self.settings.get_rdir(), *args)

    def get_rdir_fhandle(self, *args):
        """
        returns a filehandle to the fname in the rdir
        """
        out_path = self.rdir_path(*args)
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return tps_utils.aopen(out_path, 'w')

    def perform_qc(self):
        qc_engine = tps_qc.TPS_qc(self, self.settings, self.threads)
        if self.settings.get_property('collapse_identical_reads'):
            qc_engine.plot_pcr_bias()
        qc_engine.identify_contaminating_sequences()
        qc_engine.print_library_count_concordances()
        qc_engine.plot_average_read_positions()
        qc_engine.plot_count_distributions()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    parser.add_argument("--make-tables",
                        help="Makes tables.",
                        action='store_true')
    parser.add_argument("--perform-qc",
                        help="performs quality control analysis.",
                        action='store_true')
    parser.add_argument("--make-plots",
                        help="Makes plots.",
                        action='store_true')
    parser.add_argument("--comparisons",
                        help="Does comparisons to other experiments",
                        action='store_true')
    parser.add_argument("--all-tasks",
                        help="Makes plots, tables, folding and comparisons",
                        action='store_true')
    parser.add_argument("--threads",
                        help="Max number of processes to use",
                        type = int, default = 8)
    args = parser.parse_args()

    return args

def main():
    """
    """
    args = parse_args()
    settings = tps_settings.TPS_settings(args.settings_file)
    tps_experiment = TPSe(settings, args.threads)
    print 'TPSe ready'
    if args.perform_qc or args.all_tasks:
        print 'QC'
        settings.write_to_log('performing QC')
        tps_experiment.perform_qc()
        settings.write_to_log('done performing QC')
    if args.make_tables or args.all_tasks:
        print 'tables'
        settings.write_to_log('making tables')
        tps_experiment.make_tables()
        settings.write_to_log('done making tables')
    if args.make_plots or args.all_tasks:
        print 'plots'
        settings.write_to_log('making plots')
        tps_experiment.make_plots()
        settings.write_to_log('done making plots')

    if args.comparisons or args.all_tasks:
        settings.write_to_log('doing comparisons')
        tps_experiment.compare_all_other_experiments()

main()