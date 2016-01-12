import operator

__author__ = 'Boris Zinshteyn'
"""
Intended for processing of DMS or other chemical probing data on rRNA
Based on Alex Robertson's original RBNS pipeline, available on github
"""
import sys
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42  #leaves most text as actual text in PDFs, not outlines
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
import mod_settings
import mod_utils
import mod_lib
import mod_qc
import count_reads_and_mismatches


class mod_seq_run:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.remove_adaptor()
        self.trim_reads()
        self.create_shapemapper_settings()
        self.run_shapemapper()

    def remove_adaptor(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.adaptorless_reads_exist():
                    break
            else:
                return

        if self.settings.get_property('trim_adaptor'):
            mod_utils.make_dir(self.rdir_path('adaptor_removed'))
            bzUtils.parmap(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting),
                           self.settings.iter_lib_settings(), nprocs=self.threads)

    def remove_adaptor_one_lib(self, lib_settings):
        lib_settings.write_to_log('adaptor trimming')
        if self.settings.get_property('discard_untrimmed'):
            command_to_run = 'cutadapt --adapter %s --overlap 3 --discard-untrimmed --minimum-length %d %s --output %s 1>>%s 2>>%s' % (self.settings.get_property('adaptor_sequence'), self.settings.get_property('min_post_adaptor_length'),
                               lib_settings.get_fastq_file(), lib_settings.get_adaptor_trimmed_reads(), lib_settings.get_log(),
                               lib_settings.get_log())
        else:
            command_to_run = 'cutadapt --adapter %s --overlap 3 --minimum-length %d %s --output %s 1>>%s 2>>%s' % (self.settings.get_property('adaptor_sequence'), self.settings.get_property('min_post_adaptor_length'),
                   lib_settings.get_fastq_file(), lib_settings.get_adaptor_trimmed_reads(), lib_settings.get_log(),
                   lib_settings.get_log())
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')

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
        mod_utils.make_dir(self.rdir_path('trimmed_reads'))
        bzUtils.parmap(lambda lib_setting: self.trim_one_lib(lib_setting), self.settings.iter_lib_settings(),
                       nprocs = self.threads)
        self.settings.write_to_log('trimming reads complete')

    def trim_one_lib(self, lib_settings):
        lib_settings.write_to_log('trimming_reads')
        first_base_to_keep = self.settings.get_property('first_base_to_keep')  #the trimmer is 1-indexed. 1 means keep
                                                                               #  every base
        last_base_to_keep = self.settings.get_property('last_base_to_keep')  #Will keep entire 3' end if this is greater
                                                                             #than or equal to the read length
        if self.settings.get_property('trim_adaptor'):
            subprocess.Popen('gunzip -c %s | fastx_trimmer -f %d -l %d -o %s >>%s 2>>%s' % (lib_settings.get_adaptor_trimmed_reads(),
                                                                                      first_base_to_keep, last_base_to_keep,
                                                                                      lib_settings.get_trimmed_reads(),
                                                                                      lib_settings.get_log(),
                                                                                      lib_settings.get_log()), shell=True).wait()
        else:
            subprocess.Popen('gunzip -c %s | fastx_trimmer -f %d -l %d -o %s >>%s 2>>%s' % (lib_settings.get_fastq_file(),
                                                                                      first_base_to_keep, last_base_to_keep,
                                                                                      lib_settings.get_trimmed_reads(),
                                                                                      lib_settings.get_log(),
                                                                                      lib_settings.get_log()), shell=True).wait()
        lib_settings.write_to_log('trimming_reads done')

    def create_shapemapper_settings(self):
        """
        Create a preferences file from the default that will run shapemapper on all datasets.
        The comparisons are essentially unimportant, as we joust want to get the mutation
        rates that shapemapper spits out
        :return:
        """
        self.settings.write_to_log('creating shapemapper config file and fasta files')
        reference_config_file = open(self.settings.get_property('shapemapper_ref_file'))
        output_config_file = open(self.settings.get_shapemapper_config_file(), 'w')
        rRNA_seqs = mod_utils.convertFastaToDict(self.settings.get_rRNA_fasta())
        all_chromsomes = ', '.join(sorted(rRNA_seqs.keys()))
        for line in reference_config_file:
            if line.startswith("<chromosome Identifiers go here>"):
                #this is where we map all of the library names to which chromosomes we want to map to
                for lib_settings in self.settings.iter_lib_settings():
                    output_config_file.write('%s: %s = %s\n' % (lib_settings.sample_name,
                                                                os.path.basename(lib_settings.get_trimmed_reads()), all_chromsomes))
            elif line.startswith('<profiles go here>'):
                for i in range(len(self.settings.get_property('experimentals'))):
                    output_config_file.write('name = %s_%s\n' % (sorted(rRNA_seqs.keys())[0],
                                                                 self.settings.get_property('experimentals')[i]))
                    output_config_file.write('target = %s\n' % (sorted(rRNA_seqs.keys())[0]))
                    output_config_file.write('plus_reagent = %s\n' % (self.settings.get_property('experimentals')[i]))
                    output_config_file.write('minus_reagent = %s\n' % (self.settings.get_property('no_mod_controls')[i]))
                    output_config_file.write('denat_control = %s\n\n' % (self.settings.get_property('with_mod_controls')[i]))
            else:
                output_config_file.write(line)
        reference_config_file.close()
        output_config_file.close()
        #shapemapper needs an individual FASTA file for each RNA seq that's being mapped to
        for rna_name in rRNA_seqs:
            f = open(os.path.join(os.path.dirname(lib_settings.get_trimmed_reads()), rna_name+'.fa'), 'w')
            f.write('>%s\n' % rna_name)
            f.write(rRNA_seqs[rna_name])
            f.close()
        self.settings.write_to_log('done creating shapemapper config file and fasta files')

    def need_to_run_shapemapper(self):
        if self.settings.get_property('force_shapemapper'):
            return True
        else:
            rRNA_seqs = mod_utils.convertFastaToDict(self.settings.get_rRNA_fasta())
            for sample_name in len(self.settings.get_property('experimentals'))):


    def run_shapemapper(self):
        """
        runs shapemapper from the preferences file created above
        :return:
        """
        self.settings.write_to_log('running shapemapper')
        os.chdir(os.path.dirname(self.settings.get_shapemapper_config_file()))
        if self.need_to_run_shapemapper():
            subprocess.Popen('ShapeMapper.py %s' % (self.settings.get_shapemapper_config_file()), shell=True).wait()
        self.settings.write_to_log('done running shapemapper')

    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
        mod_utils.make_dir(self.rdir_path('sequence_counts'))
        self.libs = []
        map(lambda lib_settings: self.initialize_lib(lib_settings), self.settings.iter_lib_settings())
        self.settings.write_to_log('initializing libraries, counting reads, done')


    def initialize_lib(self, lib_settings):
        lib = mod_lib.ModLib(self.settings, lib_settings)
        self.libs.append(lib)


    def make_tables(self):
        mod_utils.make_dir(self.rdir_path('tables'))

    def make_plots(self):
        mod_utils.make_dir(self.rdir_path('plots'))



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
        mod_utils.make_dir(self.rdir_path('collapsed_reads'))
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

    def plot_rRNA_read_distributions(self, sequence_name):
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

    def get_barcode_match(self, barcode, barcodes):
        """
        takes a barcode and returns the one it matches (hamming <= 1)
        else
        empty string
        """
        if barcode in barcodes:
            return barcode
        for barcode_j in barcodes:
            if mod_utils.hamming_N(barcode, barcode_j) <= self.settings.get_property('mismatches_allowed_in_barcode'):
                return barcode_j
        return ''

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
        return mod_utils.aopen(out_path, 'w')

    def perform_qc(self):
        qc_engine = mod_qc.TPS_qc(self, self.settings, self.threads)
        #qc_engine.identify_contaminating_sequences()
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
    settings = mod_settings.mod_settings(args.settings_file)
    all_datasets = mod_seq_run(settings, args.threads)
    if args.perform_qc or args.all_tasks:
        print 'QC'
        settings.write_to_log('performing QC')
        all_datasets.perform_qc()
        settings.write_to_log('done performing QC')
    if args.make_tables or args.all_tasks:
        print 'tables'
        settings.write_to_log('making tables')
        all_datasets.make_tables()
        settings.write_to_log('done making tables')
    if args.make_plots or args.all_tasks:
        print 'plots'
        settings.write_to_log('making plots')
        all_datasets.make_plots()
        settings.write_to_log('done making plots')

    if args.comparisons or args.all_tasks:
        settings.write_to_log('doing comparisons')
        all_datasets.compare_all_other_experiments()

main()