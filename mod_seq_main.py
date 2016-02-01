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
import subprocess

import mod_settings
import mod_utils
import mod_lib
import mod_qc
import mod_plotting


class mod_seq_run:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.remove_adaptor()
        self.trim_reads()
        self.create_shapemapper_settings()
        self.run_shapemapper()
        self.initialize_libs()
        #self.make_plots()
        self.make_tables()

    def remove_adaptor(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.adaptorless_reads_exist():
                    break
            else:
                return

        if self.settings.get_property('trim_adaptor'):
            mod_utils.make_dir(self.rdir_path('adaptor_removed'))
            mod_utils.parmap(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting),
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
        mod_utils.parmap(lambda lib_setting: self.trim_one_lib(lib_setting), self.settings.iter_lib_settings(),
                       nprocs = self.threads)
        self.settings.write_to_log('trimming reads complete')

    def trim_one_lib(self, lib_settings):
        lib_settings.write_to_log('trimming_reads')
        first_base_to_keep = self.settings.get_property('first_base_to_keep')  #the trimmer is 1-indexed. 1 means keep
                                                                               #  every base
        last_base_to_keep = self.settings.get_property('last_base_to_keep')  #Will keep entire 3' end if this is greater
                                                                             #than or equal to the read length
        if self.settings.get_property('trim_adaptor'):
            subprocess.Popen('gunzip -c %s | fastx_trimmer -f %d -Q33 -l %d -o %s >>%s 2>>%s' % (lib_settings.get_adaptor_trimmed_reads(),
                                                                                      first_base_to_keep, last_base_to_keep,
                                                                                      lib_settings.get_trimmed_reads(),
                                                                                      lib_settings.get_log(),
                                                                                      lib_settings.get_log()), shell=True).wait()
        else:
            subprocess.Popen('gunzip -c %s | fastx_trimmer -f %d -Q33 -l %d -o %s >>%s 2>>%s' % (lib_settings.get_fastq_file(),
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
        all_chromsomes = ', '.join(sorted(self.settings.rRNA_seqs.keys()))
        for line in reference_config_file:
            if line.startswith("<chromosome Identifiers go here>"):
                #this is where we map all of the library names to which chromosomes we want to map to
                for lib_settings in self.settings.iter_lib_settings():
                    output_config_file.write('%s: %s = %s\n' % (lib_settings.sample_name,
                                                                os.path.basename(lib_settings.get_trimmed_reads()), all_chromsomes))
            elif line.startswith('<profiles go here>'):
                for i in range(len(self.settings.get_property('experimentals'))):
                    output_config_file.write('name = %s_%s\n' % (sorted(self.settings.rRNA_seqs.keys())[0],
                                                                 self.settings.get_property('experimentals')[i]))
                    output_config_file.write('target = %s\n' % (sorted(self.settings.rRNA_seqs.keys())[0]))
                    output_config_file.write('plus_reagent = %s\n' % (self.settings.get_property('experimentals')[i]))
                    output_config_file.write('minus_reagent = %s\n' % (self.settings.get_property('no_mod_controls')[i]))
                    output_config_file.write('denat_control = %s\n\n' % (self.settings.get_property('with_mod_controls')[i]))
            else:
                output_config_file.write(line)
        reference_config_file.close()
        output_config_file.close()
        #shapemapper needs an individual FASTA file for each RNA seq that's being mapped to
        for rna_name in self.settings.rRNA_seqs:
            f = open(os.path.join(os.path.dirname(lib_settings.get_trimmed_reads()), rna_name+'.fa'), 'w')
            f.write('>%s\n' % rna_name)
            f.write(self.settings.rRNA_seqs[rna_name])
            f.close()
        self.settings.write_to_log('done creating shapemapper config file and fasta files')

    def need_to_run_shapemapper(self):
        if self.settings.get_property('force_shapemapper'):
            return True
        else:
            shapemapper_output_dir = os.path.join(os.path.dirname(self.settings.get_shapemapper_config_file()),
                                                  'output', 'counted_mutations_columns')
            for sample_name in self.settings.get_property('experimentals') + self.settings.get_property(
                    'no_mod_controls')+ self.settings.get_property('with_mod_controls'):
                for rRNA_name in self.settings.rRNA_seqs:
                    expected_file_name = os.path.join(shapemapper_output_dir, sample_name+'_'+rRNA_name+'.csv')
                    if not mod_utils.file_exists(expected_file_name):
                        return True
            return False

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
        lib = mod_lib.ModLib(self, self.settings, lib_settings)
        self.libs.append(lib)

    def get_lib_from_name(self, normalizing_lib_name):
        for lib in self.libs:
            if lib.lib_settings.sample_name == normalizing_lib_name:
                return lib
        return None

    def get_normalizable_libs(self):
        normalizeable_libs = []
        for lib in self.libs:
            if lib.lib_settings.sample_name in self.settings.get_property('experimentals'):
                normalizeable_libs.append(lib)
        return normalizeable_libs

    def make_tables(self):
        mod_utils.make_dir(self.rdir_path('tables'))
        self.pickle_mutation_rates('mutation_rates.pkl')
        self.pickle_mutation_rates('back_subtracted_mutation_rates.pkl', subtract_background=True)
        self.pickle_mutation_rates('control_subtracted_mutation_rates.pkl', subtract_control=True)
        self.write_wigs('')
        self.write_wigs('back_subtract', subtract_background=True)
        self.write_wigs('control_subtract', subtract_control=True)
        self.write_mutation_rates_tsv('mutation_rates.tsv')
        self.write_mutation_rates_tsv('back_subtracted_mutation_rates.tsv', subtract_background=True)
        self.write_mutation_rates_tsv('control_subtracted_mutation_rates.tsv', subtract_control=True)


    def write_mutation_rates_tsv(self, suffix, subtract_background=False, subtract_control=False):

        if subtract_background or subtract_control:
            libs_to_write = self.get_normalizable_libs()
        else:
            libs_to_write = self.libs

        for lib in libs_to_write:
            lib.write_tsv_tables(os.path.join(self.rdir_path('tables'), lib.lib_settings.sample_name+'_'+suffix),
                                 subtract_background=subtract_background, subtract_control=subtract_control)





    def write_wigs(self, suffix, subtract_background=False, subtract_control=False):
        mod_utils.make_dir(self.rdir_path('wigs'))
        if subtract_background or subtract_control:
            libs_to_write = self.get_normalizable_libs()
        else:
            libs_to_write = self.libs
        #will also write a file to make batch import into mochiview easier
        f = open(os.path.join(self.rdir_path('wigs'), 'mochi_batch_'+suffix+'.txt'), 'w')
        f.write('SEQUENCE_SET\tFILE_NAME\tDATA_TYPE\tNAME\n')
        for lib in libs_to_write:
            f.write('<replace>\t%s\t<replace>\t%s\n' % (lib.lib_settings.sample_name+'_'+suffix+'.wig.gz', lib.lib_settings.sample_name+'_'+suffix))
            lib.write_mutation_rates_to_wig(os.path.join(self.rdir_path('wigs'), lib.lib_settings.sample_name+'_'+suffix),
                                      subtract_background=subtract_background, subtract_control=subtract_control)
        f.close()

    def pickle_mutation_rates(self, suffix, subtract_background=False, subtract_control=False):
        if subtract_background or subtract_control:
            libs_to_pickle = self.get_normalizable_libs()
        else:
            libs_to_pickle = self.libs
        for lib in libs_to_pickle:
            lib.pickle_mutation_rates(os.path.join(self.rdir_path('tables'), lib.lib_settings.sample_name+'_'+suffix),
                                      subtract_background=subtract_background, subtract_control=subtract_control)

    def make_plots(self):
        mod_utils.make_dir(self.rdir_path('plots'))
        mod_plotting.plot_mutated_nts_pie(self.libs, os.path.join(self.rdir_path('plots'), 'raw_mutation_fractions'))
        mod_plotting.plot_mutated_nts_pie(self.libs,
                                          os.path.join(self.rdir_path('plots'),
                                                       'background_sub_mutation_fractions'), subtract_background = True)
        mod_plotting.plot_mutated_nts_pie(self.libs,
                                          os.path.join(self.rdir_path('plots'),
                                                       'control_sub_mutation_fractions'), subtract_control = True)
        mod_plotting.plot_mutation_rate_cdfs(self.libs, os.path.join(self.rdir_path('plots'), 'mutation_rate_cdf'),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'))




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
            mod_utils.parmap(lambda lib_setting: self.collapse_one_fastq_file(lib_setting), self.settings.iter_lib_settings(),
                             nprocs = self.threads)
        else:
            mod_utils.parmap(lambda lib_setting: self.fastq_to_fasta(lib_setting), self.settings.iter_lib_settings(),
                             nprocs = self.threads)
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