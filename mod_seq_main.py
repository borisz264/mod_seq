
__author__ = 'Boris Zinshteyn'
"""
Intended for processing of DMS or other chemical probing data on rRNA
Based on Alex Robertson's original RNA Bind n Seq pipeline, available on github
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42  #leaves most text as actual text in PDFs, not outlines
import os
import argparse
import subprocess

import mod_settings
import mod_utils
import mod_lib
import mod_plotting


class mod_seq_run:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.remove_adaptor()
        self.trim_reads()
        self.run_shapemapper()
        self.generate_mapping_index()
        self.map_reads()
        self.initialize_libs()
        self.make_plots()
        self.make_plots(exclude_constitutive=True)
        self.make_tables()
        self.make_tables(exclude_constitutive=True)
        #self.annotate_structures()
        #self.annotate_structures(exclude_constitutive=True)

    def remove_adaptor(self):
        self.settings.write_to_log('removing adaptors with skewer')
        for lib_settings in self.settings.iter_lib_settings():
            if not lib_settings.adaptorless_reads_exist():
                break
        else:
            self.settings.write_to_log('using existing adaptor-trimmed reads')
            return
        mod_utils.make_dir(self.rdir_path('adaptor_removed'))
        num_datasets = len([lib for lib in self.settings.iter_lib_settings()])
        num_instances = min(num_datasets, self.threads)
        threads_per_instance = self.threads/num_instances
        mod_utils.parmap(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting, threads_per_instance), self.settings.iter_lib_settings(), nprocs=num_instances)
        self.settings.write_to_log('removing adaptors done')

    def remove_adaptor_one_lib(self, lib_settings, threads):
        lib_settings.write_to_log('adaptor trimming')
        """
        -x specifies the 3' adaptor to trim from the forward read
        -Q specifies the lowest acceptable mean read quality before trimming
        -l specifies the minimum post-trimming read length
        -L specifies the maximum post-trimming read length
        -o is the output prefix
        --threads specifies number of threads to use
        """
        command_to_run = 'skewer -x %s  -k 1 -o %s --quiet --threads %s %s 1>>%s 2>>%s' % (
            self.settings.get_property('adaptor_sequence'),
            lib_settings.get_adaptor_trimmed_reads(prefix_only=True),
            threads,
            lib_settings.get_fastq_file(),
            lib_settings.get_log(), lib_settings.get_log())
        subprocess.Popen(command_to_run, shell=True).wait()
        subprocess.Popen('gzip %s-trimmed.fastq' % (lib_settings.get_adaptor_trimmed_reads(prefix_only=True)), shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')

    def trim_reads(self):
        """
        Trim reads by given amount, removing potential random barcoding sequences from 5' end
        Trimming from 3' end can also help if mapping is problematic by reducing chance for indels to prevent mapping
        :return:
        """
        self.settings.write_to_log( 'trimming reads with seqtk')
        for lib_settings in self.settings.iter_lib_settings():
            if not lib_settings.trimmed_reads_exist():
                break
        else:
            self.settings.write_to_log('using existing trimmed reads')
            return
        mod_utils.make_dir(self.rdir_path('trimmed_reads'))
        mod_utils.parmap(lambda lib_setting: self.trim_one_lib(lib_setting), self.settings.iter_lib_settings(),
                       nprocs = self.threads)
        self.settings.write_to_log('trimming reads complete')

    def trim_one_lib(self, lib_settings):
        lib_settings.write_to_log('trimming_reads')
        bases_to_trim = self.settings.get_property('first_base_to_keep')-1
        subprocess.Popen('seqtk trimfq -b %d -e 0 %s | gzip > %s 2>>%s' % (bases_to_trim,
                                                                                 lib_settings.get_adaptor_trimmed_reads(),
                                                                                 lib_settings.get_trimmed_reads(),
                                                                                 lib_settings.get_log()), shell=True).wait()
        lib_settings.write_to_log('trimming_reads done')

    def need_to_run_shapemapper(self):
        for lib_setting in self.settings.iter_lib_settings():
            for rRNA_name in self.settings.rRNA_seqs:
                expected_file_name = os.path.join(lib_setting.get_shapemapper_out_dir(), 'Pipeline_Modified_'+rRNA_name+'_mutation_counts.txt')
                if not mod_utils.file_exists(expected_file_name):
                    return True
        return False

    def run_shapemapper(self):
        """
        runs shapemapper2.0 on the samples in batches
        :return:
        """
        self.settings.write_to_log('running shapemapper')
        if self.need_to_run_shapemapper():
            mod_utils.make_dir(self.rdir_path('shapemapper'))
            all_settings = [lib_setting for lib_setting in self.settings.iter_lib_settings()]
            num_datasets = len(all_settings)
            num_instances = min(num_datasets, self.threads)
            threads_per_instance = self.threads/num_instances
            mod_utils.parmap(lambda lib_setting: self.run_single_shapemapper(lib_setting, threads_per_instance), all_settings, nprocs=num_instances)
        else:
            self.settings.write_to_log('using existing shapemapper output')
        self.settings.write_to_log('done running shapemapper')

    def run_single_shapemapper(self, lib_setting, threads_per_instance):
        '''
        :param lib_setting: 
        :param threads_per_instance: 
        :return: 
        '''
        command_to_run = 'shapemapper --target %s --name %s --nproc %d --output-counted-mutations --out %s --temp %s --log %s --min-depth 1 --overwrite --modified --U %s' %\
                         (self.settings.get_property('rrna_fasta'),
                          lib_setting.sample_name,
                          threads_per_instance,
                          lib_setting.get_shapemapper_out_dir(),
                          lib_setting.get_shapemapper_temp_dir(),
                          lib_setting.get_shapemapper_log(),
                          lib_setting.get_trimmed_reads())
        subprocess.Popen(command_to_run, shell=True).wait()

    def generate_mapping_index(self):
        """
        builds a STAR index from the input fasta file
        """
        self.settings.write_to_log('building STAR index')
        if not self.settings.star_index_exists():
            mod_utils.make_dir(self.settings.get_star_index())
            subprocess.Popen('STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --genomeSAindexNbases 4 1>>%s 2>>%s' %
                             (self.threads, self.settings.get_star_index(), self.settings.get_rRNA_fasta(), self.settings.get_log(), self.settings.get_log()), shell=True).wait()
        self.settings.write_to_log('building STAR index complete')

    def map_reads(self):
        """
        map all reads using STAR
        :return:
        """
        self.settings.write_to_log('mapping reads')
        for lib_settings in self.settings.iter_lib_settings():
            if not lib_settings.mapped_reads_exist():
                break
        else:
            return
        mod_utils.make_dir(self.rdir_path('mapped_reads'))
        all_settings = [lib_setting for lib_setting in self.settings.iter_lib_settings()]
        num_datasets = len(all_settings)
        num_instances = min(num_datasets, self.threads)
        threads_per_instance = self.threads/num_instances
        mod_utils.parmap(lambda lib_setting: self.map_one_library(lib_setting, threads_per_instance), all_settings, nprocs=num_instances)
        self.settings.write_to_log( 'finished mapping reads')

    def map_one_library(self, lib_settings, threads):
        lib_settings.write_to_log('mapping_reads')
        command_to_run = 'STAR --limitBAMsortRAM 20000000000 --runThreadN %d --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c --outWigNorm None --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 1 --outWigType wiggle read1_5p --outFileNamePrefix %s --outReadsUnmapped FastX 1>>%s 2>>%s' %\
                         (threads, self.settings.get_star_index(), lib_settings.get_trimmed_reads(),
                          lib_settings.get_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())

        subprocess.Popen(command_to_run, shell=True).wait()
        subprocess.Popen('samtools index %s' % (lib_settings.get_mapped_reads()), shell=True).wait()
        lib_settings.write_to_log('mapping_reads done')

    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
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

    def get_modified_libs(self):
        modified_libs = []
        for lib in self.libs:
            if (lib.lib_settings.sample_name in
                    self.settings.get_property('experimentals')) or (lib.lib_settings.sample_name in self.settings.get_property('with_mod_controls')):
                modified_libs.append(lib)
        return modified_libs

    def make_tables(self, exclude_constitutive=False):
        #subfolders = ['raw', 'background_subtracted', 'control_subtracted', 'fold_change']
        subfolders = ['raw', 'fold_change']
        for subfolder in subfolders:
            mod_utils.make_dir(self.rdir_path('tables', subfolder))
            mod_utils.make_dir(self.rdir_path('pickles', subfolder))
            mod_utils.make_dir(self.rdir_path('tables', subfolder, 'exclude_constitutive'))
            mod_utils.make_dir(self.rdir_path('pickles', subfolder, 'exclude_constitutive'))
        self.pickle_mutation_rates('mutation_rates.pkl', exclude_constitutive=exclude_constitutive)
        #self.pickle_mutation_rates('back_subtracted_mutation_rates.pkl', subtract_background=True, exclude_constitutive=exclude_constitutive)
        #self.pickle_mutation_rates('control_subtracted_mutation_rates.pkl', subtract_control=True, exclude_constitutive=exclude_constitutive)
        #self.pickle_fold_changes('mutation_rate_fold_changes.pkl', exclude_constitutive=True)
        self.write_wigs('')
        self.write_wigs('back_subtract', subtract_background=True)
        self.write_wigs('control_subtract', subtract_control=True)
        self.write_mutation_rates_tsv('mutation_rates.tsv', exclude_constitutive=exclude_constitutive)
        #self.write_mutation_rates_tsv('back_subtracted_mutation_rates.tsv', subtract_background=True, exclude_constitutive=exclude_constitutive)
        self.write_mutation_rates_tsv('control_subtracted_mutation_rates_lowess.tsv', subtract_control=True, exclude_constitutive=exclude_constitutive, lowess_correct = True)
        #self.write_mutation_rates_tsv('lowess_control_subtracted_mutation_rates.tsv', subtract_control=True,
        #                              exclude_constitutive=exclude_constitutive, lowess_correct=True)
        self.write_combined_mutation_rates_tsv()
        self.write_combined_mutation_rates_tsv(exclude_constitutive=True)
        self.write_combined_mutation_counts_tsv()
        self.write_combined_mutation_counts_tsv(exclude_constitutive=True)

    def write_mutation_rates_tsv(self, suffix, subtract_background=False, subtract_control=False, exclude_constitutive=False, lowess_correct=False):
        if subtract_background or subtract_control:
            libs_to_write = self.get_normalizable_libs()
        else:
            libs_to_write = self.libs
        if subtract_background == False and subtract_control == False:
            prefix = 'raw'
        elif subtract_background == True and subtract_control == False:
            prefix = 'background_subtracted'
        elif subtract_background == False and subtract_control == True:
            prefix = 'fold_change'

        if exclude_constitutive:
            for lib in libs_to_write:
                lib.write_tsv_tables(os.path.join(self.rdir_path('tables', prefix, 'exclude_constitutive'),
                                                  lib.lib_settings.sample_name+'_'+suffix[:-4]+'_exclude_constitutive'+suffix[-4:]),
                                     subtract_background=subtract_background, subtract_control=subtract_control, exclude_constitutive=exclude_constitutive, lowess_correct=lowess_correct)
        else:
            for lib in libs_to_write:
                lib.write_tsv_tables(os.path.join(self.rdir_path('tables', prefix), lib.lib_settings.sample_name+'_'+suffix),
                                     subtract_background=subtract_background, subtract_control=subtract_control, exclude_constitutive=exclude_constitutive, lowess_correct=lowess_correct)

    def write_combined_mutation_rates_tsv(self, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        if subtract_background and subtract_control:
            raise SyntaxError('Cannot subtract background and control simultaneously')

        if subtract_background or subtract_control:
            libs_to_write = list(self.get_normalizable_libs())
        else:
            libs_to_write = list(self.libs)

        if subtract_background == False and subtract_control == False:
            prefix = 'raw_'
        elif subtract_background == True and subtract_control == False:
            prefix = 'background_subtracted_'
        elif subtract_background == False and subtract_control == True:
            prefix = 'control_subtracted_'
        if exclude_constitutive:
            f = open(self.rdir_path('tables', prefix+'all_datasets_exclude_constitutive.tsv'), 'w')

        else:
            f = open(self.rdir_path('tables', prefix+'all_datasets.tsv'), 'w')
        f.write('rRNA\tposition\tnucleotide\t%s\n' % ('\t'.join([lib.lib_settings.sample_name for lib in libs_to_write])))
        for rRNA_name in sorted(self.settings.rRNA_seqs.keys()):
            for position in range(len(self.settings.rRNA_seqs[rRNA_name])):
                nuc_identity = self.settings.rRNA_seqs[rRNA_name][position]
                nuc_values = []
                for lib in libs_to_write:
                    nucleotide = lib.get_nucleotide(rRNA_name, position+1)
                    assert nucleotide.identity == nuc_identity
                    if not subtract_background and not subtract_control:
                        nuc_values.append(nucleotide.mutation_rate)
                    elif subtract_background:
                        nuc_values.append(nucleotide.get_back_sub_mutation_rate())
                    elif subtract_control:
                        nuc_values.append(nucleotide.get_control_sub_mutation_rate())
                assert len(nuc_values) == len(libs_to_write)
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    f.write('%s\t%d\t%s\t%s\n' % (rRNA_name, position+1, nuc_identity, '\t'.join(['' for nuc_value in nuc_values])))
                else:
                    f.write('%s\t%d\t%s\t%s\n' % (rRNA_name, position+1, nuc_identity, '\t'.join([str(nuc_value) for nuc_value in nuc_values])))
        f.close()

    def write_combined_mutation_counts_tsv(self, exclude_constitutive=False):
        libs_to_write = list(self.libs)
        prefix = 'raw_'

        if exclude_constitutive:
            f = open(self.rdir_path('tables', prefix+'mutation_counts_exclude_constitutive.tsv'), 'w')
        else:
            f = open(self.rdir_path('tables', prefix+'mutation_counts.tsv'), 'w')
        f.write('rRNA\tposition\tnucleotide\t%s\n' % ('\t'.join([lib.lib_settings.sample_name for lib in libs_to_write])))
        for rRNA_name in sorted(self.settings.rRNA_seqs.keys()):
            for position in range(len(self.settings.rRNA_seqs[rRNA_name])):
                nuc_identity = self.settings.rRNA_seqs[rRNA_name][position]
                nuc_values = []
                for lib in libs_to_write:
                    nucleotide = lib.get_nucleotide(rRNA_name, position+1)
                    assert nucleotide.identity == nuc_identity
                    nuc_values.append(nucleotide.total_mutation_counts)
                assert len(nuc_values) == len(libs_to_write)
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    f.write('%s\t%d\t%s\t%s\n' % (rRNA_name, position+1, nuc_identity, '\t'.join(['' for nuc_value in nuc_values])))
                else:
                    f.write('%s\t%d\t%s\t%s\n' % (rRNA_name, position+1, nuc_identity, '\t'.join([str(nuc_value) for nuc_value in nuc_values])))
        f.close()

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

    def pickle_mutation_rates(self, suffix, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        if subtract_background:
            libs_to_pickle = self.get_modified_libs()
        elif subtract_control:
            libs_to_pickle = self.get_normalizable_libs()
        else:
            libs_to_pickle = self.libs

        if subtract_background == False and subtract_control == False:
            prefix = 'raw'
        elif subtract_background == True and subtract_control == False:
            prefix = 'background_subtracted'
        elif subtract_background == False and subtract_control == True:
            prefix = 'control_subtracted'

        if exclude_constitutive:
            for lib in libs_to_pickle:
                lib.pickle_mutation_rates(os.path.join(self.rdir_path('pickles', prefix, 'exclude_constitutive'),
                                                       lib.lib_settings.sample_name+'_'+suffix[:-4]+'_exclude_constitutive'+suffix[-4:]),
                                                       subtract_background=subtract_background, subtract_control=subtract_control,
                                                       exclude_constitutive=exclude_constitutive)
        else:
            for lib in libs_to_pickle:
                lib.pickle_mutation_rates(os.path.join(self.rdir_path('pickles', prefix), lib.lib_settings.sample_name+'_'+suffix),
                                          subtract_background=subtract_background, subtract_control=subtract_control, exclude_constitutive=exclude_constitutive)

    def pickle_fold_changes(self, suffix, exclude_constitutive=False):
        libs_to_pickle = self.get_normalizable_libs()

        if exclude_constitutive:
            for lib in libs_to_pickle:
                lib.pickle_mutation_fold_change(os.path.join(self.rdir_path('pickles', 'fold_change', 'exclude_constitutive'),
                                                       lib.lib_settings.sample_name + '_' + suffix[
                                                                                            :-4] + '_exclude_constitutive' + suffix[
                                                                                                                             -4:]), exclude_constitutive=exclude_constitutive)
        else:
            for lib in libs_to_pickle:
                lib.pickle_mutation_fold_change(
                    os.path.join(self.rdir_path('pickles', 'fold_change'), lib.lib_settings.sample_name + '_' + suffix), exclude_constitutive=exclude_constitutive)


    def make_plots(self, exclude_constitutive=False):
        if exclude_constitutive:
            mod_utils.make_dir(self.rdir_path('plots', 'exclude_constitutive'))
            mod_utils.make_dir(self.rdir_path('plots', 'exclude_constitutive', 'functional_groups'))
            mod_utils.make_dir(self.rdir_path('plots', 'exclude_constitutive', 'interactive'))
            rdir = self.rdir_path('plots','exclude_constitutive')
            file_tag = '_exclude_constitutive'
            #TODO: the names for the ROC curve chromosomes are hard coded and need to be changed between samples
            #mod_plotting.generate_roc_curves(self.settings.get_property('tptn_file_25s'), self.settings.rRNA_seqs, os.path.join(rdir, '23S_ROC_curves'), self.get_modified_libs(), 'E.c.23S_rRNA', self.settings.get_property('affected_nucleotides'))
            #mod_plotting.generate_roc_curves(self.settings.get_property('tptn_file_18s'), self.settings.rRNA_seqs, os.path.join(rdir, '16S_ROC_curves'), self.get_modified_libs(), 'E.c.16S_rRNA', self.settings.get_property('affected_nucleotides'))
        else:
            mod_utils.make_dir(self.rdir_path('plots'))
            mod_utils.make_dir(self.rdir_path('plots', 'interactive'))
            rdir = self.rdir_path('plots')
            file_tag = ''

        mod_plotting.plot_mutated_nts_pie(self.libs, os.path.join(rdir, 'raw_mutation_fractions'+file_tag), exclude_constitutive=exclude_constitutive)
        mod_plotting.plot_rt_stop_pie(self.libs, os.path.join(rdir, 'raw_rt_stops'+file_tag), exclude_constitutive=exclude_constitutive)
        mod_plotting.plot_mutation_breakdown_pie(self.libs, os.path.join(rdir, 'raw_mutation_types'+file_tag), exclude_constitutive=exclude_constitutive)

        mod_plotting.plot_mutated_nts_pie(self.libs,
                                          os.path.join(rdir, 'background_sub_mutation_fractions'+file_tag),
                                          subtract_background = True, exclude_constitutive=exclude_constitutive)
        mod_plotting.plot_rt_stop_pie(self.libs, os.path.join(rdir, 'back_sub_rt_stops'+file_tag), subtract_background = True, exclude_constitutive=exclude_constitutive)

        mod_plotting.plot_mutation_rate_cdfs(self.libs, os.path.join(rdir, 'mutation_rate_cdf'+file_tag),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)

        mod_plotting.plot_mutation_rate_violins(self.libs, os.path.join(rdir, 'mutation_rate_cdf'+file_tag),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)

        mod_plotting.ma_plots(self.get_normalizable_libs(), os.path.join(rdir, 'MA'+file_tag),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)
        mod_plotting.ma_plots_by_count(self.get_normalizable_libs(), os.path.join(rdir, 'MA_raw_counts'+file_tag),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)
        mod_plotting.ma_plots_by_count(self.get_normalizable_libs(), os.path.join(rdir, 'MA_raw_counts_lowess'+file_tag),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive, lowess_correct=True)
        mod_plotting.mutation_rate_scatter(self.get_normalizable_libs(), os.path.join(rdir, 'scatter_mismatch_rate'+file_tag),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)

        if self.settings.get_property('make_interactive_plots'):
                mod_plotting.ma_plots_interactive(self.get_normalizable_libs(), os.path.join(rdir, 'interactive', 'MA'+file_tag),
                                                         nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                                         exclude_constitutive=False)
                mod_plotting.scatter_interactive(self.get_normalizable_libs(),
                                                  os.path.join(rdir, 'interactive', 'scatter' + file_tag),
                                                  nucleotides_to_count=self.settings.get_property(
                                                      'affected_nucleotides'),
                                                  exclude_constitutive=False)
                mod_plotting.ma_plots_interactive_by_count(self.get_normalizable_libs(),
                                                  os.path.join(rdir, 'interactive', 'MA_counts' + file_tag),
                                                  nucleotides_to_count=self.settings.get_property(
                                                      'affected_nucleotides'),
                                                  exclude_constitutive=False)
                mod_plotting.ma_plots_interactive_by_count(self.get_normalizable_libs(),
                                                  os.path.join(rdir, 'interactive', 'MA_counts_lowess' + file_tag),
                                                  nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                                           exclude_constitutive=False, lowess_correct=True)

    def annotate_structures(self, exclude_constitutive=False):
        if exclude_constitutive:
            mod_utils.make_dir(self.rdir_path('structures', 'protections_highlighted', 'exclude_constitutive'))
            mod_utils.make_dir(self.rdir_path('structures', 'colored_by_change', 'exclude_constitutive'))
            file_tag = '_exclude_constitutive'
        else:
            mod_utils.make_dir(self.rdir_path('structures', 'protections_highlighted'))
            mod_utils.make_dir(self.rdir_path('structures', 'colored_by_change'))
            file_tag = ''
        if exclude_constitutive:
            mod_plotting.highlight_structure(self.get_normalizable_libs(), self.rdir_path('structures', 'protections_highlighted', 'exclude_constitutive'),
                                             nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                             exclude_constitutive=exclude_constitutive)
            # mod_plotting.color_by_change(self.get_normalizable_libs(), self.rdir_path('structures', 'colored_by_change', 'exclude_constitutive'),
            #                              nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
            #                              exclude_constitutive=exclude_constitutive)
        else:
            mod_plotting.highlight_structure(self.get_normalizable_libs(), self.rdir_path('structures', 'protections_highlighted'),
                                 nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
                                 exclude_constitutive=exclude_constitutive)
            # mod_plotting.color_by_change(self.get_normalizable_libs(), self.rdir_path('structures', 'colored_by_change'),
            #                              nucleotides_to_count=self.settings.get_property('affected_nucleotides'),
            #                              exclude_constitutive=exclude_constitutive)

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

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
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


main()