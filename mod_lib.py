from collections import defaultdict
import os
import mod_utils
import gzip
import numpy as np
import math
from statsmodels.nonparametric.smoothers_lowess import lowess
import scipy.stats.mstats as mstats


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
        self.assign_rt_stops()
        self.winsorize_rt_stops()

    def parse_shapemapper_output_files(self):
        for rRNA_name in self.experiment_settings.rRNA_seqs:
            shapemapper_output_file = os.path.join(self.lib_settings.get_shapemapper_out_dir(),
                                                   'Pipeline_Modified_'+rRNA_name+'_mutation_counts.txt')
            assert mod_utils.file_exists(shapemapper_output_file)
            self.rRNA_mutation_data[rRNA_name] = rRNA_mutations(self, self.lib_settings, self.experiment_settings, shapemapper_output_file, rRNA_name)

    def assign_rt_stops(self):
        read_counts = mod_utils.parse_wig(self.lib_settings.get_5p_count_wig())
        for rRNA_name in self.rRNA_mutation_data:
            for position in read_counts[rRNA_name]:
                #the rt stop is BEFORE it adds the nucleotide across from the modified one
                #so the signal in read_counts[rRNA_name][position] is caused by the mod at self.nucleotides[position-1]
                if position-1 in self.rRNA_mutation_data[rRNA_name].nucleotides:
                    nuc = self.rRNA_mutation_data[rRNA_name].nucleotides[position-1]
                    nuc.rt_stops = read_counts[rRNA_name][position]
                    self.rRNA_mutation_data[rRNA_name].total_rt_stops += nuc.rt_stops

    def winsorize_rt_stops(self):
        #winsorize the data  by setting all RT stop counts above "winsorization_upper_percentile" to the value of the RT stops at that percentile
        #collect list of all values, winsorize these to get max value, then use the max to filter RT stop counts and save new value
        #then divide all by the max to get an RT stop or reactivity score.
        all_rt_stop_counts = self.list_rt_stop_counts()
        winsorized_counts = mstats.winsorize(all_rt_stop_counts, limits=(0, 1.-self.get_property('winsorization_upper_limit')), inplace=False)
        winsorized_max = max(winsorized_counts)
        for rRNA in self.rRNA_mutation_data.values():
            for nucleotide in rRNA.nucleotides.values():
                if nucleotide.rt_stops > winsorized_max:
                    nucleotide.winsorized_rt_stops = winsorized_max
                else:
                    nucleotide.winsorized_rt_stops = nucleotide.rt_stops
                nucleotide.rt_stop_score = float(nucleotide.winsorized_rt_stops)/float(winsorized_max)

    def count_mutation_rates_by_nucleotide(self, subtract_background = False, subtract_control = False, exclude_constitutive=False):
        """
        counts, over all RNAs, the total number of mutation rates at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.
        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        total_counts = defaultdict(int)

        for rRNA_name in self.rRNA_mutation_data:
            rRNA_counts = self.rRNA_mutation_data[rRNA_name].count_mutation_rates_by_nucleotide(subtract_background=subtract_background,
                                                                                                subtract_control=subtract_control, exclude_constitutive=exclude_constitutive)
            for nucleotide_type in rRNA_counts:
                total_counts[nucleotide_type] += rRNA_counts[nucleotide_type]
        return total_counts

    def count_rt_stop_rpm_by_nucleotide(self, subtract_background = False, subtract_control = False, exclude_constitutive=False):
        """
        counts, over all RNAs, the total number of RT stops at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.
        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        total_counts = defaultdict(int)

        for rRNA_name in self.rRNA_mutation_data:
            rRNA_counts = self.rRNA_mutation_data[rRNA_name].count_rt_stop_rpm_by_nucleotide(subtract_background=subtract_background, exclude_constitutive=exclude_constitutive)
            for nucleotide_type in rRNA_counts:
                total_counts[nucleotide_type] += rRNA_counts[nucleotide_type]
        return total_counts

    def count_mutation_types_by_nucleotide(self, subtract_background = False, subtract_control = False, exclude_constitutive=False):
        """
        counts, over all RNAs, the total number of each type ofmutation rates at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.
        :return: a dict like {A: {A->G:1054}, T:{T->G:1054}, C: {C->G:1054}, G:{G->C:1054}}
        """
        total_counts = defaultdict((lambda : defaultdict(int)))

        for rRNA_name in self.rRNA_mutation_data:
            rRNA_counts = self.rRNA_mutation_data[rRNA_name].count_mutation_types_by_nucleotide(subtract_background=subtract_background,
                                                                                                subtract_control=subtract_control, exclude_constitutive=exclude_constitutive)
            for nucleotide_type in rRNA_counts:
                for mutation_type in rRNA_counts[nucleotide_type]:
                    total_counts[nucleotide_type][mutation_type] += rRNA_counts[nucleotide_type][mutation_type]
        return total_counts


    def count_mutation_rates_by_type(self, subtract_background = False, subtract_control = False, exclude_constitutive=False):
        """
        counts, over all RNAs, the total number of mutation rates at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.
        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        total_counts = defaultdict(int)

        for rRNA_name in self.rRNA_mutation_data:
            rRNA_counts = self.rRNA_mutation_data[rRNA_name].count_mutation_rates_by_nucleotide(subtract_background=subtract_background,
                                                                                                subtract_control=subtract_control, exclude_constitutive=exclude_constitutive)
            for nucleotide_type in rRNA_counts:
                total_counts[nucleotide_type] += rRNA_counts[nucleotide_type]
        return total_counts

    def list_mutation_rates(self, subtract_background = False, subtract_control = False, nucleotides_to_count = 'ATCG', exclude_constitutive=False):
        all_mutation_rates = []
        for rRNA_name in self.rRNA_mutation_data:
            all_mutation_rates.extend(self.rRNA_mutation_data[rRNA_name].
                                      list_mutation_rates(subtract_background = subtract_background, subtract_control = subtract_control,
                                                          nucleotides_to_count = nucleotides_to_count, exclude_constitutive=exclude_constitutive))
        return all_mutation_rates

    def list_rt_stop_rpms(self, subtract_background = False, subtract_control = False, nucleotides_to_count = 'ATCG', exclude_constitutive=False):
        all_rt_stop_rpms = []
        for rRNA_name in self.rRNA_mutation_data:
            all_rt_stop_rpms.extend(self.rRNA_mutation_data[rRNA_name].
                                      list_rt_stop_rpms(subtract_background = subtract_background, subtract_control = subtract_control,
                                                          nucleotides_to_count = nucleotides_to_count, exclude_constitutive=exclude_constitutive))
        return all_rt_stop_rpms

    def list_rt_stop_counts(self, nucleotides_to_count = 'ATCG'):
        all_rt_stop_counts = []
        for rRNA_name in self.rRNA_mutation_data:
            all_rt_stop_counts.extend(self.rRNA_mutation_data[rRNA_name].list_rt_stop_counts(nucleotides_to_count=nucleotides_to_count))
        return all_rt_stop_counts


    def list_fold_changes(self, nucleotides_to_count = 'ATCG', exclude_constitutive=False):
        all_mutation_rates = []
        for rRNA_name in self.rRNA_mutation_data:
            all_mutation_rates.extend(self.rRNA_mutation_data[rRNA_name].
                                      list_mutation_fold_changes(nucleotides_to_count = nucleotides_to_count, exclude_constitutive=exclude_constitutive))
        return all_mutation_rates

    def get_normalizing_lib(self):
        """
        #returns the library that is the normalization for this one (no-modification control)
        """
        if self.lib_settings.sample_name in self.experiment_settings.get_property('experimentals'):
            lib_index = self.experiment_settings.get_property('experimentals').index(self.lib_settings.sample_name)
            normalizing_lib_name = self.experiment_settings.get_property('no_mod_controls')[lib_index]
            return self.experiment.get_lib_from_name(normalizing_lib_name)
        elif self.lib_settings.sample_name in self.experiment_settings.get_property('with_mod_controls'):
            lib_index = self.experiment_settings.get_property('with_mod_controls').index(self.lib_settings.sample_name)
            normalizing_lib_name = self.experiment_settings.get_property('no_mod_controls')[lib_index]
            return self.experiment.get_lib_from_name(normalizing_lib_name)
        else:
            return None
    def get_normalizing_lib_with_mod(self):
        """
        #returns the library that is the normalization for this one (with-modification control)
        """
        if self.lib_settings.sample_name in self.experiment_settings.get_property('experimentals'):
            lib_index = self.experiment_settings.get_property('experimentals').index(self.lib_settings.sample_name)
            normalizing_lib_name = self.experiment_settings.get_property('with_mod_controls')[lib_index]
            return self.experiment.get_lib_from_name(normalizing_lib_name)
        else:
            return None

    def get_nucleotide(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position]

    def get_mutation_count_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].total_mutation_counts

    def get_coverage_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].sequencing_depth

    def get_mutation_rate_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].mutation_rate

    def get_rt_stop_rpm_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].get_rt_stop_rpm()

    def get_rt_stop_score_at_position(self, rRNA_name, position):
        return self.rRNA_mutation_data[rRNA_name].nucleotides[position].rt_stop_score


    def write_tsv_tables(self, tsv_filename, subtract_background=False, subtract_control=False, exclude_constitutive=False,
                         lowess_correct = False):

        if subtract_background and subtract_control:
            raise SyntaxError('Cannot subtract background and control simultaneously')

        f = open(tsv_filename, 'w')
        if subtract_background:
                f.write('CHROMOSOME\tPOSITION\tMUTATION_RATE\tBKGD_SUB_MUT_RATE\tBKGD_SUB_ERROR\n')
        elif subtract_control:
            if lowess_correct:
                nucleotides_to_count = self.experiment_settings.get_property('affected_nucleotides')
                self.lowess_correct_mutation_fold_changes(nucleotides_to_count=nucleotides_to_count,
                                                          exclude_constitutive=exclude_constitutive)
            f.write('CHROMOSOME\tPOSITION\tNUC\tEXP_MUTATION_RATE\tEXP_MUTATION_COUNTS\tEXP_99%_min\tEXP_99%_max\tCTRL_MUT_RATE\tCTRL_MUT_counts'
                    '\tCTRL_99%_min\tCTRL_99%_max\tEXP-CTRL\tCTRL_POISSON_SUB_ERROR\tFOLD_CHANGE\tPROTECTION_CALL\n')

        elif not subtract_background and not subtract_control:
                f.write('CHROMOSOME\tPOSITION\tMUTATION_RATE\tERROR\n')


        for rRNA_name in self.rRNA_mutation_data:
            for position in self.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = self.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if lowess_correct and nucleotide.identity not in nucleotides_to_count:
                    continue
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    if subtract_background:
                        f.write(self.rRNA_mutation_data[rRNA_name].rRNA_name+'\t'+str(nucleotide.position)+'\t'
                                +'0'+'\t'+'0'+'\t'
                                +'0'+'\n')
                    elif subtract_control:
                        f.write(self.rRNA_mutation_data[rRNA_name].rRNA_name+'\t'+str(nucleotide.position)+'\t'+
                               str(nucleotide.identity)+'\t\t\t\t\t\t\t\t\t\t\t\t\n')
                    elif not subtract_background and not subtract_control:
                        f.write(self.rRNA_mutation_data[rRNA_name].rRNA_name+'\t'+str(nucleotide.position)+'\t'
                                +'0'+'\t'+'0'+'\n')
                else:

                    if subtract_background:
                        f.write(self.rRNA_mutation_data[rRNA_name].rRNA_name+'\t'+str(nucleotide.position)+'\t'
                                +str(nucleotide.mutation_rate)+'\t'+str(nucleotide.get_back_sub_mutation_rate())+'\t'
                                +str(nucleotide.get_back_sub_error())+'\n')
                    elif subtract_control:


                        ctrl_nuc = nucleotide.get_control_nucleotide()
                        if lowess_correct:
                            fold_change = nucleotide.lowess_fc
                        else:
                            fold_change = nucleotide.get_control_fold_change_in_mutation_rate()
                        exp_wil_bottom, exp_wil_top = nucleotide.get_wilson_approximate_score_interval()
                        ctrl_wil_bottom, ctrl_wil_top = ctrl_nuc.get_wilson_approximate_score_interval()
                        f.write('%s\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n' %
                                (rRNA_name, nucleotide.position, nucleotide.identity, nucleotide.mutation_rate, nucleotide.total_mutation_counts,
                                exp_wil_bottom, exp_wil_top, ctrl_nuc.mutation_rate, ctrl_nuc.total_mutation_counts,
                                ctrl_wil_bottom, ctrl_wil_top, nucleotide.get_control_sub_mutation_rate(),
                                nucleotide.get_control_sub_error(), fold_change,
                                nucleotide.determine_protection_status(confidence_interval=self.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                       fold_change_cutoff=self.experiment_settings.get_property('fold_change_cutoff'), lowess_correct = lowess_correct)))
                    elif not subtract_background and not subtract_control:
                        f.write(self.rRNA_mutation_data[rRNA_name].rRNA_name+'\t'+str(nucleotide.position)+'\t'
                                +str(nucleotide.mutation_rate)+'\t'+str(nucleotide.get_error())+'\n')

        f.close()

    def pickle_mutation_rates(self, output_name, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        """
        stores mutation rates as a simple pickle, of {rRNA_name:{position:mutation rate}}
        :param subtract_background:
        :return:
        """
        output_dict = {}
        for rRNA in self.rRNA_mutation_data:
            output_dict[rRNA] = {}
            for position in self.rRNA_mutation_data[rRNA].nucleotides:
                nucleotide = self.rRNA_mutation_data[rRNA].nucleotides[position]
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    output_dict[rRNA][position] = 0
                else:
                    if subtract_background and subtract_control:
                        raise SyntaxError('Cannot subtract background and control simultaneously')
                    if subtract_background:
                        output_dict[rRNA][position] = max((nucleotide.mutation_rate - self.get_normalizing_lib().
                                                    get_mutation_rate_at_position(rRNA, nucleotide.position)), 0.)
                    elif subtract_control:
                        output_dict[rRNA][position] = nucleotide.mutation_rate - self.get_normalizing_lib_with_mod().get_mutation_rate_at_position(rRNA, nucleotide.position)
                    else:
                        output_dict[rRNA][position] = nucleotide.mutation_rate
        mod_utils.makePickle(output_dict, output_name)

    def pickle_mutation_fold_change(self, output_name, exclude_constitutive=False):
        """
        stores mutation rates as a simple pickle, of {rRNA_name:{position:mutation rate}}
        :param subtract_background:
        :return:
        """
        output_dict = {}
        for rRNA in self.rRNA_mutation_data:
            output_dict[rRNA] = {}
            for position in self.rRNA_mutation_data[rRNA].nucleotides:
                nucleotide = self.rRNA_mutation_data[rRNA].nucleotides[position]
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    output_dict[rRNA][position] = 1.0
                else:
                    try:
                        output_dict[rRNA][position] = nucleotide.mutation_rate/self.get_normalizing_lib_with_mod().get_mutation_rate_at_position(rRNA, nucleotide.position)
                    except:
                        output_dict[rRNA][position] = float('inf')*nucleotide.mutation_rate
        mod_utils.makePickle(output_dict, output_name)

    def write_mutation_rates_to_wig(self, output_prefix, subtract_background = False, subtract_control = False):
        """
        write out mutation rates to a wig file that can be opened with a program like IGV or mochiview,
        given the corresponding rRNA fasta as a genome, of course
        :param output_prefix:
        :param subtract_background:
        :param subtract_control
        :return:
        """
        wig = gzip.open(output_prefix+'.wig.gz', 'w')

        if subtract_background:
            wig.write('track type=wiggle_0 name=%s\n' % (self.lib_settings.sample_name+'_back_sub'))
        elif subtract_control:
            wig.write('track type=wiggle_0 name=%s\n' % (self.lib_settings.sample_name+'_control_sub'))
        elif not subtract_background and not subtract_control:
            wig.write('track type=wiggle_0 name=%s\n' % (self.lib_settings.sample_name))
        for rRNA_name in self.rRNA_mutation_data:
            if subtract_background and subtract_control:
                    raise SyntaxError('Cannot subtract background and control simultaneously')

            if subtract_background:
                wig.write('variableStep chrom=%s\n' % (rRNA_name))
                for position in sorted(self.rRNA_mutation_data[rRNA_name].nucleotides.keys()):
                    if subtract_background:
                        wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                                nucleotides[position].get_back_sub_mutation_rate()))
                    else:
                        wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                                nucleotides[position].get_back_sub_mutation_rate()))
            elif subtract_control:
                wig.write('variableStep chrom=%s\n' % (rRNA_name))
                for position in sorted(self.rRNA_mutation_data[rRNA_name].nucleotides.keys()):
                    if subtract_control:
                        wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                                nucleotides[position].get_control_sub_mutation_rate()))
                    else:
                        wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                                nucleotides[position].get_control_sub_mutation_rate()))
            else:
                wig.write('variableStep chrom=%s\n' % (rRNA_name))
                for position in sorted(self.rRNA_mutation_data[rRNA_name].nucleotides.keys()):
                    if subtract_background:
                        wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                                nucleotides[position].mutation_rate))
                    else:
                        wig.write('%d\t%f\n' % (position, self.rRNA_mutation_data[rRNA_name].
                                                nucleotides[position].mutation_rate))
        wig.close()

    def write_rt_stops_to_wig(self, output_prefix):
        """
        write out mutation rates to a wig file that can be opened with a program like IGV or mochiview,
        given the corresponding rRNA fasta as a genome, of course
        :param output_prefix:
        :param subtract_background:
        :param subtract_control
        :return:
        """
        wig = gzip.open(output_prefix+'.wig.gz', 'w')

        wig.write('track type=wiggle_0 name=%s\n' % (self.lib_settings.sample_name))
        for rRNA_name in self.rRNA_mutation_data:
                wig.write('variableStep chrom=%s\n' % (rRNA_name))
                for position in sorted(self.rRNA_mutation_data[rRNA_name].nucleotides.keys()):
                        wig.write('%d\t%f\n' % (position+1, self.get_rt_stop_rpm_at_position(rRNA_name, position)))
        wig.close()

    def write_rt_stop_scores_to_wig(self, output_prefix):
        """
        write out mutation rates to a wig file that can be opened with a program like IGV or mochiview,
        given the corresponding rRNA fasta as a genome, of course
        :param output_prefix:
        :param subtract_background:
        :param subtract_control
        :return:
        """
        wig = gzip.open(output_prefix+'.wig.gz', 'w')

        wig.write('track type=wiggle_0 name=%s\n' % (self.lib_settings.sample_name))
        for rRNA_name in self.rRNA_mutation_data:
                wig.write('variableStep chrom=%s\n' % (rRNA_name))
                for position in sorted(self.rRNA_mutation_data[rRNA_name].nucleotides.keys()):
                        wig.write('%d\t%f\n' % (position+1, self.get_rt_stop_score_at_position(rRNA_name, position)))
        wig.close()

    def get_changed_nucleotides(self, change_type, nucleotides_to_count='ATCG', exclude_constitutive=False,
                                  confidence_interval = 0.99, fold_change_cutoff = 3, subtract_background=False):
        changed_nucleotides = {}
        for rRNA_name in self.rRNA_mutation_data:
            changed_nucleotides[rRNA_name] = self.rRNA_mutation_data[rRNA_name].\
                get_changed_nucleotides(change_type, nucleotides_to_count=nucleotides_to_count,
                                        exclude_constitutive=exclude_constitutive,
                                        confidence_interval = confidence_interval,
                                        fold_change_cutoff = fold_change_cutoff,
                                        subtract_background=subtract_background)
        return changed_nucleotides

    def get_nucleotides_from_list(self, nucleotide_list, nucleotides_to_count = 'ATCG', exclude_constitutive=False):
        """

        :param nucleotide_list: a list of nucleotide-identifying strings like: 'S.c.18S_rRNA 2125 A'
        :return: a list of the nucleotide objects matching those strings
        """
        nucleotides = []
        for nucleotide_string in nucleotide_list:
            rRNA_name, position, identity = nucleotide_string.strip().split(' ')
            position = int(position)
            identity = identity.upper().replace('U', 'T')
            #print nucleotide_string, identity
            assert identity in 'ATCGU'
            if identity in nucleotides_to_count:
                nucleotide_match = self.get_nucleotide(rRNA_name, position)
                assert nucleotide_match.identity == identity
                if not (exclude_constitutive and nucleotide_match.exclude_constitutive):
                    nucleotides.append(nucleotide_match)
        return nucleotides

    def get_all_nucleotides(self, nucleotides_to_count = 'ATCG', exclude_constitutive=False):
        """
        return a list of all nucleotides subject to the optional parameters
        :param nucleotides_to_count:
        :param exclude_constitutive:
        :return:
        """
        nucleotides = []
        for rRNA in self.rRNA_mutation_data.values():
            nucleotides += rRNA.get_all_nucleotides(nucleotides_to_count = nucleotides_to_count,
                                                    exclude_constitutive=exclude_constitutive)
        return nucleotides

    def lowess_correct_mutation_fold_changes(self, nucleotides_to_count ='ATCG', exclude_constitutive=False, max_fold_reduction=0.001, max_fold_increase=100):
        """
        Add a lowess regression corrected fold change, by regressing on the
        :param nucleotides_to_count:
        :param exclude_constitutive:
        :return:
        """
        nucs = self.get_all_nucleotides(nucleotides_to_count = nucleotides_to_count,
                                   exclude_constitutive=exclude_constitutive)
        nuc_fold_changes = [nuc.get_control_fold_change_in_mutation_rate() for nuc in nucs]
        for i in range(len(nuc_fold_changes)):
            if nuc_fold_changes[i] == 0:
                nuc_fold_changes[i] = max_fold_reduction
            elif nuc_fold_changes[i] == float('inf'):
                nuc_fold_changes[i] = max_fold_increase
        nuc_avg_counts = [(nuc.total_mutation_counts+nuc.get_control_nucleotide().total_mutation_counts)/2.0
                          for nuc in nucs]
        fc_log = [math.log(fc, 10) for fc in nuc_fold_changes]
        mag_log = [math.log(m, 10) if m > 0 else -1. for m in nuc_avg_counts]
        lowess_fc_log = lowess(fc_log, mag_log, return_sorted=False)
        lowess_fc = 10 ** lowess_fc_log


        for i in range(len(nucs)):
            nucs[i].lowess_fc = nucs[i].get_control_fold_change_in_mutation_rate()/lowess_fc[i]

class rRNA_mutations:
    def __init__(self, lib, lib_settings, experiment_settings, mutation_filename, rRNA_name):
        self.lib = lib
        self.lib_settings = lib_settings
        self.experiment_settings = experiment_settings
        self.nucleotides = {}
        self.rRNA_name = rRNA_name
        self.parse_mutations_columns(mutation_filename)
        self.total_rt_stops = 0.0

    def parse_mutations_columns(self, filename):
        f= open(filename, 'rU')
        lines =  f.readlines()
        self.sequence = self.experiment_settings.rRNA_seqs[self.rRNA_name]
        headers = lines[0].strip().split('\t')
        position = 1
        for line in lines[1:]:
            if line.strip().strip('\t') != '':
                nucleotide_data = Nucleotide(self, position, headers, line, self.lib_settings)
                self.nucleotides[nucleotide_data.position] = nucleotide_data
                position += 1
        f.close()

    def count_mutation_rates_by_nucleotide(self, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        """
        counts, over this RNA, the total number of mutations at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.

        NOTE that this will set any background-subtracted rate of less than zero to zero

        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        counts = defaultdict(int)
        for nucleotide in self.nucleotides.values():
            if exclude_constitutive and nucleotide.exclude_constitutive:
                continue
            else:
                if subtract_background and subtract_control:
                    raise SyntaxError('Cannot subtract background and control simultaneously')

                if subtract_background:
                    counts[nucleotide.identity] += max((nucleotide.mutation_rate - self.lib.get_normalizing_lib().
                                                    get_mutation_rate_at_position(self.rRNA_name, nucleotide.position)), 0.)
                elif subtract_control:
                    counts[nucleotide.identity] += nucleotide.mutation_rate - self.lib.get_normalizing_lib_with_mod().get_mutation_rate_at_position(self.rRNA_name, nucleotide.position)
                else:
                    counts[nucleotide.identity] += nucleotide.mutation_rate
        return counts

    def count_rt_stop_rpm_by_nucleotide(self, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        """
        counts, over this RNA, the total number of mutations at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a modification.

        NOTE that this will set any background-subtracted rate of less than zero to zero

        :return: a dict like {A: 1054, T:32, C: 604, G:99}
        """
        counts = defaultdict(int)
        for nucleotide in self.nucleotides.values():
            if exclude_constitutive and nucleotide.exclude_constitutive:
                continue
            else:
                if subtract_background and subtract_control:
                    raise SyntaxError('Cannot subtract background and control simultaneously')

                if subtract_background:
                    counts[nucleotide.identity] += max((nucleotide.get_rt_stop_rpm() - self.lib.get_normalizing_lib().
                                                    get_rt_stop_rpm_at_position(self.rRNA_name, nucleotide.position)), 0.)
                elif subtract_control:
                    counts[nucleotide.identity] += nucleotide.get_rt_stop_rpm() - self.lib.get_normalizing_lib_with_mod().get_rt_stop_rpm_at_position(self.rRNA_name, nucleotide.position)
                else:
                    counts[nucleotide.identity] += nucleotide.get_rt_stop_rpm()
        return counts

    def count_mutation_types_by_nucleotide(self, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        """
        counts, over this RNA, the total number of mutation of each type at each of A, T, C, G
        This is to get an idea of which nucleotides are being affected by a particular mutation

        NOTE that this will set any background-subtracted rate of less than zero to zero
        """
        counts = defaultdict((lambda : defaultdict(int)))
        for nucleotide in self.nucleotides.values():
            if exclude_constitutive and nucleotide.exclude_constitutive:
                continue
            else:
                for mutation_type in nucleotide.mutations_by_type:
                    counts[nucleotide.identity][mutation_type] += nucleotide.mutations_by_type[mutation_type]
        return counts

    def list_mutation_rates(self, subtract_background=False, subtract_control = False, nucleotides_to_count='ATCG', exclude_constitutive=False):
        """
        #note that these values may be less than zero when background is subtracted
        :param subtract_background:
        :return:
        """
        rates = []
        for nucleotide in self.nucleotides.values():
            if nucleotide.identity in nucleotides_to_count:
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    continue
                else:
                    if subtract_background and subtract_control:
                        raise SyntaxError('Cannot subtract background and control simultaneously')

                    if subtract_background:

                        rates.append((nucleotide.mutation_rate - self.lib.get_normalizing_lib().
                                                        get_mutation_rate_at_position(self.rRNA_name, nucleotide.position)))
                    elif subtract_control:

                        rates.append((nucleotide.mutation_rate - self.lib.get_normalizing_lib_with_mod().
                                                        get_mutation_rate_at_position(self.rRNA_name, nucleotide.position)))
                    else:
                        rates.append(nucleotide.mutation_rate)
        return rates

    def list_rt_stop_rpms(self, subtract_background=False, subtract_control = False, nucleotides_to_count='ATCG', exclude_constitutive=False):
        """
        #note that these values may be less than zero when background is subtracted
        :param subtract_background:
        :return:
        """
        rates = []
        for nucleotide in self.nucleotides.values():
            if nucleotide.identity in nucleotides_to_count:
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    continue
                else:
                    if subtract_background and subtract_control:
                        raise SyntaxError('Cannot subtract background and control simultaneously')

                    if subtract_background:

                        rates.append((nucleotide.get_rt_stop_rpm() - self.lib.get_normalizing_lib().
                                                        get_rt_stop_rpm_at_position(self.rRNA_name, nucleotide.position)))
                    elif subtract_control:

                        rates.append((nucleotide.get_rt_stop_rpm() - self.lib.get_normalizing_lib_with_mod().
                                                        get_rt_stop_rpm_at_position(self.rRNA_name, nucleotide.position)))
                    else:
                        rates.append(nucleotide.get_rt_stop_rpm())
        return rates

    def list_rt_stop_counts(self, nucleotides_to_count='ATCG'):
        """
        #note that these values may be less than zero when background is subtracted
        :param subtract_background:
        :return:
        """
        rates = []
        for nucleotide in self.nucleotides.values():
            if nucleotide.identity in nucleotides_to_count:
                rates.append(nucleotide.rt_stops)
        return rates


    def list_mutation_fold_changes(self, nucleotides_to_count='ATCG', exclude_constitutive=False):
        """
        #note that these values may be less than zero when background is subtracted
        :param subtract_background:
        :return:
        """
        rates = []
        for nucleotide in self.nucleotides.values():
            if nucleotide.identity in nucleotides_to_count:
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    continue
                elif nucleotide.get_control_fold_change_in_mutation_rate() == 0.0 or \
                                nucleotide.get_control_fold_change_in_mutation_rate() == float('inf'):
                    continue
                else:
                    rates.append(nucleotide.get_control_fold_change_in_mutation_rate())
        return rates

    def get_changed_nucleotides(self, change_type, nucleotides_to_count='ATCG', exclude_constitutive=False,
                                  confidence_interval = 0.99, fold_change_cutoff = 3, subtract_background=False):
        nucleotides = []
        for nucleotide in self.nucleotides.values():
            if nucleotide.identity in nucleotides_to_count:
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    continue
                else:
                    prot_call = nucleotide.determine_protection_status(confidence_interval=confidence_interval,
                                                           fold_change_cutoff=fold_change_cutoff,
                                                           subtract_background=subtract_background)
                    if prot_call == change_type:
                        nucleotides.append(nucleotide)
        return nucleotides


    def get_all_nucleotides(self, nucleotides_to_count='ATCG', exclude_constitutive=False):
        """
        return a list of all nucleotides subject to the optional parameters
        :param nucleotides_to_count:
        :param exclude_constitutive:
        :return:
        """
        nucleotides = []
        for nucleotide in self.nucleotides.values():
            if nucleotide.identity in nucleotides_to_count:
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    continue
                else:
                    nucleotides.append(nucleotide)
        return nucleotides


class Nucleotide:
    def __init__(self, rRNA, position, headers, mutation_data_line, lib_settings):
        self.rRNA = rRNA
        self.position = position
        self.identity = rRNA.sequence[self.position-1]
        self.mutations_by_type = {} #will map each type of mutation to the number of such mutations detected
        self.rt_stops = 0
        self.lib_settings = lib_settings
        self.parse_mutation_data_line(headers, mutation_data_line)
        self.set_exclusion_flag()



    def __str__(self):
        return "%s%d in %s of %s" % (self.identity, self.position, self.rRNA.rRNA_name, self.lib_settings.sample_name)


    def parse_mutation_data_line(self, headers, mutation_data_line):
        ll = mutation_data_line.strip().split('\t')
        self.total_mutation_counts = sum([float(ll[i]) for i in range(0, 26)])
        self.sequencing_depth = float(ll[26])
        self.effective_sequencing_depth = float(ll[27])
        try:
            self.mutation_rate = self.total_mutation_counts/self.sequencing_depth
        except:
            self.mutation_rate = 0      
        for i in range(0, 26):
            self.mutations_by_type[headers[i]] = float(ll[i])

    def set_exclusion_flag(self):
            try:
                exclusions = self.lib_settings.experiment_settings.exclude_constitutive[self.rRNA.rRNA_name]
                if self.position in exclusions:
                    self.exclude_constitutive = True
                else:
                    self.exclude_constitutive = False
            except KeyError:
                self.exclude_constitutive = False

    def get_back_sub_mutation_rate(self):
        return (self.mutation_rate - self.get_background_nucleotide().mutation_rate)

    def get_control_sub_mutation_rate(self, subtract_background=False):
        if subtract_background:
            return (self.get_back_sub_mutation_rate() - self.get_control_nucleotide().get_back_sub_mutation_rate())
        else:
            return (self.mutation_rate - self.get_control_nucleotide().mutation_rate)

    def get_control_fold_change_in_mutation_rate(self, subtract_background = False):
        try:
            if subtract_background:
                return self.get_back_sub_mutation_rate()/self.get_control_nucleotide().get_back_sub_mutation_rate()
            else:
                return self.mutation_rate/self.get_control_nucleotide().mutation_rate
        except ZeroDivisionError:
            return float('inf')

    def get_fold_signal_over_background(self, background_nuc):
        try:
            return (self.mutation_rate/background_nuc.mutation_rate)
        except ZeroDivisionError:
            return float('inf')

    def get_control_fold_change_error(self, subtract_background=False, max_fold_reduction=0.001, max_fold_increase=100,
                                      lowess_correct=False):
        try:
            if lowess_correct:
                ratio = self.lowess_fc
            else:
                ratio = self.get_control_fold_change_in_mutation_rate(subtract_background=subtract_background)
            if ratio == float('inf') or ratio == -1*float('inf'):
                ratio = max_fold_increase
            elif ratio<=0:
                ratio = max_fold_reduction
            if subtract_background:
                num = self.get_back_sub_mutation_rate()
                num_error = self.get_back_sub_error()
                denom = self.get_control_nucleotide().get_back_sub_mutation_rate()
                denom_error = self.get_control_nucleotide().get_back_sub_error()
            else:
                num = self.mutation_rate
                num_error = self.get_error()
                denom = self.get_control_nucleotide().mutation_rate
                denom_error = self.get_control_nucleotide().get_error()
            return ratio*math.sqrt((num_error/num)**2+(denom_error/denom)**2)
        except ZeroDivisionError:
            return float('inf')

    def get_signal_error(self, background_nuc, max_fold_reduction=0.0001, max_fold_increase=10000):
        '''
        return the counting error for the signal of test_nuc over background_nuc
        '''

        try:
            ratio = self.get_fold_signal_over_background(background_nuc)
            if ratio == float('inf') or ratio == -1*float('inf'):
                ratio = max_fold_increase
            elif ratio<=0:
                ratio = max_fold_reduction
            num = self.mutation_rate
            num_error = self.get_error()
            denom = background_nuc.mutation_rate
            denom_error = background_nuc.get_error()
            return ratio*math.sqrt((num_error/num)**2+(denom_error/denom)**2)
        except ZeroDivisionError:
            return float('inf')

    def get_control_mutation_rate(self):
            return self.rRNA.lib.get_normalizing_lib_with_mod().\
                get_mutation_rate_at_position(self.rRNA.rRNA_name, self.position)

    def get_control_nucleotide(self):
        return self.rRNA.lib.get_normalizing_lib_with_mod().rRNA_mutation_data[self.rRNA.rRNA_name].nucleotides[self.position]

    def get_background_nucleotide(self):
        return self.rRNA.lib.get_normalizing_lib().rRNA_mutation_data[self.rRNA.rRNA_name].nucleotides[self.position]

    def get_wilson_approximate_score_interval(self, confidence_interval = 0.99):
        """
        Computes the wilson score interval, which APPROXIMATES the confidence interval for the mean of the binomial
        distribution, given a sampling of the distribution.
        :return:
        """
        alpha = (1.0-confidence_interval)
        z = 1.0-(alpha/2.0)
        n = self.sequencing_depth
        p = self.mutation_rate
        #breaking up equation
        a = 1.0/(1.0+(z**2)/n)
        b = p+((z**2)/(2.0*n))
        c = z*math.sqrt((p*(1.0-p))/n + (z**2)/(4*(n**2)))
        interval_bottom = a*(b-c)
        interval_top = a*(b+c)
        return interval_bottom, interval_top

    def determine_protection_status(self, confidence_interval = 0.99, fold_change_cutoff = 5, subtract_background=False,
                                    max_fold_reduction=0.001, max_fold_increase=100, lowess_correct=False):
        #self_min, self_max = self.get_wilson_approximate_score_interval(confidence_interval=confidence_interval)
        #control_min, control_max = self.get_control_nucleotide().\
        #    get_wilson_approximate_score_interval(confidence_interval=confidence_interval)
        #if mod_utils.ranges_overlap(self_min, self_max, control_min, control_max) \
        #        or (self.get_control_fold_change_in_mutation_rate()<fold_change_cutoff
        #            and self.get_control_fold_change_in_mutation_rate()>1.0/fold_change_cutoff) or self.identity not in \
        #        self.lib_settings.experiment_settings.get_property('affected_nucleotides'):
        #    return "no_change"
        if lowess_correct:
            fold_change = self.lowess_fc
            if subtract_background:
                print "WARNING: lowess correction overrides background subtraction"
        else:
            fold_change = self.get_control_fold_change_in_mutation_rate(subtract_background=subtract_background)
        #these outliers are always on the edge of the rRNA, so they're probably spurious low-coverage events
        if not (self.signal_above_background(self.get_control_nucleotide(), self.get_background_nucleotide(), confidence_interval=0.9) or
                self.signal_above_background(self, self.get_background_nucleotide(), confidence_interval=0.9)):
            return "no_change"
        elif fold_change == float('inf') or fold_change == -1*float('inf'):
            #fold_change = max_fold_increase
            return "no_change"
        elif fold_change<=0:
            #fold_change = max_fold_reduction
            return "no_change"

        mean = math.log(fold_change) #natural log to make dist more gaussian
        standard_deviation = self.get_control_fold_change_error(subtract_background=subtract_background,
                                                                lowess_correct=lowess_correct)/fold_change #error propogation for natural log
        p, z = mod_utils.computePfromMeanAndStDevZscore(mean, standard_deviation, 0) #what is the chance that no change could come from this dist?
        if (p > 1.0-confidence_interval and p<confidence_interval)or (self.get_control_fold_change_in_mutation_rate(subtract_background=subtract_background)<fold_change_cutoff
                                             and self.get_control_fold_change_in_mutation_rate(subtract_background=subtract_background)>1.0/fold_change_cutoff)\
                or self.identity not in self.lib_settings.experiment_settings.get_property('affected_nucleotides'):
            return "no_change"
        elif self.get_control_sub_mutation_rate(subtract_background=subtract_background)<0:
            return "protected"
        elif self.get_control_sub_mutation_rate(subtract_background=subtract_background)>0:
            return "deprotected"
        else:
            return "something_is_wrong_change_zero"

    def signal_above_background(self, test_nuc, background_nuc, confidence_interval = 0.9, max_fold_reduction=0.0001, max_fold_increase=10000):
        '''
        return True if signal in test dataset is statistically significantly above the  background dataset
        Must provide a nucleotide object for test_lib and background_lib
        '''
        fold_change = test_nuc.get_fold_signal_over_background(background_nuc)
        if fold_change == float('inf') or fold_change == -1 * float('inf'):
            fold_change = max_fold_increase
        elif fold_change <= 0:
            fold_change = max_fold_reduction
        mean = math.log(fold_change) #natural log to make dist more gaussian
        standard_deviation = test_nuc.get_signal_error(background_nuc)/fold_change #error propogation for natural log
        p, z = mod_utils.computePfromMeanAndStDevZscore(mean, standard_deviation, 0) #what is the chance that no change could come from this dist?
        if (p<confidence_interval):
            return False
        elif test_nuc.mutation_rate-background_nuc.mutation_rate<0:
            return False
        elif test_nuc.mutation_rate-background_nuc.mutation_rate>0:
            return True
        else:
            return False


    def get_error(self):
        try:
            return(np.sqrt(self.mutation_rate/self.sequencing_depth))
        except ZeroDivisionError:
            return float('inf')

    def get_back_sub_error(self):
        mutation_rate = self.get_back_sub_mutation_rate()
        if mutation_rate < 0:
            mutation_rate = 0
        try:
            return(np.sqrt(mutation_rate/self.sequencing_depth))
        except ZeroDivisionError:
            return float('inf')

    def get_control_sub_error(self):
        mutation_rate = self.get_back_sub_mutation_rate()
        if mutation_rate < 0:
            mutation_rate = 0
        try:
            return(np.sqrt(mutation_rate/self.sequencing_depth))
        except ZeroDivisionError:
            return float('inf')

    def get_rt_stop_rpm(self):
        return self.rt_stops/(self.rRNA.total_rt_stops/1000000.0)

