from collections import defaultdict
import matplotlib.pyplot as plt
import re
import scipy.stats
import subprocess
import os
import cPickle
import mod_utils
import gzip
import numpy as np
from collections import Counter
import math

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

    def list_mutation_rates(self, subtract_background = False, subtract_control = False, nucleotides_to_count = 'ATCG', exclude_constitutive=False):
        all_mutation_rates = []
        for rRNA_name in self.rRNA_mutation_data:
            all_mutation_rates.extend(self.rRNA_mutation_data[rRNA_name].
                                      list_mutation_rates(subtract_background = subtract_background, subtract_control = subtract_control,
                                                          nucleotides_to_count = nucleotides_to_count, exclude_constitutive=exclude_constitutive))
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

    def write_tsv_tables(self, tsv_filename, subtract_background=False, subtract_control=False, exclude_constitutive=False):
        if subtract_background and subtract_control:
            raise SyntaxError('Cannot subtract background and control simultaneously')

        f = open(tsv_filename, 'w')

        if subtract_background:
                f.write('CHROMOSOME\tPOSITION\tMUTATION_RATE\tBKGD_SUB_MUT_RATE\tBKGD_SUB_ERROR\n')
        elif subtract_control:
                f.write('CHROMOSOME\tPOSITION\tNUC\tEXP_MUTATION_RATE\tEXP_99%_min\tEXP_99%_max\tCTRL_MUT_RATE'
                        '\tCTRL_99%_min\tCTRL_99%_max\tEXP-CTRL\tCTRL_POISSON_SUB_ERROR\tFOLD_CHANGE\tPROTECTION_CALL\n')
        elif not subtract_background and not subtract_control:
                f.write('CHROMOSOME\tPOSITION\tMUTATION_RATE\tERROR\n')


        for rRNA_name in self.rRNA_mutation_data:
            for position in self.rRNA_mutation_data[rRNA_name].nucleotides:
                nucleotide = self.rRNA_mutation_data[rRNA_name].nucleotides[position]
                if exclude_constitutive and nucleotide.exclude_constitutive:
                    if subtract_background:
                        f.write(self.rRNA_mutation_data[rRNA_name].rRNA_name+'\t'+str(nucleotide.position)+'\t'
                                +'0'+'\t'+'0'+'\t'
                                +'0'+'\n')
                    elif subtract_control:
                        f.write(self.rRNA_mutation_data[rRNA_name].rRNA_name+'\t'+str(nucleotide.position)+
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
                        exp_wil_bottom, exp_wil_top = nucleotide.get_wilson_approximate_score_interval()
                        ctrl_wil_bottom, ctrl_wil_top = ctrl_nuc.get_wilson_approximate_score_interval()
                        f.write('%s\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n' %
                                (rRNA_name, nucleotide.position, nucleotide.identity, nucleotide.mutation_rate,
                                exp_wil_bottom, exp_wil_top, ctrl_nuc.mutation_rate,
                                ctrl_wil_bottom, ctrl_wil_top, nucleotide.get_control_sub_mutation_rate(),
                                nucleotide.get_control_sub_error(), nucleotide.get_control_fold_change_in_mutation_rate(),
                                nucleotide.determine_protection_status(confidence_interval=self.experiment_settings.get_property('confidence_interval_cutoff'),
                                                                       fold_change_cutoff=self.experiment_settings.get_property('fold_change_cutoff'))))
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
                nucleotide_data = Nucleotide(self, headers, line, self.lib_settings)
                self.nucleotides[nucleotide_data.position] = nucleotide_data
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
                pass
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
                    pass
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

class Nucleotide:
    def __init__(self, rRNA, headers, mutation_data_line, lib_settings):
        self.rRNA = rRNA
        self.mutations_by_type = {} #will map each type of mutation to the number of such mutations detected
        self.lib_settings = lib_settings
        self.parse_mutation_data_line(headers, mutation_data_line)
        self.set_exclusion_flag()


    def __str__(self):
        return "%s%d in %s of %s" % (self.identity, self.position, self.rRNA.rRNA_name, self.lib_settings.sample_name)


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
        return (self.mutation_rate - self.rRNA.lib.get_normalizing_lib().\
            get_mutation_rate_at_position(self.rRNA.rRNA_name, self.position))

    def get_control_sub_mutation_rate(self):
        return (self.mutation_rate - self.rRNA.lib.get_normalizing_lib_with_mod().\
            get_mutation_rate_at_position(self.rRNA.rRNA_name, self.position))

    def get_control_fold_change_in_mutation_rate(self):
        try:
            return (self.mutation_rate/self.rRNA.lib.get_normalizing_lib_with_mod().\
                get_mutation_rate_at_position(self.rRNA.rRNA_name, self.position))
        except ZeroDivisionError:
            return float('inf')
    def get_control_mutation_rate(self):
            return self.rRNA.lib.get_normalizing_lib_with_mod().\
                get_mutation_rate_at_position(self.rRNA.rRNA_name, self.position)

    def get_control_nucleotide(self):
        return self.rRNA.lib.get_normalizing_lib_with_mod().rRNA_mutation_data[self.rRNA.rRNA_name].nucleotides[self.position]

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

    def determine_protection_status(self, confidence_interval = 0.99, fold_change_cutoff = 5):
        self_min, self_max = self.get_wilson_approximate_score_interval(confidence_interval=confidence_interval)
        control_min, control_max = self.get_control_nucleotide().\
            get_wilson_approximate_score_interval(confidence_interval=confidence_interval)
        if mod_utils.ranges_overlap(self_min, self_max, control_min, control_max) \
                or (self.get_control_fold_change_in_mutation_rate()<fold_change_cutoff
                    and self.get_control_fold_change_in_mutation_rate()>1.0/fold_change_cutoff) or self.identity not in \
                self.lib_settings.experiment_settings.get_property('affected_nucleotides'):
            return "no_change"
        elif self.get_control_sub_mutation_rate()<0:
            return "protected"
        elif self.get_control_sub_mutation_rate()>0:
            return "deprotected"
        else:
            return "something_is_wrong_change_zero"

    def get_error(self):
        try:
            return(np.sqrt(self.mutation_rate/self.sequencing_depth))
        except ZeroDivisionError:
            return(0)

    def get_back_sub_error(self):
        mutation_rate = self.get_back_sub_mutation_rate()
        if mutation_rate < 0:
            mutation_rate = 0
        try:
            return(np.sqrt(mutation_rate/self.sequencing_depth))
        except ZeroDivisionError:
            return(0)

    def get_control_sub_error(self):
        mutation_rate = self.get_back_sub_mutation_rate()
        if mutation_rate < 0:
            mutation_rate = 0
        try:
            return(np.sqrt(mutation_rate/self.sequencing_depth))
        except ZeroDivisionError:
            return(0)

