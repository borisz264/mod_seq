__author__ = 'boris'
"""
Based on the Rouskin DMS-seq paper:
    True Positives: Bases that are unpaired in the secondary structure, and the reactive atom has a solvent accessible
        surface area (to a 3A radius sphere) of greater than 2A squared.
    True Negatives: are Watson-crick paired (A-U or C-G) in the secondary structure model

"""

import sys, mod_utils
from rna import rna

class rRNA:
    def __init__(self, bpseq_filename):
        self.nucleotides = {}
        self.parse_bpseq(bpseq_filename)


    def parse_bpseq(self, bpseq_filename):
        """

        :param bpseq_filename:
        :return: a dictionary mapping position
        """
        f = open(bpseq_filename, 'rU')
        lines = f.readlines()
        for line in lines[4:]:
            nt_num, nt_ident, paired_to = line.rstrip().split(' ')
            self.nucleotides[int(nt_num)] = nucleotide(self, nt_num, nt_ident, paired_to)
        f.close()

    def parse_rouskin_accessibility(self, filenames):
        for filename in filenames:
            f = open(filename, 'rU')
            for line in f:
                ll = line.rstrip().split('\t')
                self.nucleotides[int(ll[0])].sasa = float(ll[1])
            f.close()

    def parse_zinshteyn_accessibility(self, filename):
        f = open(filename)
        for line in f:
            ll = line.rstrip().split('\t')
            assert rna(ll[1]) == self.nucleotides[int(ll[0])].identity
            self.nucleotides[int(ll[0])].sasa = float(ll[2])
        f.close()

    def write_to_file(self, sasa_cutoff, nucs_to_test, output_prefix):
        f = open('%s_%.1f_%s.txt' % (output_prefix, sasa_cutoff, nucs_to_test), 'w')

        for nt_pos in self.nucleotides:
            nuc = self.nucleotides[nt_pos]
            assert nt_pos == nuc.position
            assert not (nuc.is_true_positive(sasa_cutoff) and nuc.is_true_negative())
            if nuc.identity in nucs_to_test:
                if nuc.is_true_positive(sasa_cutoff):
                    result_string = 'TP'
                elif nuc.is_true_negative():
                    result_string = 'TN'
                else:
                    result_string = 'X'
            else:
                result_string = 'X'
            f.write('%d\t%s\t%s\n' % (nuc.position, nuc.identity, result_string))

class nucleotide:
    def __init__(self, rRNA, position, identity, paired_to):
        self.rRNA = rRNA
        self.identity = identity
        self.position = int(position)
        self.paired_to = int(paired_to)
        self.sasa = None #solvent accessible surface ares

    def is_WC_paired(self):
        if not self.paired_to == 0:
            if mod_utils.revComp(self.identity, isRNA = True) == self.rRNA.nucleotides[self.paired_to].identity:
                #if this nucleotide is the WC complement of the one it's paired to, return True
                return True
        return False

    def is_true_positive(self, sasa_cutoff):
        if (not self.is_WC_paired()) and (not self.sasa == None) and self.sasa > sasa_cutoff:
            return True
        else:
            return False

    def is_true_negative(self):
        if self.is_WC_paired():
            return True
        else:
            return False

def main():
    bpseq_filename, outprefix = sys.argv[1:3]
    #rouskin_sasa_files = sys.argv[3:]
    zin_sasa_file = sys.argv[3]
    rRNA_data = rRNA(bpseq_filename)
    #rRNA_data.parse_rouskin_accessibility(rouskin_sasa_files)
    rRNA_data.parse_zinshteyn_accessibility(zin_sasa_file)
    rRNA_data.write_to_file(2.0, 'AC', outprefix)

if __name__ == '__main__':
    main()