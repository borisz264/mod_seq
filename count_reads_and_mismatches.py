#Author: MKT
#purpose: to parse bowtie output files (after they have been weighted based on multiple mappings with weight_repeats.py and convert to pickled python dictionary. This one hashes by the start position of the read so it will be fast to classify them into exons, introns, etc. with the annotation file.
#version history:
#created 9/?/09
#updated 1/26/10- changed it to hash by start with respect to strand and not the middle of the read. For - strand reads this means adding 34 to the start w.r.t the + strand for a 35 bp read.
#4/1/10	Pavan: commented out superfluous print commands
#5/26/10 Boris: adjusting antisense read position by 20 instead of 34, since the mapping only uses 21 bp, this needs to be changed if the mapping size is changed. If read sizes aren't constant, then OH BOY!!! Also, added M1killer chromosome
#6/24/10 Boris: Corrected for bowtie vs genome indexing
#6/23/2013 Boris: Removed arbitrary hard-coding chromosome names, and accounted for split reads (indels, splice junctions)
#7/21/2015 Thomas: Cleaned up code a bit, in addition to the mapping positions of the read 5p ends, it now outputs in the same format 1) number of times each position was covered by a read, 2) number of times each position was mutated relatove to reference sequence
strands= ['+','-']

import sys, cPickle as pickle, gzip, re, numpy, os
import mod_settings, mod_utils
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
from collections import defaultdict

def createStrandDict(strands):
    """
    """
    d = {}
    for strand in strands:
        d[strand] = {}
    return d

def pickleDict(d, fn, suffix):
    
    fn = fn.split('.')
    out_fn = ''.join([fn[0]+suffix+'.pkl'])
    with open(out_fn, 'w') as g:
        pickle.dump(d, g)

def readGenomicCoverage(rel_cov, strand, read_start):
    """
    Assumes:
        rel_cov, a list relating position in read (i) to relative genomic postion rel_cov[i]
        strand, the mapping strand of the read
        read_start, the starting genomic position of the read
    Returns:
        pos_l, a list relating position in read (i) to absolute genomic position pos_l[i]
    """
    
    
    if strand == '+':
        pos_l = [p+read_start for p in rel_cov]
    elif strand == '-':
        pos_l = [read_start-p for p in rel_cov]
    return pos_l

def parse_MDz_and_cigar(cigarString, MDzString, mappingLength, seq):
    """
    Assumes:
        cigarString - a string of [#][A-Z][#][A-Z] etc... that described the alignment of a read to the genome
        MDzString, the MD:z string for the current read
        genome_cov, a list relating position in read (i) to absolute genomic position pos_l[i]
        seq, a string representation of the read
    Returns:
        A list of positions with mismatches in absolute genomic position
        the genomic distance spanned by this read


    CIGAR STRING:
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch
    H can only be present as the
    rst and/or last operation.
    S may only have H operations between them and the ends of the CIGAR string.
    For mRNA-to-genome alignment, an N operation represents an intron. For other types of
    alignments, the interpretation of N is not de
    ned.
    4 Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
    """

    #Parse cigar strings into numbers & tags
    numbers = re.split('[A-Z,=]', cigarString)[:-1] #not sure why, but this seems to produce an extra blank entry thay I'm stripping off
    tags = re.split('[0-9]*', cigarString)[1:]
    numbers = [int(n) for n in numbers]

    assert len(numbers) == len(tags)

    genomeMulitpliers = {'M':1, 'I':0, 'D':1, 'N':1, '=':1, 'X':1}
    readMulitpliers = {'M':1, 'I':1, 'D':0, 'N':0, '=':1, 'X':1}

    # Initialize counters
    readMappingSpan, genomeMappingSpan, genome_pos = 0, 0, 0
    genomeCoverage = []

    #print 'New Read'
    insertions = {}# need to know insertion positions and sizes in read for later parsing
    cigar_tups = zip(tags, numbers)
    for cigar in cigar_tups:
        tag, tag_len = cigar
        if tag == 'I':
            insertions[readMappingSpan] = tag_len
        genome_multiplier = genomeMulitpliers[tag]
        read_multiplier = readMulitpliers[tag]
        readMappingSpan += read_multiplier*tag_len
        genomeMappingSpan += genome_multiplier*tag_len

        # Determines which  genomic positions relative to the read 5' end are covered by a read
        # For M, =, X: positions covered by read: so add to list, and increment counter
        # For D (deletion) and N (skipped region, similar to deletion, meant for introns): positions covered by read:
        #   so add to list, and increment counter, so the length of the array
        #   can be longer than the length of the read, but it should match the genome mapping span
        #   I include N, since for my application of mapping to rRNA, anything defined as N should be an RT deletion that
        #   I consider spanned by the read. If there are true intronic reads, then the counter ought be incremented without
        #   adding to the array.
        if genome_multiplier: # Increment counter if genome multiplier = 1, ie position spanned by read
            if tag in 'M=XDN':
                genomeCoverage += range(genome_pos,genome_pos+tag_len)
                genome_pos += tag_len

    assert readMappingSpan == mappingLength
    assert genomeMappingSpan == len(genomeCoverage)

    """
        test all cases at:
    http://davetang.org/muse/2011/01/28/perl-and-sam/

    The MD field aims to achieve SNP/indel calling without looking at the reference.
    For example, a string `10A5^AC6' means from the leftmost reference base in the
    alignment, there are 10 matches followed by an A on the reference which is
    different from the aligned read base; the next 5 reference bases are matches
    followed by a 2bp deletion from the reference; the deleted sequence is AC;
    the last 6 bases are matches. The MD field ought to match the CIGAR string.
    The string '0T0C37T' indicates the first base of the read is a mismatch from
     the reference sequence of T, the second base is a mismatch from the reference
     sequence C, followed by 37 matches, with a final mismatch from the reference A

    Boris 20151116 - this is getting an overhaul to deal with indels tags
    THomas's `10A5^AC6' example above should identify a mismatch at position 10 (zero-indexed), a deletion of length 2
     at position 16
    all indels will be assigned to the 3' end of the event. Only deletions will be used for the final count

    """



    MDz = MDzString.split(':')[2]
    MD_tags = re.findall('[0-9]+|[A-Z,^]+',MDz)
    MD_tags = [int(x) if re.match('[0-9]+',x) else x for x in MD_tags]

    mutations_rel_read = ['M']*readMappingSpan #will keep track of mutations along the read, first nuc of read is position 0
    mutations_rel_genome = ['M']*genomeMappingSpan #will keep track of mutations along the genome, first nuc of read is position 0
    genomic_event_positions = [] # A list to store positions of mismatches
    read_counter = 0 # A counter to store the current position of the read (0 indexed)
    genome_counter = 0 # A counter to store the current position along the genome (0 indexed with respect to read start)
    for tag in MD_tags:
        #insertion relative to genome?
        if read_counter in insertions:
            assert insertions[read_counter] != 0
            insertion_size = insertions[read_counter]
            read_counter += insertions[read_counter]
            assert mutations_rel_read[read_counter-1] =='M'
            mutations_rel_read[read_counter-1] = ('I', insertion_size)

        if tag in ['A','T','C','G']: #If the tag is a base it indicates a mismatch
            if seq[read_counter] != 'N': # We want to ignore 'N' bases in sequence
                assert mutations_rel_read[read_counter] =='M'
                assert mutations_rel_genome[genome_counter] =='M'
                mutations_rel_read[read_counter] = (tag, seq[read_counter])
                mutations_rel_genome[genome_counter] = (tag, seq[read_counter])
                genomic_event_positions.append(read_counter) # Append the current position in the read
            read_counter += 1 #Increment the counter
            genome_counter += 1
        elif isinstance(tag, int): # If the tag is an int, this represents matches
            genome_counter += tag
            read_counter += tag # Increment counter by number of matches
        elif tag.startswith('^'):#the read contains a deletion relative to the genome
            deletion_length = len(tag)-1
            genome_counter += deletion_length
            assert mutations_rel_genome[genome_counter-1] =='M'
            mutations_rel_genome[genome_counter-1] = ('D', deletion_length)
            genomic_event_positions.append(genome_counter-1)

    #print mutations_rel_read
    #print mutations_rel_genome
    #print read_counter, genome_counter
    #print genomeCoverage
    return genomic_event_positions, genomeCoverage, mutations_rel_genome, mutations_rel_read, readMappingSpan, genomeMappingSpan

def checkTag(tag, fields):
    """
    Assumes:
        tag, a str specifying a SAMtools tag
        fields, a SAMfile line in list format
    Does:
        Determines if the tag is at the expected index, if not, searches for the tag in fields
    Returns:
        fields[i], containing the tag
    """
    tag_fields = [field for field in fields[7:] if tag in field]
    if len(tag_fields) == 0:#uh oh, no tag found
        if tag == 'NH:i:': #my alignments seem to lack this tag
            return 'NH:i:1'
    return tag_fields[0]

def pie_read_5p_ends(read_5p_ends, genome_dict, out_prefix):

    fig = plt.figure(figsize=(8,17))

    plot = fig.add_subplot(211)
    nuc_counts = defaultdict(int)
    for chromosome in read_5p_ends['+']:
        for position in read_5p_ends['+'][chromosome]:
            if position-1 > 0 :
                nuc = genome_dict[chromosome][position-1]
                nuc_counts[nuc] += read_5p_ends['+'][chromosome][position]
    labels = sorted(nuc_counts.keys())
    sizes = [nuc_counts[nt] for nt in labels]
    plot.pie(sizes, labels = labels, colors = mod_utils.rainbow)
    plot.set_title('nt exactly at read 5p ends across rRNA')

    plot = fig.add_subplot(212)
    nuc_counts = defaultdict(int)
    for chromosome in read_5p_ends['+']:
        for position in read_5p_ends['+'][chromosome]:
            if position-2 > 0 :
                nuc = genome_dict[chromosome][position-2]
                nuc_counts[nuc] += read_5p_ends['+'][chromosome][position]
    labels = sorted(nuc_counts.keys())
    sizes = [nuc_counts[nt] for nt in labels]
    plot.pie(sizes, labels = labels, colors = mod_utils.rainbow)
    plot.set_title('1nt upstream of read 5p ends across rRNA')

    plt.savefig(out_prefix + '_nt_5p_ends_pie.pdf', transparent='True', format='pdf')
    plt.clf()



def plot_mutated_nts_pie(mutated_nts_count, title, out_prefix):
    fig = plt.figure(figsize=(8,8))

    plot = fig.add_subplot(111)#first a pie chart of mutated nts
    labels = sorted(mutated_nts_count.keys())
    sizes = [mutated_nts_count[nt] for nt in labels]
    total = float(sum(sizes))
    merged_labels = ['%s %.3f' % (labels[i], sizes[i]/total) for i in range(len(sizes))]
    plot.pie(sizes, labels = merged_labels, colors = mod_utils.rainbow)
    plot.set_title(title)

    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def plot_full_mutation_stats(mutations_counts, indel_distribution, mutations_by_position, positional_coverage, title, x_label, out_prefix):
    fig = plt.figure(figsize=(16,16))
    fig.suptitle(title)

    plot = fig.add_subplot(221)#first a pie chart of mutated nts
    labels = sorted(mutations_counts.keys())
    sizes = [mutations_counts[nt] for nt in labels]
    plot.pie(sizes, labels = labels, colors = mod_utils.rainbow)

    plot = fig.add_subplot(222)#second a histogram of indel sizes
    bins = range(1, 20)
    bins.append(100)
    plot.hist(indel_distribution, color = mod_utils.black, bins = bins)
    plot.set_xlim(0,10)
    plot.set_xticks(numpy.arange(0,10)+0.5)
    plot.set_xticklabels(numpy.arange(0,10))
    plot.set_xlabel(x_label)
    plot.set_ylabel("# events")

    bar_width = 0.5
    all_events = sorted(mutations_counts, key=mutations_counts.get, reverse=True)
    top_events = all_events[:len(mod_utils.rainbow)]
    mut_positions = sorted(mutations_by_position.keys())
    cov_positions = sorted(positional_coverage.keys())

    plot = fig.add_subplot(223)#a stacked bar graph of positional mutation rates
    bottoms = [0]*len(mut_positions)
    bottoms = numpy.array(bottoms)
    plot_layers = []
    color_index = 0
    for event in top_events:
        event_amounts = [mutations_by_position[position][event] if event in mutations_by_position[position] else 0 for position in mut_positions]
        plot_layers.append(plot.bar(mut_positions, event_amounts, bar_width, bottom = bottoms, color = mod_utils.rainbow[color_index%len(mod_utils.rainbow)], label=event, lw = 0))
        color_index += 1
        bottoms = bottoms + numpy.array(event_amounts)
    plot.set_ylabel("mutation counts")
    plot.set_xlabel("read position")
    plot.set_xticks(numpy.array(mut_positions)[::5]+bar_width/2.0)
    plot.set_xticklabels(mut_positions[::5])

    plot = fig.add_subplot(224)#a stacked bar graph of positional mutation rates normalized to coverage
    bottoms = [0]*len(mut_positions)
    bottoms = numpy.array(bottoms)
    plot_layers = []
    color_index = 0
    for event in top_events:
        event_amounts = [mutations_by_position[position][event]/positional_coverage[position] if event in mutations_by_position[position] else 0 for position in mut_positions]
        plot_layers.append(plot.bar(mut_positions, event_amounts, bar_width, bottom = bottoms, color = mod_utils.rainbow[color_index%len(mod_utils.rainbow)], label=event, lw = 0))
        color_index += 1
        bottoms = bottoms + numpy.array(event_amounts)
    plot.set_ylabel("mutation counts/coverage")
    plot.set_xlabel("read position")
    plot.set_xticks(numpy.array(mut_positions)[::5]+bar_width/2.0)
    plot.set_xticklabels(mut_positions[::5])
    lg=plt.legend(loc=2,prop={'size':10}, labelspacing=0.2)
    lg.draw_frame(False)





    plt.savefig(out_prefix + '.pdf', transparent='True', format='pdf')
    plt.clf()

def normalized_mutation_rates(mutation_counts, coverage_counts):
    normalized_mutations = {}
    for strand in mutation_counts:
        if not strand in normalized_mutations:
            normalized_mutations[strand] = {}
        for chromosome in normalized_mutations[strand]:
            if not chromosome in normalized_mutations[strand]:
                normalized_mutations[strand][chromosome] = {}
            for position in normalized_mutations[strand][chromosome][position]:
                normalized_mutations[strand][chromosome][position] = \
                    float(mutation_counts[strand][chromosome][position])/float(coverage_counts[strand][chromosome][position])
    return normalized_mutations


def count_reads(lib_settings):
    """
    """

    # Create empty dicts for storing counts data
    srt_dict = createStrandDict(strands) # Counts for 5' end of read our standard data format
    cov_dict = createStrandDict(strands) # Counts of times covered by a read
    mut_dict = createStrandDict(strands) # Counts of mismatches at a position
    read_mutations = defaultdict(int) #counts different types of mutations relative to read
    genome_mutations = defaultdict(int) #counts different types of mutations relative to genome
    mutations_by_read_position = defaultdict(dict)
    read_position_coverage = defaultdict(float)
    mutations_by_genome_position = defaultdict(dict)
    genome_position_coverage = defaultdict(float)
    mutated_nts = defaultdict(float)
    read_insertion_sizes = []
    genomic_deletion_sizes = []
    with gzip.open(lib_settings.get_mapped_reads_sam_gz(), 'r') as f:
        for line in f: # Iterate through SAM file lines
            if not line.startswith('@'):
                # Parse line into relevant strings

                fields = line.strip().split('\t')
                ID = fields[0] #the first field in the mapped file corresponds to a unique id number for that read- these should correspond to the names in the raw_seqs dictionary
                flag = int(fields[1])
                '''
                The flag field provides a lot of info about the read, it is the decimal representation of a bit string, each digit of which is true or false

                Bit 0 = The read was part of a pair during sequencing
                Bit 1 = The read is mapped in a pair
                Bit 2 = The query sequence is unmapped
                Bit 3 = The mate is unmapped
                Bit 4 = Strand of query (0=forward 1=reverse)
                So, to see if a flag represents a read on the  - strand, we evaluate (16 & 'flag'), where & is the bitwise and operator,
                which will be non-zero (True) only if this read is on the - strand
                '''
                if (4&flag):#if this is an unmapped read, don't bother
                    continue

                if (16&flag):
                    strand = '-'
                else:
                    strand = '+'
                chrom = fields[2]
                MAPQ = int(fields[4])
                if int(MAPQ) >= lib_settings.get_property('min_mapping_quality'):
                    cigarString = fields[5]
                    seq = fields[9]
                    mappingLength = len(seq)
                    qScores = fields[10]
                    # Some lines seem to lack some strings this throws of indexing of NM:i, MD:Z, and NH:i strings
                    NHstr = checkTag('NH:i:',fields)
                    NMstr = checkTag('NM:i:',fields)
                    MDstr = checkTag('MD:Z:',fields)
                    assert 'NM:i' in NMstr
                    assert 'MD:Z' in MDstr
                    assert 'NH:i' in NHstr
                    multiplicity = float(NHstr.split(':')[2])

                    fields = line.strip().split('\t')
                    counts = float(1.0/multiplicity) # Weight of read
                    MDzString = MDstr

                    # Add subdicts for chromosome if needed
                    if chrom not in srt_dict[strand]:
                        srt_dict[strand][chrom] = defaultdict(float)
                    if chrom not in cov_dict[strand]:
                        cov_dict[strand][chrom] = defaultdict(float)
                    if chrom not in mut_dict[strand]:
                        mut_dict[strand][chrom] = defaultdict(float)

                    # Parse cigar string, get genome mapping span, and relative genomic positions covered by read
                    rel_genomic_event_positions, rel_genome_coverage, mutations_rel_genome, mutations_rel_read, readMappingSpan, genomeMappingSpan = parse_MDz_and_cigar(cigarString, MDzString, mappingLength, seq)

                    for pos in range(len(mutations_rel_genome)):
                        genome_position_coverage[pos] += counts
                        event = mutations_rel_genome[pos]
                        if not event == 'M': #count if it's not a match
                            assert event[0] != 'I'
                            if event[0] == 'D':
                                genomic_deletion_sizes.append(event[1])
                                event = event[0]
                            if event not in mutations_by_genome_position[pos]:
                                 mutations_by_genome_position[pos][event] = 0
                            mutations_by_genome_position[pos][event] += counts
                            genome_mutations[event] += counts
                            if event[0] in 'ATCG':
                                mutated_nts[event[0]] += counts

                    for pos in range(len(mutations_rel_read)):
                        read_position_coverage[pos] += counts
                        event = mutations_rel_read[pos]
                        if not event == 'M': #count if it's not a match
                            assert event[0] != 'D'
                            if event[0] == 'I':
                                read_insertion_sizes.append(event[1])
                                event = event[0]
                            if event not in mutations_by_read_position[pos]:
                                 mutations_by_read_position[pos][event] = 0
                            mutations_by_read_position[pos][event] += counts
                            read_mutations[event] += counts

                    # Set start position of read
                    if strand== '+':
                        start=int(fields[3])
                    else:
                        #When a read maps to the minus strand, bowtie returns the reverse complement, and indicates
                        # where this reverse mapped on the + strand. Thus the original 5' end of the read actually
                        # was x nt downstream on the + strand
                        start=int(fields[3])+genomeMappingSpan-1

                    # translate relative positions to absolute positions
                    genome_cov = readGenomicCoverage(rel_genome_coverage, strand, start) # get genome coverage

                    srt_dict[strand][chrom][start] += counts #just add the number of counts to that start position
                    for pos in [p for p in genome_cov]: # Increment positions for coverage dict
                        cov_dict[strand][chrom][pos] += counts

                    # If mismatches need to parse, get the absolute genomic pos, and increment counters
                    genMismatches = readGenomicCoverage(rel_genomic_event_positions, strand, start)
                    for event_position in genMismatches:
                        mut_dict[strand][chrom][event_position] += counts

    normalized_mutations = normalized_mutation_rates(mut_dict, cov_dict)
    mod_utils.makePickle(normalized_mutations, lib_settings.get_normalized_mutation_counts())


    mod_utils.makePickle(srt_dict, lib_settings.get_read_5p_counts())
    mod_utils.makePickle(cov_dict, lib_settings.get_positional_coverage())
    mod_utils.makePickle(mut_dict, lib_settings.get_mutation_counts())

    mod_utils.makePickle(genome_mutations, lib_settings.get_counting_prefix() + '.genome_mutations.pkl')
    mod_utils.makePickle(mutations_by_genome_position, lib_settings.get_counting_prefix() + '.genome_position_mutations.pkl')
    mod_utils.makePickle(genome_position_coverage, lib_settings.get_counting_prefix() + '.genome_position_coverage.pkl')

    mod_utils.makePickle(mutated_nts, lib_settings.get_counting_prefix() + '.nt_mutations.pkl')

    mod_utils.makePickle(read_mutations, lib_settings.get_counting_prefix() + '.read_mutations.pkl')
    mod_utils.makePickle(mutations_by_read_position, lib_settings.get_counting_prefix() + '.read_position_mutations.pkl')
    mod_utils.makePickle(read_position_coverage, lib_settings.get_counting_prefix() + '.read_position_coverage.pkl')

    mod_utils.makePickle(genomic_deletion_sizes, lib_settings.get_counting_prefix() + '.deletion_sizes.pkl')

    mod_utils.makePickle(read_insertion_sizes, lib_settings.get_counting_prefix() + '.insertion_sizes.pkl')

    plot_mutated_nts_pie(mod_utils.unPickle(lib_settings.get_counting_prefix() + '.nt_mutations.pkl'), 'mutated rRNA nts in ' + lib_settings.sample_name, lib_settings.get_counting_prefix()+'.mutated_nts' )
    plot_full_mutation_stats(mod_utils.unPickle(lib_settings.get_counting_prefix() + '.read_mutations.pkl'), mod_utils.unPickle(lib_settings.get_counting_prefix() + '.insertion_sizes.pkl'),
                             mod_utils.unPickle(lib_settings.get_counting_prefix() + '.read_position_mutations.pkl'), mod_utils.unPickle(lib_settings.get_counting_prefix() + '.read_position_coverage.pkl'), 'mutations wrt reads', "insertion size",
                             lib_settings.get_counting_prefix()+'.read_mutations')
    plot_full_mutation_stats(mod_utils.unPickle(lib_settings.get_counting_prefix() + '.genome_mutations.pkl'), mod_utils.unPickle(lib_settings.get_counting_prefix() + '.deletion_sizes.pkl'), mod_utils.unPickle(lib_settings.get_counting_prefix() + '.genome_position_mutations.pkl'),
                             mod_utils.unPickle(lib_settings.get_counting_prefix() + '.genome_position_coverage.pkl'), 'mutations wrt genome', "deletion size",
                             lib_settings.get_counting_prefix()+'.genome_mutations')
    pie_read_5p_ends(srt_dict, mod_utils.convertFastaToDict(lib_settings.experiment_settings.get_rRNA_fasta()), lib_settings.get_counting_prefix())

def test():
    """

    with uncommented print statements should get:
    parse_MDz_and_cigar('36M', 'MD:Z:1A0C0C0C1T0C0T27', 36, 'CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG')
    ['M', ('A', 'G'), ('C', 'A'), ('C', 'T'), ('C', 'A'), 'M', ('T', 'G'), ('C', 'G'), ('T', 'G'), 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M']
    ['M', ('A', 'G'), ('C', 'A'), ('C', 'T'), ('C', 'A'), 'M', ('T', 'G'), ('C', 'G'), ('T', 'G'), 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M']
    36 36
    parse_MDz_and_cigar('6M1I29M', 'MD:Z:0C1C0C1C0T0C27', 36, 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT')
    [('C', 'G'), 'M', ('C', 'G'), ('C', 'A'), 'M', ('C', 'G'), ('I', 1), ('T', 'G'), ('C', 'G'), 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M']
    [('C', 'G'), 'M', ('C', 'G'), ('C', 'A'), 'M', ('C', 'G'), ('T', 'G'), ('C', 'G'), 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M']
    36 35
    parse_MDz_and_cigar('9M9D27M', 'MD:Z:2G0A5^ATGATGTCA27', 36, 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC')
    ['M', 'M', ('G', 'T'), ('A', 'G'), 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M']
    ['M', 'M', ('G', 'T'), ('A', 'G'), 'M', 'M', 'M', 'M', 'M', ('D', 9), 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M']
    36 45
    """


    print "parse_MDz_and_cigar('36M', 'MD:Z:1A0C0C0C1T0C0T27', 36, 'CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG')"
    parse_MDz_and_cigar('36M', 'MD:Z:1A0C0C0C1T0C0T27', 36, 'CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG')
    print "parse_MDz_and_cigar('6M1I29M', 'MD:Z:0C1C0C1C0T0C27', 36, 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT')"
    parse_MDz_and_cigar('6M1I29M', 'MD:Z:0C1C0C1C0T0C27', 36, 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT')
    print "parse_MDz_and_cigar('9M9D27M', 'MD:Z:2G0A5^ATGATGTCA27', 36, 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC')"
    parse_MDz_and_cigar('9M9D27M', 'MD:Z:2G0A5^ATGATGTCA27', 36, 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC')

#test()