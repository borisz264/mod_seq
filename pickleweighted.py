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

import sys, cPickle as pickle, gzip, re
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
    out_fn = ''.join([fn[0]+suffix+'.p'])
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
                mutations_rel_genome[genome_counter] = (tag, seq[genome_counter])
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
    return genomic_event_positions, genomeCoverage, readMappingSpan, genomeMappingSpan

    
def main(args):
    """
    """
    sortedfile, fn = args
    
    # Create empty dicts for storing counts data
    srt_dict = createStrandDict(strands) # Counts for 5' end of read our standard data format
    cov_dict = createStrandDict(strands) # Counts of times covered by a read
    mut_dict = createStrandDict(strands) # Counts of mismatches at a position
    read_cov_dict = defaultdict(float) #count coverage at each position in read( basically, how many ready were at least this long
    read_mut_dict = defaultdict(float) #count mismatches or indels assigned to each read position
    types_of_mutations = defaultdict(float) #counts different types of mutations

    with gzip.open(sortedfile, 'r') as f:
        for line in f: # Iterate through SAM file lines
            # Parse line into relevant strings
            fields = line.strip().split('\t')
            counts = float(fields[0].split('&')[1]) # Weight of read
            strand = fields[1]
            chrom = fields[2].strip()
            seq = fields[4].strip()
            cigarString = fields[6].strip()
            NMi = int(fields[8].split(':')[2])
            MDzString = fields[9].strip()
            mappingLength = int(fields[7].strip())

            # Add subdicts for chromosome if needed
            if chrom not in srt_dict[strand]:
                srt_dict[strand][chrom] = defaultdict(float)
            if chrom not in cov_dict[strand]:
                cov_dict[strand][chrom] = defaultdict(float)
            if chrom not in mut_dict[strand]:
                mut_dict[strand][chrom] = defaultdict(float)
            
            # Parse cigar string, get genome mapping span, and relative genomic positions covered by read
            rel_genomic_event_positions, rel_genome_coverage, readMappingSpan, genomeMappingSpan = parse_MDz_and_cigar(cigarString, MDzString, mappingLength, seq)
            
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
            if NMi > 0:
                genMismatches = readGenomicCoverage(rel_genomic_event_positions, strand, start)
                for event_position in genMismatches:
                    mut_dict[strand][chrom][event_position] += counts
    
    pickleDict(srt_dict, fn, '_5p')
    pickleDict(cov_dict, fn, '_cov')
    pickleDict(mut_dict, fn, '_mut')
    
#main (sys.argv[1:])

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

test()