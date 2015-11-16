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

def parseCigarString(cigarString, mappingLength):
    '''
    arguments:
        cigarString - a string of [#][A-Z][#][A-Z] etc... that described the alignment of a read to the genome
        mappingLength - an int of the r
    returns:
        the genomic distance spanned by this read
        
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
    '''
    
    #Parse cigar strings into numbers & tags
    numbers = re.split('[A-Z,=]', cigarString)[:-1] #not sure why, bu this seems to produce an extra blank entry thay I'm stripping off
    tags = re.split('[0-9]*', cigarString)[1:] 
    numbers = [int(n) for n in numbers]
    
    assert len(numbers) == len(tags)
    
    genomeMulitpliers = {'M':1, 'I':0, 'D':1, 'N':1, '=':1, 'X':1}
    readMulitpliers = {'M':1, 'I':1, 'D':0, 'N':0, '=':1, 'X':1}
    
    # Initialize counters
    readMappingSpan, genomeMappingSpan, genome_pos = 0, 0, 0
    genomeCoverage = []
    
    #print 'New Read'
    cigar_tups = zip(tags, numbers)
    for cigar in cigar_tups:
        tag, tag_len = cigar
        genome_multiplier = genomeMulitpliers[tag]
        read_multiplier = readMulitpliers[tag]
        readMappingSpan += read_multiplier*tag_len
        genomeMappingSpan += genome_multiplier*tag_len
        
        # Determines which  genomic positions relative to the read 5' end are covered by a read
        # For M: positions covered by read: so add to list, and increment counter
        # For N: positions spanned by read: increment counter 
        # Insertions and deletions are ignored,  a read with an indel never sees this function
        if genome_multiplier: # Increment counter if genome multiplier = 1, ie position spanned by read
            if tag == 'M' or tag =='=' or tag == 'X':
                genomeCoverage += range(genome_pos,genome_pos+tag_len)
            genome_pos += tag_len
    
    assert readMappingSpan == mappingLength
    assert mappingLength == len(genomeCoverage)

    return readMappingSpan, genomeMappingSpan, genomeCoverage

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
        pos_l = [p+read_start if p != 'I' else 'I' for p in rel_cov]
    elif strand == '-':
        pos_l = [read_start-p if p != 'I' else 'I' for p in rel_cov]
    return pos_l

def parseMDz(MDzString, NMi, genome_cov, seq):
    """
    Assumes:
        MDzString, the MD:z string for the current read
        genome_cov, a list relating position in read (i) to absolute genomic position pos_l[i]
        seq, a string representation of the read
    Returns:
        A list of positions with mismatches in absolute genomic position
    """
    
    #The MD field aims to achieve SNP/indel calling without looking at the reference.
    #For example, a string `10A5^AC6' means from the leftmost reference base in the
    #alignment, there are 10 matches followed by an A on the reference which is 
    #different from the aligned read base; the next 5 reference bases are matches
    #followed by a 2bp deletion from the reference; the deleted sequence is AC;
    #the last 6 bases are matches. The MD field ought to match the CIGAR string.
    #The string '0T0C37T' indicates the first base of the read is a mismatch from
    # the reference sequence of T, the second base is a mismatch from the reference
    # sequence C, followed by 37 matches, with a final mismatch from the reference A
    
    MDz = MDzString.split(':')[2]
    MD_tags = re.findall('[0-9]+|[A-Z,^]+',MDz)
    MD_tags = [int(x) if re.match('[0-9]+',x) else x for x in MD_tags]
    
    mm_pos = [] # A list to store positions of mismatches
    read_counter = 0 # A counter to store the current position of the read (0 indexed)
    for tag in MD_tags:
        if tag in ['A','T','C','G']: #If the tag is a base it indicates a mismatch
            if seq[read_counter] != 'N': # We want to ignore 'N' bases in sequence
                mm_pos.append(read_counter) # Append the current position in the read
            read_counter += 1 #Increment the counter
        elif isinstance(tag, int): # If the tag is an int, this represents matches
            read_counter += tag # Increment counter by number of matches
    
    return mm_pos
    
def main(args):
    """
    """
    sortedfile, fn = args
    
    # Create empty dicts for storing counts data
    srt_dict = createStrandDict(strands) # Counts for 5' end of read our standard data format
    cov_dict = createStrandDict(strands) # Counts of times covered by a read
    mut_dict = createStrandDict(strands) # Counts of mismatches at a position
    types_of_mutations = defaultdict(dict) #counts different types of mutations

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
            
            # Skip reads with insertions or deletions. Note, some code still accomodates
            # insertions, represented as 'I', but not mismatch parsing. Deletions haven't
            # been addressed
            #Boris- I want indels to be accommodated

            #if 'I' in cigarString or 'D' in cigarString:
            #    continue
            #
            '''
            TODO:

            '''
            # Add subdicts for chromosome if needed
            if chrom not in srt_dict[strand]:
                srt_dict[strand][chrom] = defaultdict(float)
            if chrom not in cov_dict[strand]:
                cov_dict[strand][chrom] = defaultdict(float)
            if chrom not in mut_dict[strand]:
                mut_dict[strand][chrom] = defaultdict(float)
            
            # Parse cigar string, get genome mapping span, and relative genomic positions covered by read
            readMappingSpan, genomeMappingSpan, relGenomeCoverage = parseCigarString(cigarString, mappingLength)
            
            # Set start position of read
            if strand== '+':
                start=int(fields[3])
            else:
                #When a read maps to the minus strand, bowtie returns the reverse complement, and indicates
                # where this reverse mapped on the + strand. Thus the original 5' end of the read actually
                # was x nt downstream on the + strand
                start=int(fields[3])+genomeMappingSpan-1
            
            # translate relative positions to absolute positions
            genome_cov = readGenomicCoverage(relGenomeCoverage, strand, start) # get genome coverage
               
            srt_dict[strand][chrom][start] += counts #just add the number of counts to that start position
            for pos in [p for p in genome_cov if p != 'I']: # Increment positions for coverage dict
                cov_dict[strand][chrom][pos] += counts
            
            # If mismatches need to parse, get the absolute genomic pos, and increment counters
            if NMi > 0:
                relMismatches = parseMDz(MDzString, NMi, genome_cov, seq)
                genMismatches = readGenomicCoverage(relMismatches, strand, start)
                for mm in genMismatches:
                    mut_dict[strand][chrom][pos] += counts
    
    pickleDict(srt_dict, fn, '_5p')
    pickleDict(cov_dict, fn, '_cov')
    pickleDict(mut_dict, fn, '_mut')
    
main (sys.argv[1:])




