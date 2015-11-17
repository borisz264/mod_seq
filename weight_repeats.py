#to convert from the mapped reads in the bowtie output to the old Bowtie1 format, and also computes mapping multiplicity
#The SAM file is read from stdin (piped in from SAMtools)
#2015-05-01 TMC: Now use NH:i string, not MAPQ score to calculate multiplicity of mapping.
#2015-05-01 TMC: Now writes NM:i, and MD:z strings to output, this will allow mapping of mismatches/mutations in count_reads_and_mismatches.py

import sys, gzip

outFile = sys.argv[1]

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
    return [field for field in fields[7:] if tag in field][0]

h=gzip.open(outFile, 'w')
print "weighting file opened"
for line in sys.stdin:#this is the SAM file
    print "sam file opened"
    if not line.startswith('@'):
        fields = line.strip().split('\t')
        oldID = fields[0] #the first field in the mapped file corresponds to a unique id number for that read- these should correspond to the names in the raw_seqs dictionary
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
        if (16&flag):
            strand = '-'
        else:
            strand = '+'
        chrom = fields[2]
        position = fields[3]
        MAPQ = int(fields[4])
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
        #multip = float(mapQ_to_multiplicity[MAPQ])
        multip = float(NHstr.split(':')[2])
        weighted_id = '%s&%f' % (oldID, 1.0/multip)
        newLine = '\t'.join([weighted_id, strand, chrom, position, seq, qScores, cigarString, str(mappingLength), NMstr, MDstr])
        h.write(newLine+'\n')

h.close()