import os
import operator
import itertools
import gzip
import numpy as np
from scipy import stats
import cPickle as pickle
import math
import multiprocessing
'''
Colorblind safe colors from Bang Wong, Nature Methods 8. 441 (2011)
'''
black = (0,0,0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)
colors = [black, orange, skyBlue, bluishGreen, vermillion, blue, reddishPurple, yellow]
rainbow = [black, vermillion, orange, bluishGreen, blue, reddishPurple, 'violet']
markers = ['.', 'o', 'v', 's', '^', 'p', 'x', '+']
line_styles = ['solid', 'dashed', 'dotted']

bokeh_black = (0,0,0)
bokeh_orange = (230,159,0)
bokeh_skyBlue = (86,180,233)
bokeh_bluishGreen = (0,158,115)
bokeh_yellow = (240,228,66)
bokeh_blue = (0,114,178)
bokeh_vermillion = (213,94,0)
bokeh_reddishPurple = (204,121,167)

#parralellization code from
# http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class
def spawn(f):
    def fun(q_in,q_out):
        while True:
            i,x = q_in.get()
            if i == None:
                break
            q_out.put((i,f(x)))
    return fun

def parmap(f, X, nprocs = multiprocessing.cpu_count()):
    q_in   = multiprocessing.Queue(1)
    q_out  = multiprocessing.Queue()
    proc = [multiprocessing.Process(target=spawn(f),args=(q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]
    [p.join() for p in proc]
    return [x for i,x in sorted(res)]

def unPickle(fileName):
    #returns the pickled object stored in a pickle file
    f = open(fileName, 'r')
    o = pickle.load(f)
    f.close()
    return o

def makePickle(o, fileName):
    f = open(fileName, 'w')
    pickle.dump(o, f)
    f.close()


def make_dir(dirname):
    """
    Makes the directory; doesn't throw an error if it exists.
    """
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except:
            print 'The directory was made by another thread extremely recently.'

##########
#MATH
##########
def divideWithError(num, stdDevNum, denom, stdDevDenom):
    '''
    divides the two values with provided standard deviations, and returns the mean and error of the ratio using standard error propogation
    '''
    num = float(num)
    denom = float(denom)
    stdDevNum = float(stdDevNum)
    stdDevDenom = float(stdDevDenom)

    ratio = num/denom
    ratioError = ratio*math.sqrt((stdDevNum/num)**2+(stdDevDenom/denom)**2)

    return ratio, ratioError

def subtractWithError(num, stdDevNum, denom, stdDevDenom):
    '''
    divides the two values with provided standard deviations, and returns the mean and error of the ratio using standard error propogation
    '''
    num = float(num)
    denom = float(denom)
    stdDevNum = float(stdDevNum)
    stdDevDenom = float(stdDevDenom)

    ratio = num/denom
    ratioError = ratio*math.sqrt((stdDevNum/num)**2+(stdDevDenom/denom)**2)

    return ratio, ratioError

def next_square_number(number):
    return int(math.ceil(math.sqrt(number)))**2

def computePfromMeanAndStDevZscore(mean, standard_deviation, testValue):
    #computes probability that test value came from a gaussian with the given mean and standard deviation
    try:
        z = (float(testValue)-mean)/standard_deviation
        p = stats.norm.sf(z)
        return p, z
    except ZeroDivisionError:
        return 0.5, 0

def filter_x_y_pairs(x, y, filter_list = [float('inf'), -1*float('inf')]):
    '''
    takes 2 paired arrays, and returns matched copies of them with any positions with values in
    filter_list removed from both arrays, to keep them synced.
    alos removes NaN (defined by testing if the entry equals itself, which fails for NaN)
    :param filter_list: list of values to remove
    :return:
    '''
    filtered_x, filtered_y = [], []
    assert len(x) == len(y)
    for i in range(len(x)):
        if x[i] not in filter_list and y[i] not in filter_list and x[i]==x[i] and y[i]==y[i]:
            filtered_x.append(x[i])
            filtered_y.append(y[i])
    return np.array(filtered_x), np.array(filtered_y)

def ranges_overlap(min1, max1, min2, max2):
    """

    :param min1:
    :param max1:
    :param min2:
    :param max2:
    :return: return True if the 2 ranges overlap (edge inclusive), else False
    """

    if min1 <= max2 and min2 <= max1:
        return True
    return False

def file_exists(fname):
    """
    makes sure a given file exists
    """
    if not os.path.exists(fname):
        return False
    fstats = os.stat(fname)
    if not fstats[6]:
        return False
    if not os.access(fname, os.R_OK):
        raise ValueError('Input File %s cannot be read' % fname)
    return True


def getBinIndex_soft_upper(v, bins):
    for i in range(len(bins) - 1):
        if v > bins[i] and v <= bins[i + 1]:
            return i
    return -1


def get_barcode(line):
    """
    - Extracts the barcode from the first line of a fastq quartet
        - Assumes the first line is of the form:
            @D5FF8JN1:4:1101:1220:2099#ACTTGA/1
    """
    return line.split('#')[-1].split('/')[0]


def get_index_from_kmer(kmer):
    """
    returns the base10 version of the base 4 DNA representation
    """
    index = 0
    base2face = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, base in enumerate(kmer):
        if not base in 'ACGT':
            return -1
        power = len(kmer) - 1 - i
        index += base2face[base] * (4 ** power)
    return index

def convertFastaToDict(fastaFile):
    '''
    converts a fasta file to a dict of {sequenceName:sequence}
    '''
    currentName = None
    currentSequence = None
    seqDict = {}
    f = open(fastaFile)
    for line in f:
        if not line.strip() == '' and not line.startswith('#'):#ignore empty lines and commented out lines
            if line.startswith('>'):#> marks the start of a new sequence
                if not currentName == None: #after we've reached the first > line, we know what the sequence corresponds to
                    seqDict[currentName] = rna_to_dna(currentSequence.upper())
                currentName = line.strip()[1:]
                currentSequence = ''
            else:
                currentSequence += line.strip()
    f.close()
    seqDict[currentName] = rna_to_dna(currentSequence.upper())

    return seqDict

def get_kmer_from_index(kmax, index):
    """
    takes a number (essentially base 4)
    and returns the kmer it corresponds to in alphabetical order
    eg.
    AAAA = 0*1
    CA = 4*4 + 0*1
    GC = 3*4 + 1 * 1
    """
    bases = 'ACGT'
    out = ''
    for k in range(kmax - 1, -1, -1):
        face, index = divmod(index, 4 ** k)
        out += bases[face]
    return out


def yield_kmers(k):
    """
    An iterater to all kmers of length k in alphabetical order
    """
    bases = 'ACGT'
    for kmer in itertools.product(bases, repeat=k):
        yield ''.join(kmer)


def aopen(file, mode='r'):
    if file[-3:] == '.gz':
        return gzip.open(file, mode + 'b')
    else:
        return open(file, mode)


def hamming_N(str1, str2):
    if not len(str1) == len(str2):
        raise(ValueError, 'lengths don\'t match')
    str1 = str1.upper()
    str2 = str2.upper()
    str1 = str1.replace('N', '#')
    return sum(itertools.imap(operator.ne, str1, str2))


# from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))


def close_float_value(a, b, max_percent=1.0):
    if a == 0 and b == 0:
        return True
    if not (a > 0 and b > 0):
        return False
    ratio = float(max(a, b)) / float(min(a, b))
    percent_increase = (ratio - 1.0) * 100.0
    return percent_increase < max_percent


def significantly_enriched(xs, zthresh=2., scale='linear'):
    assert scale in ['linear', 'log']
    if scale =='log':
        xs = np.log2(xs)
    xs = stats.zscore(xs)
    return [x > zthresh for x in xs]

def getAllMismatchedSeqs(kmer, mismatchPositions):
    nucs = ['A', 'C', 'G', 'T']
    #generate tuples of allowed nucs at each mismatch position using a recursive algorithm
    allowedNucs = {}
    mismatchPositions = np.array(mismatchPositions)
    assert len(set(mismatchPositions)) == len(mismatchPositions)
    if len(mismatchPositions) == 0:
        yield kmer
    else:
        mismatchNucs = [] + nucs
        #print kmer
        #print mismatchPositions
        #print mismatchPositions[0]
        #print kmer[mismatchPositions[0]]
        mismatchNucs.remove(kmer[mismatchPositions[0]])
        downstreamMismatchSeqs = getAllMismatchedSeqs(kmer[mismatchPositions[0]+1:], mismatchPositions[1:]-(mismatchPositions[0]+1))
        for mismatchNuc in mismatchNucs:
            for downstreamMismatchSeq in downstreamMismatchSeqs:
                returnSeq = kmer[:mismatchPositions[0]] + mismatchNuc +downstreamMismatchSeq
                assert len(returnSeq) == len(kmer)
                yield returnSeq

def getPaddedMismatchedAdjacentKmers(kmerSequence, padding, numMismatches):
    '''
    Yield all sequences of length (len(kmerSequence)+padding )that contain the given kmer, with exactly the given number of mismatches.
    The order yielded is as follows:
        First mismatches are allowed at position 0 to (numMismatches-1)
            For each register:
                Iterate through all possible nucs at mismatch position in alphabetical order
                    Iterate through each nucleotide in padding positions in alphabetical order.
                Shift to next register
            Move most 3' mismatch position down by one, but not past the end of the kmerSequence if end of KmerSequence
            is reached, shift secondmost 3' mismatch 1 nt 3', and reset most 3' mismatch to 1nt 3' of that one
    '''

    # for troubleshooting, want to check that no repeats are generated, so will assert that size of this list and set
    # must be the same
    kmer_set = set()
    kmer_list =[]
    upper_to_combined = {}
    nucs = 'ACGT'
    #initialize mismatchPositions
    #print numMismatches
    if numMismatches == 0:
        for mismatchedKmer in [kmerSequence]:
            for shift in range(padding+1):
                #generate all possible mismatches to the kmer
                for leftPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = shift)]:
                    for rightPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = padding-shift)]:
                        paddedSeq = leftPaddingSeq+mismatchedKmer+rightPaddingSeq
                        if paddedSeq not in kmer_set:
                            kmer_list.append(paddedSeq)
                        kmer_set.add(paddedSeq)
    else:
        mismatchPositionsList = itertools.combinations(range(len(kmerSequence)), numMismatches)
        for mismatchPositions in mismatchPositionsList:
            #print mismatchPositions
            for mismatchedKmer in getAllMismatchedSeqs(kmerSequence, mismatchPositions):
                for shift in range(padding+1):
                    #generate all possible mismatches to the kmer
                    for leftPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = shift)]:
                        for rightPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = padding-shift)]:
                            paddedSeq = leftPaddingSeq+mismatchedKmer+rightPaddingSeq
                            paddedUpper = paddedSeq.upper()
                            if paddedUpper not in kmer_set:
                                kmer_list.append(paddedUpper)
                            kmer_set.add(paddedUpper)

    #print kmer_list
    #print kmer_set
    #print len(kmer_list), len(kmer_set)
    #assert len(kmer_list) == len(kmer_set)

    return kmer_list


def revComp(seq, isRNA = False):
    seq = seq.upper()
    compDict = {'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C', 'N':'N', '-':'-', '.':'.', '*':'*'}
    revComp = ''.join([compDict[c] for c in seq[::-1]])
    if isRNA:
        return revComp.replace('T', 'U')
    return revComp

def rna_to_dna(dna_seq):
    return dna_seq.replace('U','T').replace('u', 't')

def dna_to_rna(dna_seq):
    return dna_seq.replace('T','U').replace('t', 'u')

def parse_wig(filename):
    mydict = {}
    for line in open(filename):
        if "variableStep" in line:
            chrom = line.strip().split('=')[-1]
            mydict[chrom] = {}

        else:
            ll = line.strip().split('\t')
            mydict[chrom][int(ll[0])] = float(ll[1])

    return mydict