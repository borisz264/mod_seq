__author__ = 'boris'
"""
takes the info (key) file from the hopkins sequencing core, and prints formats convenient to this pipeline
"""

import sys

infile = sys.argv[1]
f = open(infile, 'rU')
lib_names = []
file_names = []
for line in f.readlines()[1:]:
    ll = line.strip('\n').split(',')
    Project, FCID, Lane, Index, Library_Name, File_Name = ll
    lib_names.append('"'+Library_Name+'"')
    file_names.append('"'+File_Name+'_1.fastq.gz'+'"')
f.close()

print 'fastq_gz_files = [%s]' % (','.join(file_names))
print 'sample_names = [%s]' % (','.join(lib_names))
