
"""
aligner1.py
Author:    Wei Lu
Year:        2016
"""
from __future__ import division
import sys
from Bio import SeqIO

#Support Class to represent alignments and check relationships between alignments
class Alignment(object):
    """Represents an alignment result information"""
    
    def __init__(self, chr, pos, strand):
        """Constructor creates an alignment with reference chromosome name, left-most position
        of exact read match and the matched strand."""
        self._chr = chr
        self._pos = pos
        self._strand = strand
    
    def getPos(self):
        return self._pos
    
    def is_concordant_with(self, other):
        """Method to check if paired-end reads are mapped correctly, which means if one read
        maps on forward strand and the other read should map on reverve strand of reference
        genome"""
    
        return False

    def __repr__(self):
        """Return the string representation of the student."""
        return self._chr+"\t"+str(self._pos)+"\t"+self._strand

#Support class to execute the handiwork of read alignment to a reference sequence
class Aligner1(object):
    '''
    Represents a process of aligning read to a reference sequence.
    Aligner1 represents exact match of read to reference sequence
    '''
    
    def __init__(self, ref):
        """Constructor create an aligner that contain reference sequence and id. """
        self._refname = ref.id
        self._refseq = ref.seq
        self._rrefseq = ref.seq.reverse_complement()
        
    def align_forward(self, read):
        """Method to find exact matches on forward reference sequence. Return a list of integer
        that represent position of all matches"""
        alignments = []
        location = -1
        while True:
            location = self._refseq.find(read.seq, location+1)
            if location == -1: break
            alignments.append(location)
        return alignments
    
    def align_reverse(self, read):
        """Method to find exact matches on reverse strand of reference sequence"""
        alignments = []
        location = -1
        while True:
            location = self._refseq.find(read.seq.reverse_complement(), location+1)
            if location == -1: break
            alignments.append(location)
        return alignments


#Check the command line arguments
if len(sys.argv) < 4:
    print "Usage: <reference file (fasta)> <read file (fasta)> <out file>"
    sys.exit(0)

#Read the reference sequence and initiate the aligner
try:
    for s in SeqIO.parse(sys.argv[1], "fasta"):
        aligner = Aligner1(s)
        break #Stop after the fist sequence in the reference
except IOError as e:
    print "Could not read reference sequence file (see below)!"
    print e
    sys.exit(1)

output = open(sys.argv[3], 'w')
output.write("READ_NAME\tREF_NAME\tPOS\tSTRAND\tNUMBER_OF_ALIGNMENTS\n")

#Support function to sort alignment object by their postion
def alignmentSort(alignments):
    for fillslot in range(len(alignments)-1,0,-1):
        positionOfMax=0
        for location in range(1,fillslot+1):
            if alignments[location].getPos()>alignments[positionOfMax].getPos():
                positionOfMax = location

        temp = alignments[fillslot]
        alignments[fillslot] = alignments[positionOfMax]
        alignments[positionOfMax] = temp

totalReads = 0  
match1 = 0
match2 = 0
match3 = 0
for read in SeqIO.parse(sys.argv[2], "fasta"):
    alignments = []
    positionsf = aligner.align_forward(read)
    for pos in positionsf:
        alignments.append(Alignment(aligner._refname,  pos, "+"))
    positionsr = aligner.align_reverse(read)
    for pos in positionsr:
        alignments.append(Alignment(aligner._refname, pos, "-"))
    totalReads += 1
    if len(alignments) == 1:
        match1 += 1
        output.write("%s\t%s\t%d\n" % 
                     (read.name, str(alignments[0]), len(alignments)))
    elif len(alignments) > 1:
        match2 += 1
        alignmentSort(alignments)
        output.write("%s\t%s\t%d\n" % 
                     (read.name, str(alignments[0]), len(alignments)))
    else:
        match3 += 1
        output.write("%s\t*\t%s\t*\t%s\n" % (read.name , str(0) ,  str(0)))
    
output.close()

print "%d reads have been aligned to reference sequence" % totalReads
print "%d reads have failed to align to reference sequence %f %%" % (match3, (match3/totalReads)*100)
print "%d reads has matched exactly once to reference sequence %f %%" % (match1, (match1/totalReads)*100)
print "%d reads has matched more than once to reference sequence %f %%" % (match2, (match2/totalReads)*100)


print "Done"