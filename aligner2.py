"""
aligner2.py
Author:     Wei Lu
Year:        2016
"""
from __future__ import division
import sys
from Bio import SeqIO

#Support class to execute the handiwork of read alignment to a reference sequence
class Aligner2:
    """Represents a process of aligning read to a reference sequence.
    Aligner1 represents exact match of read to reference sequence"""
    
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

#Support Class to represent alignments and check relationships between alignments
class Alignment:
    """Represents an alignment result information"""
    
    def __init__(self, chr, pos, strand):
        """Constructor creates an alignment with reference chromosome name, left-most position
        of exact read match and the matched strand."""
        self._chr = chr
        self._pos = pos
        self._strand = strand
    
    def getPos(self):
        return self._pos
    
    def getStrand(self):
        return self._strand
    
    def is_concordant_with(self, other):
        """Method to check if paired-end reads are mapped correctly, which means if one read
        maps on forward strand and the other read should map on reverve strand of reference
        genome"""
        if (self._strand == "+" and other._strand == "-") or (self._strand == "-" and other._strand == "+"):
            return True
        else: return False

    def __repr__(self):
        """Return the string representation of the student."""
        return self._chr+"\t"+str(self._pos)+"\t"+self._strand

#Support class to hold information of concordant alignment pairs
class Concordant:
    '''This class repersents the information of any aligned read pairs'''
    
    def __init__(self, alignment1, total_alignment1, alignment2, total_alignment2):
        '''Constructor creates an '''
        self._alignment1 = alignment1
        self._alignment2 = alignment2
        self._total_alignment1 = total_alignment1
        self._total_alignment2 = total_alignment1
    
    def __repr__(self):
        """Return the string representation of the concordant read pairs"""
        return str(self._alignment1)+"\t"+str(self._total_alignment1)+"\t"+str(self._alignment2)+"\t"+str(self._total_alignment2)

#Check the command line arguments
if len(sys.argv) < 5:
    print "Usage: <reference file (fasta)> <read file 1 (fasta)> <read file 2 (fasta)> <out file>"
    sys.exit(0)

#Read the reference sequence and initiate the aligner
try:
    for s in SeqIO.parse(sys.argv[1], "fasta"):
        aligner = Aligner2(s)
        break #Stop after the fist sequence in the reference
except IOError as e:
    print "Could not read reference sequence file (see below)!"
    print e
    sys.exit(1)

reads1 = [read for read in SeqIO.parse(sys.argv[2], "fasta")]
reads2 = [read for read in SeqIO.parse(sys.argv[3], "fasta")]
pairs = zip(reads1, reads2)

output = open(sys.argv[4], 'w')
output.write("READ_NAME\tREF_NAME1\tPOS1\tSTRAND1\tNUMBER_OF_ALIGNMENTS1\tREF_NAME2\tPOS2\tSTRAND2\tNUMBER_OF_ALIGNMENTS2\n")
concordant_match = 0

for pair in pairs:
    position1_f = aligner.align_forward(pair[0])
    position1_r = aligner.align_reverse(pair[0])
    position2_f = aligner.align_forward(pair[1])
    position2_r = aligner.align_reverse(pair[1])
    
    #print position1_f, position1_r
    Alignment_read1 = []
    if len(position1_f)==0 and len(position1_r)==0:
        Alignment_read1=Alignment_read1
    else:
        for pos in position1_f:
            a = Alignment(aligner._refname, pos, "+")
            Alignment_read1.append(a)
        for pos in position1_r:
            a = Alignment(aligner._refname, pos, "-")
            Alignment_read1.append(a)
                  
    Alignment_read2 = []
    if len(position2_f)==0 and len(position2_r)==0:
        Alignment_read2=Alignment_read2
    else:
        for pos in position2_f:
            a = Alignment(aligner._refname, pos, "+")
            Alignment_read2.append(a)
        for pos in position2_r:
            a = Alignment(aligner._refname, pos, "-")
            Alignment_read2.append(a)

    concordants = []
    
    for a in Alignment_read1:
        for b in Alignment_read2:
            if a.is_concordant_with(b):
                concordant = Concordant(a, len(Alignment_read1), b, len(Alignment_read2))
                concordants.append(concordant)

    if len(concordants)!=0:
        output.write("%s\t%s\n" % (pair[0].name, concordants[0]))
        concordant_match +=1
    else:
        output.write("%s\t%s\t*\t0\t%s\t%s\t*\t0\t%s\n" % (pair[0].name, aligner._refname, len(Alignment_read1), aligner._refname,len(Alignment_read2))) 

output.close() 

print "%d read pairs have been aligned to reference sequence" % len(pairs)
print "%d read pairs matched as concordants to reference sequence %f %%" % (concordant_match, (concordant_match/len(pairs))*100)
print "%d read pairs failed to match concordantly to reference sequence %f %%" % (len(pairs)-concordant_match, (len(pairs)-concordant_match)/len(pairs)*100)
 
print "Done"