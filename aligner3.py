"""
aligner3.py
Author:    Wei Lu
Year:        2016
"""
from __future__ import division
import sys
from Bio import SeqIO
import numpy as np

#Support Class to represent alignments and check relationships between alignments
class Alignment(object):
    """Represents an alignment result information"""
    
    def __init__(self, chr, pos, strand, cost):
        """Constructor creates an alignment with reference chromosome name, left-most position
        of exact read match and the matched strand."""
        self._chr = chr
        self._pos = pos
        self._strand = strand
        self._cost = cost
    
    def getPos(self):
        return self._pos
    def getStrand(self):
        return self._strand
    def getCost(self):
        return self._cost
    def getRefname(self):
        return self._chr
    
    def is_concordant_with(self, other):
        """Method to check if paired-end reads are mapped correctly, which means if one read
        maps on forward strand and the other read should map on reverse strand of reference
        genome"""
    
        return False

    def __repr__(self):
        """Return the string representation of the student."""
        return self._chr+"\t"+str(self._pos)+"\t"+self._strand+"\t"+str(self._cost)

#Penalty score dictionary for alignment
pt ={'match': 0, 'mismatch': 1, 'gap': 1}

#Supporting function to compare single letter of two string
def mch(alpha, beta):
    if alpha == beta:
        return pt['match']
    elif alpha == '-' or beta == '-':
        return pt['gap']
    else:
        return pt['mismatch']

#Function for filling alignment matrix 
def fillMatrix(s1, s2):
    m, n = len(s1), len(s2)
    score = np.zeros((m+1, n+1))
    
    #Initialization
    for i in range(m+1):
        score[i][0] = pt['gap'] * i
    for j in range(n+1):
        score[0][j] = 0 #Initialize first row as 0 

    #Fill
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = score[i-1][j-1] + mch(s1[i-1], s2[j-1])
            delete = score[i-1][j] + pt['gap']
            insert = score[i][j-1] + pt['gap']
            score[i][j] = min(diag, delete, insert)            
    return score

#Function for tracing back all possible alignments
def traceBack(score,s1,s2): 
    #Find the back tracing start points in alignment matrix
    a = np.asarray(score[len(s1)-1])
    mymin = np.min(a)
    min_positions = [i+1 for i, x in enumerate(a) if x == mymin]
    
    Alignments = []
    #Back tracing every possible alignment with lowest cost
    for start in min_positions:
        i = len(s1)
        align1, align2 = '', ''
        if start < len(s2):
            j = start
            while i > 0 and j > 0:
                score_current = score[i][j]
                score_diag = score[i-1][j-1]
                score_left = score[i][j-1]
                score_up = score[i-1][j]

                if score_current == score_diag + mch(s1[i-1], s2[j-1]):
                    a1,a2 = s1[i-1],s2[j-1]
                    i,j = i-1,j-1
                elif score_current == score_up + pt['gap']:
                    a1,a2 = s1[i-1],'-'
                    i -= 1
                elif score_current == score_left + pt['gap']:
                    a1,a2 = '-',s2[j-1]
                    j -= 1
                align1 += a1
                align2 += a2
            position = j
                   
            while i > 0:
                a1,a2 = s1[i-1],'-'
                align1 += a1
                align2 += a2
                i -= 1

            align1 = align1[::-1]
            align2 = align2[::-1]
            seqN = len(align1)
            sym = ''
            seq_score = 0

            #Calculate the local alignment cost 
            for i in range(seqN):
                a1 = align1[i]
                a2 = align2[i]
                if a1 == a2:
                    sym += a1
                    seq_score += mch(a1, a2)
                else: 
                    seq_score += mch(a1, a2)
                    sym += ' '     
            Alignments.append([position, seq_score])
    #Retrun a list of position and cost of every possible alignment for a read
    return Alignments

#Check the command line arguments

if len(sys.argv) < 4:
    print "Usage: <reference file (fasta)> <read file (fasta)> <out file>"
    sys.exit(0)

#Read the reference sequence and initiate the aligner
try:
    for s in SeqIO.parse(sys.argv[1], "fasta"):
        ref = s
        break #Stop after the fist sequence in the reference
except IOError as e:
    print "Could not read reference sequence file (see below)!"
    print e
    sys.exit(1)

output = open(sys.argv[3], 'w')
output.write("READ_NAME\tREF_NAME\tPOS\tSTRAND\tNUMBER_OF_ALIGNMENTS\tMAPPING_COST\n")

totalReads = 0  
match1 = 0
match2 = 0
match3 = 0
for read in SeqIO.parse(sys.argv[2], "fasta"):
    alignments = []
    
    score_f = fillMatrix(read.seq, ref.seq)
    loc_cost_f = [a for a in traceBack(score_f, read.seq, ref.seq)]
    for all_loc in loc_cost_f:
        alignments.append(Alignment(ref.name, all_loc[0], '+', all_loc[1]))

    score_r = fillMatrix(read.seq.reverse_complement(), ref.seq)
    loc_cost_r = [b for b in traceBack(score_r, read.seq.reverse_complement(), ref.seq)]
    for all_loc in loc_cost_r:
        alignments.append(Alignment(ref.name, all_loc[0], '-', all_loc[1]))
    
    totalReads += 1
    if len(alignments) > 1:
        match3 += 1
        sorted_align = sorted(alignments, key = lambda x: (x.getCost(), x.getPos()))
        output.write("%s\t%s\t%d\t%s\t%d\t%d\n" % 
            (read.name, sorted_align[0].getRefname(), sorted_align[0].getPos(), sorted_align[0].getStrand(), len(alignments), sorted_align[0].getCost()))

    elif len(alignments) == 1:
        match2 += 1
        output.write("%s\t%s\t%d\t%s\t%d\t%d\n" % 
            (read.name, sorted_align[0].getRefname(), sorted_align[0].getPos(), sorted_align[0].getStrand(), len(alignments), sorted_align[0].getCost()))
    else:
        match1 += 1
        output.write("%s\t*\t%s\t*%s\n" % 
                     (read.name))

output.close()

print "%d reads have been aligned to reference sequence" % totalReads
print "%d reads have failed to align to reference sequence %f %%" % (match1, (match1/totalReads)*100)
print "%d reads has matched exactly once to reference sequence %f %%" % (match2, (match2/totalReads)*100)
print "%d reads has matched more than once to reference sequence %f %%" % (match3, (match3/totalReads)*100)

print "Done"