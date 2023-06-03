"""
Utilities for spyglass
"""
import math
import numpy as np
import os
import pyfaidx
import random
import scipy.stats
import seqlogo
import sys
 
# Global vars
nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# -------------------- Error Handling --------------------
def ERROR(msg):
	"""
	Print an error message and die

	Parameters
	----------
	msg : str
	   Error message to print
	"""
	sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
	sys.exit(1)


# -------------------- Load sequences --------------------
def RetrieveFastaSeq(fasta, chromosome, start, end):
	"""
	Helper for LoadSeqs; find a specific region of given genome

	Parameters
	----------
	fasta : pyfaidx object 
	   pyfaidx object storing the faidx indexed reference sequence
	chromosome : str
	   chromosome of interest
	start : int
	   sequence start coordinate
	end : int
	   sequence end coordinate
	"""

	# return sequence on given chromosome from start coordinate to end coordinate
	return fasta[chromosome][(start - 1):end].seq	
	
def LoadSeqs(fasta, peakBed):
	"""
	Return a list of peak sequences specified in peaks file

	Parameters
	----------
	fasta : pyfaidx object 
	   pyfaidx object storing the faidx indexed reference sequence
	peakBed : str
	   BED-format file containing peak sequence regions
	"""

	seqs = []
	with open(peakBed, 'r') as pb:
		numPeaks = 0
		for line in pb:
			info = line.strip().split("\t")
			if len(info) != 6:
				ERROR("please check peak file format - see README.md for details")
			# append sequence on given chromosome beginning/ending at start/end
			seqToAdd = RetrieveFastaSeq(fasta, info[0], int(info[1]), int(info[2]))
			# convert to uppercase, if not already uppercase
			if seqToAdd.islower():
				seqToAdd = seqToAdd.upper()
			seqs.append(seqToAdd)
			numPeaks = numPeaks + 1
			seqLen = int(info[2]) - int(info[1]) + 1
	return seqs, numPeaks, seqLen

def GenerateRandomBkgSeqs(fasta, numSeqs, seqLen, log):
	"""
	Return a list of randomly generated background peak sequences from given reference genome

	Parameters
	----------
	fasta : pyfaidx object 
	   pyfaidx object storing the faidx indexed reference sequence
	numSeqs : int
	   number of background sequences
	seqLen : int
	   length of background sequences
	"""

	seqs = []
	chrs = list(fasta.keys())
	printcounter = 0
	for i in range(0, numSeqs):
		printcounter += 1
		if printcounter % 100 == 0:
			log.write("Generating background seq " + str(printcounter) + "/" + str(numSeqs))        
		# get a random chromosome
		chrom = np.random.choice(chrs, 1)
		# get a random start position on chosen chromosome
		start = random.randrange(0, len(fasta[chrom[0]][0:].seq) - seqLen)
		# append sequence on chrom beginning at start of lenth seqLen
		seqs.append(RetrieveFastaSeq(fasta, chrom[0], start, start + seqLen))         
	print(seqs)
	return seqs

	
# -------------------- Score sequences --------------------
def ScoreSeq(pwm, seq):
	"""
	Get the PWM score for a sequence

	Parameters
	----------
	pwm : 2d np.array
		position weight matrix
	seq : str
		sequence of nucleotides

	Returns
    -------
    score : float
       PWM score of seq
    """
	score = 0
	# Increment score by the corresponding A/C/T/G value for each position in the PWM
	for i in range(len(seq)):
		score = score + pwm[nucs[seq[i]]][i]
	return score
	
def ReverseComplement(seq):
	"""
	Get the reverse complement of a sequence

	Parameters
	----------
	seq : str
	   sequence of nucleotides
	
	Returns
    -------
    rev : str
       reverse complement of seq
	"""
	revcomp = ""
	revdict = {"A": "T", "C": "G", "G": "C", "T": "A"}
	# For each letter in seq, prepend its complement base to revcomp
	print(seq)
	for c in seq:
		revcomp = revdict.get(c) + revcomp
	return revcomp

def FindMaxScore(pwm, seq):
	"""
	Get the highest PWM match for a sequence
	[TODO: option to calculate without reverse complement]
	Parameters
	----------
	pwm : 2d np.array
		position weight matrix
	seq : str
		sequence of nucleotides

	Returns
    -------
    max_score : float
       top PWM score of seq
    """
	max_score = -1*np.inf
	n = pwm.shape[1]
	rev = ReverseComplement(seq)
	# Iterate through all n-length subseqs and compare the forward and reverse scores
	for i in range(len(seq)-n+1):
		max_score = max(max_score, ScoreSeq(pwm,seq[i:i+n]), ScoreSeq(pwm,rev[i:i+n]))
	return max_score


# -------------------- Set the threshold --------------------

def ComputeNucFreqs(sequences):
	"""
	Return freqs of ACGT

	Parameters
	----------
	sequences : str
	   sequence list

	Returns
	----------
	freqs : list of float
		frequencies of A, C, G, T in the sequences
	"""
	# default frequency distribution for A, C, G, T
	#freqs = [0.25, 0.25, 0.25, 0.25]
    # calculate specific frequencies
	counter = 0
	freqs = [0,0,0,0]
	for seq in sequences:
		counter += len(seq)
		for i, nuc in enumerate(["A", "C", "G", "T"]):
			freqs[i] += seq.count(nuc)
	freqs = [f / counter for f in freqs]
	return freqs
	
def RandomSequence(n, freqs):
	"""
	Generate a random sequence of length n with specified nucleotide frequencies

	Parameters
	----------
	n : int
	   length of sequence
	freqs : list of floats
	   list of frequencies [A, C, G, T]
	
	Returns
    -------
    seq : str
       random sequence
	"""
	seq = ""
	for i in range(n):
		seq += np.random.choice(["A", "C", "G", "T"], p = freqs)
	return seq

def GetThreshold(null_dist, pval):
	"""
	Score threshold for pvalue

	Parameters
	----------
	null_dist : list of floats
		scores null distribution
	pvalue : float
	   percentage of values that are above pvalue threshold

	Returns
	----------
	threshold : float
		threshold to achieve desired pvalue
	"""
	thresh = 0
	#print(null_dist)
	null_dist.sort(reverse = True)
	# set score threshold to obtain pvalue
	thresh = null_dist[int(len(null_dist) * pval)]
	return thresh


# -------------------- Test Enrichment --------------------	
def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
	"""
	Compute fisher exact test for whether motif is enriched in bound sequences

	Parameters
	----------
	peak_total : int
	   number of peaks
	peak_motif : int
	   number of peaks matching motif 
	bg_total : int
		number of background sequences
	bg_motif : int
		number of background sequences matching motif
	
	Returns
	----------
	pvalue : float
		fisher exact test pvalue
	"""
	pval = -1
	contingency_table = [[peak_motif, peak_total - peak_motif], [bg_motif, bg_total - bg_motif]]
	odds, pval = scipy.stats.fisher_exact(contingency_table)
	return pval