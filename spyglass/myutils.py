"""
Utilities for spyglass
"""
import math
import numpy as np
import pandas
import pyfaidx
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

def RetrieveFastaSeq(fasta, chromosome, start, end):
	"""
	Helper for LoadSeqs; find a specific region of given genome

	Parameters
	----------
	p1 : name
	   p1 description
	p1 : name
	   p1 description
	"""
	# CODE HERE
	
def LoadSeqs(fasta, peakBed, bgBed):
	"""
	Return a list of peak sequences specified in peaks file and a list of bg seqs if specified or randomly generated bg seqs

	Parameters
	----------
	p1 : name
	   p1 description
	p1 : name
	   p1 description
	"""
	# CODE HERE

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
		score += pwm[nucs.get(seq[i],i)]
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
def ComputeNucFreqs(fasta, seq):
	"""
	Description

	Parameters
	----------
	p1 : name
	   p1 description
	p1 : name
	   p1 description
	"""
	# CODE HERE
	
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
		p = random.uniform(0,1)
		char = "T"
		if p <= freqs[0]:
			char = "A"
		elif(p <= freqs[0]+freqs[1]):
			char = "C"
		elif(p <= freqs[0]+freqs[1]+freqs[2]):
			char = "G"
		seq += char
	return seq
	
def GetThreshold(pval):
	"""
	Description

	Parameters
	----------
	p1 : name
	   p1 description
	p1 : name
	   p1 description
	"""
	# CODE HERE
	
# -------------------- Test Enrichment --------------------	
def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
	"""
	Description

	Parameters
	----------
	p1 : name
	   p1 description
	p1 : name
	   p1 description
	"""
	# CODE HERE
