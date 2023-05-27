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
	Description

	Parameters
	----------
	p1 : name
	   p1 description
	"""
	# CODE HERE
	
def ReverseComplement(seq):
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

def FindMaxScore(pwm, seq):
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

# -------------------- Set the threshold --------------------

def ComputeNucFreqs(fasta, seq):
	"""
	Return freqs of ACGT

	Parameters
	----------
	fasta : str
	   given fasta sequence to compute frequencies from
	seq : str
	   sequence list

	Returns
	----------
	freqs : list of float
		frequencies of A, C, G, T in the sequences
	"""
	# CODE HERE
	# default frequency distribution for A, C, G, T
	freqs = [0.25, 0.25, 0.25, 0.25]
	# calculate specific frequencies 
	counter = 0
	A_dist = 0
	C_dist = 0
	G_dist = 0
	T_dist = 0
	for i in range(len(seq)):
		counter += 1
		if(seq.get(i) == "A"):
			A_dist += 1
		elif(seq.get(i) == "C"):
			C_dist += 1
		elif(seq.get(i) == "G"):
			G_dist += 1
		else 
			T_dist += 1	
	freqs[0] = A_dist 



	

# don't need to do this one	
def RandomSequence(n, seq):
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



def GetThreshold(pvalue):
	"""
	Score threshold for pvalue

	Parameters
	----------
	pvalue : float
	   percentage of values that are above pvalue threshold

	Returns
	----------
	threshold : float
		threshold to achieve desired pvalue
	"""
	# CODE HERE
	


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
	# CODE HERE
