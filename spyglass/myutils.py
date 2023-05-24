"""
Utilities for spyglass
"""
import math
import numpy as np
import os
import pandas
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
		for line in pb:
			info = line.strip().split("\t")
			# append sequence on given chromosome beginning/ending at start/end
			seqs.append(RetrieveFastaSeq(fasta, info[0], info[1], info[2]))
	return seqs

def GenerateRandomBkgSeqs(fasta, numSeqs, seqLen):
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
	chrs = fasta.keys()
	for i in range(0, numSeqs):
		# get a random chromosome
		chrom = np.random.choice(chrs, 1)
		# get a random start position on chosen chromosome
		start = random.randrange(1, len(fasta[chrom].seq) - seqLen)
		# append sequence on chrom beginning at start of lenth seqLen
		seqs.append(RetrieveFastaSeq(fasta, chrom, start, start + seqLen))
	return seqs

	

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
	Description

	Parameters
	----------
	p1 : name
	   p1 description
	p1 : name
	   p1 description
	"""
	# CODE HERE
	
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
