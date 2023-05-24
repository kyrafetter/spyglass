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
	   pyfaidx object storing the 
	chromosome : str
	   chromosome of interest
	start : int
	   p1 description
	end : int
	   p1 description
	"""

	return fasta[chromosome][(start - 1):end].seq
	
	
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
