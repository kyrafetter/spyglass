"""
Utilities for spyglass
"""
import sys

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
	Helper for LoadSeqs: returns the sequence at specified coordinates in the ref genome

	Parameters
	----------
	fasta : pyfaidx Fasta object
	   pyfaidx object storing the reference genome
	chromosome : str
	   chromosome of interest
	start : int
	   start coordinate of sequence
	end : int
	   end coordinate of sequence
	"""

	
	
	
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
