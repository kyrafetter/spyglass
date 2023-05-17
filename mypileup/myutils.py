"""
Utilities for mypileup
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

def GetQual(qualscore):
	"""
	Convert a quality score to ASCII value
	Use encoding chr(qualscore+33)

	Parameters
	----------
	qualscore : int
	   Numerical quality score

	Returns
	-------
	qualchar : str
	   Quality score ascii character
	"""
	qualchar = chr(qualscore+33)
	return qualchar

def GetReadBase(pileupread, refbase):
	"""
	Extract the mpileup base for a read
	from a single read at a particular position

	Mimics the format of "read bases" specified here:
	http://www.htslib.org/doc/samtools-mpileup.html

	Note: currently does not handle indels

	Parameters
	----------
	pileupread : pysam.PileupRead
	   The pileupread object from which 
	   to extract info
	refbase : str
	   Reference base for the position

	Returns
	-------
	nucstring : str
	   Representation of this read in
	   mpileup format
	"""
	# Extract read base
	read_base = pileupread.alignment.query_sequence[pileupread.query_position]

	# Set the base
	nucstring = ""
	if read_base == refbase:
		if pileupread.alignment.is_reverse:
			nucstring = ","
		else: nucstring = "."
	else:
		if pileupread.alignment.is_reverse:
			nucstring = read_base.lower()
		else: nucstring = read_base

	# Check if it is the first or last base of the read
	if len(pileupread.alignment.query_sequence) == (pileupread.query_position-1):
		nucstring += "$"
	if pileupread.query_position == 0:
		nucstring = "^" + GetQual(pileupread.alignment.mapping_quality) + nucstring

	return nucstring

def GetReadQual(pileupread):
	"""
	Extract the mpileup quality score for a read
	from a single read at a particular position

	Parameters
	----------
	pileupread : pysam.PileupRead
	   The pileupread object from which 
	   to extract info

	Returns
	-------
	qualstring : str
	   Quality score ascii value	
	"""
	qualscore = pileupread.alignment.query_qualities[pileupread.query_position]
	qualstring = GetQual(qualscore) 
	return qualstring
