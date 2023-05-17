#!/usr/bin/env python

"""
Command-line script to perform  pileup of BAM files

Similar to samtools mpileup
"""

import argparse
from . import myutils as myutils
from mypileup import __version__
import os
import pyfaidx
import pysam
import sys

def main():
	parser = argparse.ArgumentParser(
		prog="mypileup",
		description="Command-line script to perform  pileup of BAM files"
	)

	# Input
	parser.add_argument("bam", help="Indexed BAM file", \
		type=str)

	# Output
	parser.add_argument("-o", "--out", help="Write output to file. " \
		"Default: stdout", metavar="FILE", type=str, required=False)

	# Other options
	parser.add_argument("-f", "--fasta-ref", \
		help="faidx indexed reference sequence file", \
		type=str, metavar="FILE", required=False)
	parser.add_argument("-r", "--region", help="region in which pileup is " \
		"generated. Format chr:start-end", \
		type=str, metavar="REG", required=False)
	parser.add_argument("--version", help="Print the version and quit", \
		action="version", version = '{version}'.format(version=__version__))

	# Parse args
	args = parser.parse_args()

	# Set up output file
	if args.out is None:
		outf = sys.stdout
	else: outf = open(args.out, "w")

	# Load FASTA
	if args.fasta_ref is not None:
		if not os.path.exists(args.fasta_ref):
			myutils.ERROR("{fasta} does not exist".format(fasta=args.fasta_ref))
		reffasta = pyfaidx.Fasta(args.fasta_ref)
	else:
		reffasta = None

	# Load BAM
	bamfile = pysam.AlignmentFile(args.bam, "rb")

	if args.region is not None:
		region = args.region
	else:
		region = None

	# Peform pileup
	for pileupcolumn in bamfile.pileup(region=region):
		chrom = pileupcolumn.reference_name
		position = pileupcolumn.reference_pos
		numreads = 0
		if reffasta is not None:
			refbase = str(reffasta[chrom][position])
		else:
			refbase = "N"
		read_bases = []
		read_quals = []
		for pileupread in pileupcolumn.pileups:
			if not pileupread.is_del and not pileupread.is_refskip:
			# query position is None if is_del or is_refskip is set.
				read_bases.append(myutils.GetReadBase(pileupread, refbase))
				read_quals.append(myutils.GetReadQual(pileupread))
				numreads += 1
		outf.write("\t".join([chrom, str(position+1), refbase, str(numreads), \
				"".join(read_bases), "".join(read_quals)])+"\n")
	bamfile.close()
	outf.close()
	sys.exit(0)

if __name__ == "__main__":
    main()
