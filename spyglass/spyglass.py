#!/usr/bin/env python

"""
Command-line script to identify enriched motifs in genomic regions

Similar to Homer findMotifsGenome.pl and Meme Suite SEA
"""

import argparse
from . import myutils as myutils
from spyglass import __version__

def main():

    # create parser object
    parser = argparse.ArgumentParser(
        prog = "spyglass",
        description = "Command-line script to identify enriched motifs in genomic regions"
    )

    # required input
    parser.add_argument("reference", help = "reference genome in fasta format", type = str)
    parser.add_argument("peaks", help = "genomic regions file", type = str)
    parser.add_argument("motifs", help = "concatenated PWMs of all motifs of interest", type = str)

    # output
    parser.add_argument("-o", "--output", help = "write output to file. Default: stdout", metavar = "FILE", type = str, required = False)

    # optional options
    parser.add_argument("-b", "--background", help = "comma-separated custom frequencies for A/T/C/G. Default: 0.25,0.25,0.25,0.25", type = str, metavar = "BACKGROUND", required = False)
    parser.add_argument("--version", help = "print the version and quit", action = "version", version = '{version}'.format(version = __version__))

    # parse arguments
    args = parser.parse_args()

    # set up output file
    if args.out is None:
        outf = sys.stdout
    else:
        outf = open(args.out, "w")

    # load fasta file
    if args.reference is not None:
        if not os.path.exists(args.reference):
            myutils.ERROR("{fasta} does not exist".format(fasta = args.reference))
        reffasta = pyfaidx.Fasta(args.reference)
    else:
        reffasta = None
        myutils.ERROR("please specify a fasta file")
    






    outf.close()
    sys.exit(0)

if __name__ == "__main__":
    main()