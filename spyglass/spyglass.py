#!/usr/bin/env python

"""
Command-line script to identify enriched motifs in genomic regions

Similar to Homer findMotifsGenome.pl and Meme Suite SEA
"""

import argparse
from . import myutils as myutils
from spyglass import __version__
import os
import pyfaidx
import sys

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

    # -------------------- Load/Process --------------------

    # load fasta file
    if args.reference is not None:
        if not os.path.exists(args.reference):
            myutils.ERROR("{fasta} does not exist".format(fasta = args.reference))
        reffasta = pyfaidx.Fasta(args.reference)
    else:
        reffasta = None
        myutils.ERROR("please specify a fasta file")
    
    # load pwm file
    PWMList = []
    # load bed file(s)
    peak_seqs = LoadSeqs()
    bg_seqs = Load or generate

    # -------------------- Motif Enrichment --------------------
    reverse_seqs = [myutils.ReverseComplement(item) for item in peak_seqs] + [myutils.ReverseComplement(item) for item in bg_seqs]
    
    # initialize vars for null dist sim
    freqs = myutils.ComputeNucFreqs(peak_seqs+bg_seqs+reverse_seqs)
    numsim = 10000
    null_pval = 0.01
    enriched_pval = args.pval if args.pval is not None else 0.0000001

    for i in range(len(PWMList)):
        log.write("Motif ", i, ": ")
        # generate null dist of random seqs
        null_scores = [myutils.ScoreSeq(PWMList[i], myutils.RandomSequence(PWMList[i].shape[1], freqs)) for j in range(numsim)]
        # compute score significance threshold
        log.write("Starting thresholding... ")
        thresh = myutils.GetThreshold(null_scores, null_pval)
        # get number of matches above threshold
        log.write("finding matches... ")
        num_peak_pass = np.sum([int(myutils.FindMaxScore(pwm, seq)>thresh) for seq in peak_seqs])
        num_bg_pass = np.sum([int(myutils.FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
        # perform fisher exact test
        log.write("performing fisher's exact test...\n")
        fisher_pval = myutils.ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
        enriched = "yes" if fisher_pval < enriched_pval else "no"
        # output
        outf.write("\t".join([pwm_names[i], num_peak_pass+"/"+len(peak_seqs), num_bg_pass+"/"+len(bg_seqs), fisher_pval, enriched]))
    log.write("End time: ", datetime.datetime.now())

    log.close()
    outf.close()
    sys.exit(0)

if __name__ == "__main__":
    main()