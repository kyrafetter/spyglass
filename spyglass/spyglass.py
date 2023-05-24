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
    all_seqs = concat (peak_seqs / bg seqs/ rev comp for both)
    freqs = ComputeNucFreqs(all_seqs)
    numsim = 10000 # can change
    pval = 0.01 # default - option?

    for i in range(len(PWMList)):
        null_scores = [ScoreSeq(PWMList[i], RandomSequence(PWMList[i].shape[1], freqs)) for j in range(numsim)]
        thresh = GetThreshold(null_scores, pval)
        num_peak_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in peak_seqs])
        num_bg_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
        enriched_pval = ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
        print stuff


    # -------------------- Notes from lab exercises -------------------- 
    freqs = ComputeNucFreqs(peak_seqs+bg_seqs+[ReverseComplement(item) for item in peak_seqs] + [ReverseComplement(item) for item in bg_seqs])
    numsim = 10000
    pval = 0.01

    pwm_thresholds = []
    fig = plt.figure()
    fig.set_size_inches((10, 4))
    for i in range(3):
        null_scores = [ScoreSeq(PWMList[i], RandomSequence(PWMList[i].shape[1], freqs)) for j in range(numsim)]
        thresh = GetThreshold(null_scores, pval)
        ax = fig.add_subplot(1, 3, i+1)
        ax.hist(null_scores, bins=10);
        ax.axvline(x=thresh, color="red")
        ax.set_xlabel("Score")
        ax.set_ylabel("Frequency");
        ax.set_title(pwm_names[i])
        pwm_thresholds.append(thresh)
    fig.tight_layout()

    for i in range(len(PWMList)):
        pwm = PWMList[i]
        thresh = pwm_thresholds[i]
        num_peak_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in peak_seqs])
        num_bg_pass = np.sum([int(FindMaxScore(pwm, seq)>thresh) for seq in bg_seqs])
        pval = ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
        print("PWM: %s, %s/%s peaks, %s/%s background; p-val: %s"%(pwm_names[i], num_peak_pass, len(peak_seqs), num_bg_pass, len(bg_seqs), pval))
        





    outf.close()
    sys.exit(0)

if __name__ == "__main__":
    main()