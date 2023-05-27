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
import datetime
import numpy as np

def main():

    # create parser object
    parser = argparse.ArgumentParser(
        prog = "spyglass",
        description = "Command-line script to identify enriched motifs in genomic regions"
    )

    # required input
    parser.add_argument("reference", help = "faidx indexed reference sequence file", type = str)
    parser.add_argument("peaks", help = "BED file of genomic peak regions", type = str)
    parser.add_argument("motifs", help = "PWM file of motif PWMs of interest", type = str)

    # output
    parser.add_argument("-o", "--output", help = "write output to file. Default: stdout", metavar = "FILE", type = str, required = False)
    parser.add_argument("-l", "--logfile", help = "write log to file. Default: stderr", metavar = "FILE", type = str, required = False)

    # optional options
    parser.add_argument("-b", "--background", help = "BED file of user-specified background genomic peak regions. Default: random background implemented", type = str, metavar = "BACKGROUND", required = False)
    parser.add_argument("-p", "--pval", help = "p-value threshold for significant enrichment. Default: 0.0000001", type = float, metavar = "PVALUE", required = False)
    parser.add_argument("-r", "--reverse", help = "consider reverse complement in enrichment analysis. Default: True", type = bool, metavar = "REVERSE", required = False)
    parser.add_argument("-s", "--seqlogo", help = "generate the motif logo(s) of enriched motifs. Default: True", type = bool, metavar = "SEQLOGO", required = False)
    parser.add_argument("--version", help = "print the version and quit", action = "version", version = '{version}'.format(version = __version__))

    # parse arguments
    args = parser.parse_args()

    # set up output file
    if args.output is None:
        outf = sys.stdout
    else:
        outf = open(args.output, "w")
    
    # set up log file
    if args.log is None:
        log = sys.stderr
    else:
        log = open(args.log, "w")
    log.write("Welcome to spyglass!\n")
    log.write("Start time: ", datetime.datetime.now())
    log.write("\n\n")


    # -------------------- Load/Process --------------------

    # load fasta file
    if args.reference is not None:
        if not os.path.exists(args.reference):
            myutils.ERROR("{fasta} does not exist".format(fasta = args.reference))
        reffasta = pyfaidx.Fasta(args.reference)
        log.write("Using fasta: {fasta}".format(fasta = args.reference))
        log.write("\n")
    else:
        myutils.ERROR("please specify a fasta file")
    
    # load pwm file
    PWMList = []
    pwm_names = []
    if args.motifs is not None:
        if not os.path.exists(args.motifs):
            myutils.ERROR("{motifs} does not exist".format(motifs = args.motifs))
        m = open(args.motifs, "r")
        currPWM = []
        for line in m:
            if line[0] == ">":
                if len(currPWM) != 0:
                    PWMList.append(np.array(currPWM).transpose())
                currPWM = []
                pwm_names.append(line.strip().split(">")[0])
            else:
                weights = line.strip().split("\t")
                currPWM.append(weights)
        m.close()
        log.write("Using motifs: {motifs}".format(motifs = args.motifs))
        log.write("\n")
    else:
        myutils.ERROR("please specify a motifs pwm file")

    # load foreground bed file
    if args.peaks is not None:
        if not os.path.exists(args.peaks):
            myutils.ERROR("{peaks} does not exist".format(peaks = args.peaks))
        peak_seqs, numPeaks, seqLen = myutils.LoadSeqs(reffasta, args.peaks)
        log.write("Using peaks: {peaks}".format(peaks = args.peaks))
        log.write("\n")
    else:
        myutils.ERROR("please specify a foreground peaks bed file")
    
    # load background bed file
    if args.background is not None:
        if not os.path.exists(args.background):
            myutils.ERROR("{background} does not exist".format(background = args.background))
        bg_seqs, numPeaks, seqLen = myutils.LoadSeqs(reffasta, args.background)
        log.write("Using custom background: {background}".format(background = args.background))
        log.write("\n\n")
    else:
        bg_seqs = myutils.GenerateRandomBkgSeqs(reffasta, numPeaks, seqLen)
        log.write("Generating random background")
        log.write("\n\n")


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
    log.close()
    sys.exit(0)

if __name__ == "__main__":
    main()