#!/usr/bin/env python

"""
Command-line script to identify enriched motifs in genomic regions

Similar to Homer findMotifsGenome.pl and Meme Suite SEA
"""

from . import myutils as myutils
from spyglass import __version__
import argparse
import datetime
import numpy as np
import os
import pandas as pd
import pyfaidx
import seqlogo
import sys

def main():

    # -------------------- Parse Command Line --------------------

    # create parser object
    parser = argparse.ArgumentParser(
        prog = "spyglass",
        description = "Command-line script to identify enriched motifs in genomic regions"
    )

    # positional (required) arguments
    parser.add_argument("peaks", help = "BED file of genomic peak regions", type = str)
    parser.add_argument("reference", help = "faidx indexed reference sequence file in Fasta format", type = str)
    parser.add_argument("motifs", help = "PWM file of motif PWMs of interest", type = str)

    # optional arguments
    parser.add_argument("-o", "--output", help = "write output to directory. Default: working directory", metavar = "DIRECTORY", type = str, required = False)
    parser.add_argument("-l", "--log", help = "write log to file. Default: stdout", metavar = "FILE", type = str, required = False)
    parser.add_argument("-b", "--background", help = "BED file of user-specified background genomic peak regions. Default: background sequences will be chosen randomly from the reference genome", type = str, metavar = "BACKGROUND", required = False)
    parser.add_argument("-p", "--pval", help = "p-value threshold for significant enrichment. Default: 0.0000001", type = float, metavar = "PVALUE", required = False)
#     parser.add_argument("-r", "--not-reverse", help = "do not consider reverse complement in enrichment analysis. Default: FALSE", action = 'store_false', required = False)
    parser.add_argument("-nr", "--not-reverse", help = "do not consider reverse complement in enrichment analysis. Default: FALSE", nargs='?', const = True, default = False, required = False)
#     parser.add_argument("-s", "--seqlogo", help = "generate the motif logo(s) of enriched motifs. Default: FALSE", type = bool, metavar = "SEQLOGO", required = False)
    parser.add_argument("-s", "--seqlogo", help = "generate the motif logo(s) of enriched motifs. Default: FALSE", nargs='?', const = True, default = False, required = False)
    parser.add_argument("--version", help = "print the version and quit", action = "version", version = '{version}'.format(version = __version__))

    # parse arguments
    args = parser.parse_args()

    # -------------------- Set-Up Log and Outfile --------------------

    # set up output file
    if args.output is None:
        outdir = os.getcwd()
    else:
        outdir = args.output
    
    if outdir[-1] != "/":
        outdir = outdir + "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outf = open(outdir + "spyglass_results.tsv", "w")
    
    # set up log file
    if args.log is None:
        log = sys.stdout
    else:
        log = open(args.log, "w")

    log.write("Welcome to spyglass!\n")
    log.write("Start time: ")
    log.write(str(datetime.datetime.now()))
    log.write("\n\n")


    # -------------------- Load/Process Input Files --------------------

    log.write("Loading input files:\n")

    # load fasta file
    if args.reference is not None:
        if not os.path.exists(args.reference):
            myutils.ERROR("{fasta} does not exist".format(fasta = args.reference))
        try:
            reffasta = pyfaidx.Fasta(args.reference)
        except Exception:
            myutils.ERROR("please check fasta file format - see README.md for details")
        log.write("Using fasta: {fasta}".format(fasta = args.reference))
        log.write("\n")
    else:
        myutils.ERROR("please specify a fasta file")
    
    # load pwm file
    PWMList = []
    pwm_names = []
    if args.motifs is not None:
        # add pretty error later
        if not os.path.exists(args.motifs):
            myutils.ERROR("{motifs} does not exist".format(motifs = args.motifs))
        m = open(args.motifs, "r")
        currPWM = []
        for line in m:
            if line[0] == ">":
                if len(currPWM) != 0:
                    PWMList.append(np.array(currPWM).transpose())
                currPWM = []
                pwm_names.append(line.strip().split(">")[1])
            else:
                weights = line.strip().split("\t")
                weights_floats = [float(i) for i in weights]
                currPWM.append(weights_floats)
        PWMList.append(np.array(currPWM).transpose())
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
        log.write("\n\n")
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
        log.write("Generating random background from reference genome...\n")
        log.write("0   50   100 (% complete)\n")
        bg_seqs = myutils.GenerateRandomBkgSeqs(reffasta, numPeaks, seqLen, log)
        log.write("Done\n\n")
    

    # -------------------- Perform Motif Enrichment --------------------

       
    log.write("Starting motif enrichment analysis...\n\n")

    # calculate frequencies
    if args.not_reverse:
        log.write("note: not using reverse complement strand\n\n")
        freqs = myutils.ComputeNucFreqs(peak_seqs+bg_seqs)
    else:
        reverse_seqs = [myutils.ReverseComplement(item) for item in peak_seqs] + [myutils.ReverseComplement(item) for item in bg_seqs]
        freqs = myutils.ComputeNucFreqs(peak_seqs+bg_seqs+reverse_seqs)
    
    # initialize vars for null dist sim
    numsim = 10000
    null_pval = 0.01
    enriched_pval = args.pval if args.pval is not None else 0.0000001
    enrichment_results = []
    
    # for each motif in pwm 
    for i in range(len(PWMList)):
        log.write("Motif " + pwm_names[i] + ":\n")
        
        # generate null dist of random seqs
        null_scores = [myutils.ScoreSeq(PWMList[i], myutils.RandomSequence(PWMList[i].shape[1], freqs)) for j in range(numsim)]
        
        # compute score significance threshold
        log.write("Starting thresholding...")
        thresh = myutils.GetThreshold(null_scores, null_pval)
        log.write("Done\n")
        
        # get number of matches above threshold
        log.write("Finding motif matches in foreground and background...")
        num_peak_pass = np.sum([int(myutils.FindMaxScore(PWMList[i], seq) > thresh) for seq in peak_seqs])
        num_bg_pass = np.sum([int(myutils.FindMaxScore(PWMList[i], seq) > thresh) for seq in bg_seqs])
        log.write("Done\n")
        
        # perform fisher exact test
        log.write("Performing fisher's exact test...")
        fisher_pval = myutils.ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
        enriched = "yes" if fisher_pval < enriched_pval else "no"
        log.write("Done\n\n")
        
        # output
        enrichment_results.append([pwm_names[i], str(num_peak_pass) + "/" + str(len(peak_seqs)), str(num_bg_pass) + "/" + str(len(bg_seqs)), "{0:.{1}e}".format(fisher_pval, 3), enriched])

    # -------------------- Summarize Results --------------------

    log.write("Summarizing results...")

    if outf == sys.stdout:
        log.write("please see below\n\n")
    else:
        log.write("please see: " + outdir + "\n")

    results = pd.DataFrame(enrichment_results)
    results.columns = ["# motif_name", "motif_occurences_foreground", "motif_occurences_background", "pvalue", "enriched?"]
    results.sort_values(by = ["pvalue", "# motif_name"], ascending = False).reset_index(drop = True)
    results.to_csv(outf, sep = "\t", header = True, index = False)


    # -------------------- Generate SeqLogo --------------------
    
    # 6.03 - work in progress; does not run at this time
    #if args.seqlogo is None:
    if False:
        log.write("\nGenerating seqlogo...")
        i = 0
        for pwm in PWMList:
            print(pwm)
            print(seqlogo.Pwm(pwm.transpose()))
            #logo = seqlogo.seqlogo(seqlogo.pwm2ppm(seqlogo.Pwm(pwm)), ic_scale = True, format = 'png', size = 'medium')
            #logo.save(outdir + pwm_names[i] + "_logo.png")
            i = i + 1
        log.write("Done\n")


    # -------------------- Close Log and Outfiles and Exit --------------------

    log.write("\nCongrats! spyglass was successfully executed! Cheers!\n")
    log.write("End time: ")
    log.write(str(datetime.datetime.now()))
    log.write("\n")

    log.close()
    outf.close()
    sys.exit(0)

if __name__ == "__main__":
    main()