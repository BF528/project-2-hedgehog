#!/bin/bash

# CuffDiff for two P0 samples and two Ad samples

P0BAM=<path_to_P0.bam>

SAMPLEDIR=<path_to_sample_dir>
P0REPS=$P0BAM,$SAMPLEDIR/P0_2/accepted_hits.bam
ADREPS=$SAMPLEDIR/Ad_1/accepted_hits.bam,$SAMPLEDIR/Ad_2/accepted_hits.bam

LABEL="P0,Ad"
OUTDIR=cuffdiff_out
FASTA=<path_to_reference_fasta>

cuffdiff -p 16 -L $LABEL -u -b $FASTA -o $OUTDIR $SAMPLEDIR/merged_asm/merged.gtf $P0REPS $ADREPS
