#! /bin/bash

# Paired-end sequencing alignment
tophat -r 200 -G <path_to_ref_GTF> --segment-length=20 --segment-mismatches=1 --no-novel-juncs -o P0_1_tophat -p 16 <path_to_ref_seq> <path_to_fastq1> <path_to_fastq2>
