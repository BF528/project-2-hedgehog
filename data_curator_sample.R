module load sratoolkit

#!/bin/bash

#$ -P bf528
#$ -cwd
#$ -j y
#$ -pe mpi_16_tasks_per_node 16

echo "Running job $JOB_ID"
echo "Started: $(date +%F)"
echo "Running in directory: $PWD"
# your commands here

module load fastqc
fastq-dump -I --split-files P0_1.sra
-o --outdir  /hedgehog/project_2/data_curator/
fastqc -t 4 --outdir /projectnb/bf528/users/hedgehog/project_2/data_curator/quality/ P0_1_1.fastq
fastqc -t 4 --outdir /projectnb/bf528/users/hedgehog/project_2/data_curator/quality/ P0_1_2.fastq
echo "Job finished: $(date +%F)"
