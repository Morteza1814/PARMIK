#!/bin/bash
#SBATCH --job-name="SRR14381425"
#SBATCH --error="SRR14381425.err"
#SBATCH --output="SRR14381425.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado07
#SBATCH --partition=main
/u/rgq5aw/GIT/GenomePartialMatching/SRA_Toolkit/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR14381425 -p --fasta --outdir /bigtemp/rgq5aw/partialMatcherProject/datasets/metagenomicDB/SRR14381425/