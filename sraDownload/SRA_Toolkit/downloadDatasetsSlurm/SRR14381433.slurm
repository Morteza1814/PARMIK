#!/bin/bash
#SBATCH --job-name="SRR14381433"
#SBATCH --error="SRR14381433.err"
#SBATCH --output="SRR14381433.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado06
#SBATCH --partition=main
/u/rgq5aw/GIT/GenomePartialMatching/SRA_Toolkit/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR14381433 -p --fasta --outdir /bigtemp/rgq5aw/partialMatcherProject/datasets/metagenomicDB/SRR14381433/