#!/bin/bash
#SBATCH --job-name="SRR14381428"
#SBATCH --error="SRR14381428.err"
#SBATCH --output="SRR14381428.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado01
#SBATCH --partition=main
/u/rgq5aw/GIT/GenomePartialMatching/SRA_Toolkit/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR14381428 -p --fasta --outdir /bigtemp/rgq5aw/partialMatcherProject/datasets/metagenomicDB/SRR14381428/