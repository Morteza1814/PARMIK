#!/bin/bash
#SBATCH --job-name="SRR14381434"
#SBATCH --error="SRR14381434.err"
#SBATCH --output="SRR14381434.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado08
#SBATCH --partition=main
/u/rgq5aw/GIT/PARMIK/sraDownload/SRA_Toolkit/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR14381434 -p --fasta --outdir /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/