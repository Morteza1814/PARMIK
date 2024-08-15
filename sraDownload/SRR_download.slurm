#!/bin/bash
#SBATCH --job-name="SRR12432009"
#SBATCH --error="SRR12432009.err"
#SBATCH --output="SRR12432009.out"
#SBATCH --mem=40G
#SBATCH --time=1-00:00:00
#SBATCH --partition=cpu
/p/lava/parmik/SRA_Toolkit/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR14381418 -p --fasta --outdir /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381418/
/p/lava/parmik/SRA_Toolkit/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR14381422 -p --fasta --outdir /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/
