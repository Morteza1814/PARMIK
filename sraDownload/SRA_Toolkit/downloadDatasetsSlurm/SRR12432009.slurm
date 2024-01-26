#!/bin/bash
#SBATCH --job-name="SRR12432009"
#SBATCH --error="SRR12432009.err"
#SBATCH --output="SRR12432009.out"
#SBATCH --mem=40G
#SBATCH --time=960
#SBATCH --nodelist=doppio04
#SBATCH --partition=nolim
/u/rgq5aw/GIT/PARMIK/sraDownload/SRA_Toolkit/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR12432009 -p --fasta --outdir /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/