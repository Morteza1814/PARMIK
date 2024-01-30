#!/bin/bash
#SBATCH --job-name="Index"
#SBATCH --error="Index.err"
#SBATCH --output="Index.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado02
#SBATCH --partition=main
/p/lava/parmik/bwa/bwa index -p /p/lava/parmik/partialMatcherProject/tools/bwa/SRR14381434/index /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/SRR14381434_reads