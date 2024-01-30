#!/bin/bash
#SBATCH --job-name="Index"
#SBATCH --error="Index.err"
#SBATCH --output="Index.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado01
#SBATCH --partition=main
/p/lava/parmik/bwa/bwa index -p /p/lava/parmik/partialMatcherProject/tools/bwa/SRR14381425/index /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381425/SRR14381425_reads