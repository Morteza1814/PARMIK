#!/bin/bash
#SBATCH --job-name="Index"
#SBATCH --error="Index.err"
#SBATCH --output="Index.out"
#SBATCH --mem=200G
#SBATCH --time=960
#SBATCH --nodelist=hydro
#SBATCH --partition=main
/p/lava/parmik/bwa/bwa index -p /p/lava/parmik/partialMatcherProject/tools/bwa/SRR12432009/index /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads