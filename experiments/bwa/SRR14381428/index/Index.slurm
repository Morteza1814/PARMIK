#!/bin/bash
#SBATCH --job-name="Index"
#SBATCH --error="Index.err"
#SBATCH --output="Index.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cheetah02
#SBATCH --partition=gpu
/p/lava/parmik/bwa/bwa index -p /p/lava/parmik/partialMatcherProject/tools/bwa/SRR14381428/index /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381428/SRR14381428_reads