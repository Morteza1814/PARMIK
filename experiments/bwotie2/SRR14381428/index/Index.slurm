#!/bin/bash
#SBATCH --job-name="Index"
#SBATCH --error="Index.err"
#SBATCH --output="Index.out"
#SBATCH --mem=200G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bowtie2/bowtie2-2.5.4-linux-x86_64/bowtie2-build \
/p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381428/SRR14381428_reads \
/p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381428/index/SRR14381428
