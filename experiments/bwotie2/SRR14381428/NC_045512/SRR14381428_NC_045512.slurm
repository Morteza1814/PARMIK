#!/bin/bash
#SBATCH --job-name="SRR14381428_NC_045512"
#SBATCH --error="SRR14381428_NC_045512.err"
#SBATCH --output="SRR14381428_NC_045512.out"
#SBATCH --mem=80G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bowtie2/bowtie2-2.5.4-linux-x86_64/bowtie2 -f -a -L 11 --local \
-x /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381428/index/SRR14381428 \
-U /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-S /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381428/NC_045512/alignment_L11_a

