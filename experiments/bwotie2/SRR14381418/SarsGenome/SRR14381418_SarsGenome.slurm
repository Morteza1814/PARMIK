#!/bin/bash
#SBATCH --job-name="SRR14381418_SarsGenome"
#SBATCH --error="SRR14381418_SarsGenome.err"
#SBATCH --output="SRR14381418_SarsGenome.out"
#SBATCH --mem=80G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bowtie2/bowtie2-2.5.4-linux-x86_64/bowtie2 -f -a -L 11 --local \
-x /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381418/index/SRR14381418 \
-U /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries \
-S /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381418/SarsGenome/alignment_L11_a

