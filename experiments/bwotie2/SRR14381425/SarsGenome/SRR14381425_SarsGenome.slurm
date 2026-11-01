#!/bin/bash
#SBATCH --job-name="SRR14381425_SarsGenome"
#SBATCH --error="SRR14381425_SarsGenome.err"
#SBATCH --output="SRR14381425_SarsGenome.out"
#SBATCH --mem=80G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bowtie2/bowtie2-2.5.4-linux-x86_64/bowtie2 -f -a -L 11 --local \
-x /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381425/index/SRR14381425 \
-U /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries \
-S /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381425/SarsGenome/alignment_L11_a

