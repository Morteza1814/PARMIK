#!/bin/bash
#SBATCH --job-name="SRR12432009_SarsGenome"
#SBATCH --error="SRR12432009_SarsGenome.err"
#SBATCH --output="SRR12432009_SarsGenome.out"
#SBATCH --mem=200G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bowtie2/bowtie2-2.5.4-linux-x86_64/bowtie2 -f -a -L 11 --local \
-x /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR12432009/index/SRR12432009 \
-U /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries \
-S /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR12432009/SarsGenome/alignment_L11_a_local

