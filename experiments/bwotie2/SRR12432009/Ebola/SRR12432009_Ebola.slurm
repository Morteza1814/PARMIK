#!/bin/bash
#SBATCH --job-name="SRR12432009_Ebola"
#SBATCH --error="SRR12432009_Ebola.err"
#SBATCH --output="SRR12432009_Ebola.out"
#SBATCH --mem=80G
#SBATCH --time=1-00:00:00
#SBATCH --partition=gpu
/p/lava/parmik/bowtie2/bowtie2-2.5.4-linux-x86_64/bowtie2 -f -a -L 11 --local \
-x /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR12432009/index/SRR12432009 \
-U /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/ebola_queries \
-S /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR12432009/Ebola/alignment_L11_a_local

