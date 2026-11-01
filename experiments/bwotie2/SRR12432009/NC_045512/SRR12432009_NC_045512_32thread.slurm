#!/bin/bash
#SBATCH --job-name="SRR12432009_NC_045512_32thread"
#SBATCH --error="SRR12432009_NC_045512_32thread.err"
#SBATCH --output="SRR12432009_NC_045512_32thread.out"
#SBATCH --mem=80G
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
/p/lava/parmik/bowtie2/bowtie2-2.5.4-linux-x86_64/bowtie2 -f -a -L 11 -p 32 --local \
-x /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR12432009/index/SRR12432009 \
-U /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-S /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR12432009/NC_045512/alignment_L11_a_local

