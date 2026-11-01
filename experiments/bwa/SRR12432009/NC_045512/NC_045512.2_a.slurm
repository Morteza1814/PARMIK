#!/bin/bash
#SBATCH --job-name="NC_045512.2_a"
#SBATCH --error="NC_045512.2_a.err"
#SBATCH --output="NC_045512.2_a.out"
#SBATCH --mem=100G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bwa/bwa mem /p/lava/parmik/partialMatcherProject/tools/bwa/SRR12432009/index/index /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -a > /p/lava/parmik/partialMatcherProject/tools/bwa/SRR12432009/NC_045512/NC_045512.2_SRR12432009_a.sam