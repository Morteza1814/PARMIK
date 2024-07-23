#!/bin/bash
#SBATCH --job-name="SarsGenome"
#SBATCH --error="SarsGenome.err"
#SBATCH --output="SarsGenome.out"
#SBATCH --mem=200G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bwa/bwa mem /p/lava/parmik/partialMatcherProject/tools/bwa/SRR12432009/index/index /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries > /p/lava/parmik/partialMatcherProject/tools/bwa/SRR12432009/SarsGenome/SarsGenome_SRR12432009.sam