#!/bin/bash
#SBATCH --job-name="SarsGenome"
#SBATCH --error="SarsGenome.err"
#SBATCH --output="SarsGenome.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --partition=cpu
/p/lava/parmik/bwa/bwa mem /p/lava/parmik/partialMatcherProject/tools/bwa/SRR14381434/index/index /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries > /p/lava/parmik/partialMatcherProject/tools/bwa/SRR14381434/SarsGenome/SarsGenome.sam