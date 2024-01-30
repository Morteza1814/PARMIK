#!/bin/bash
#SBATCH --job-name="NC_045512.2"
#SBATCH --error="NC_045512.2.err"
#SBATCH --output="NC_045512.2.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cheetah02
#SBATCH --partition=gpu
/p/lava/parmik/bwa/bwa mem /p/lava/parmik/partialMatcherProject/tools/bwa/SRR12432009/index/index /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries > /p/lava/parmik/partialMatcherProject/tools/bwa/SRR12432009/NC_045512/NC_045512.2_SRR12432009.sam