#!/bin/bash
#SBATCH --job-name="parmikFN"
#SBATCH --error="parmikFN.err"
#SBATCH --output="parmikFN.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado02 
#SBATCH --partition=main
/u/rgq5aw/GIT/PARMIK/src/parmik /p/lava/parmik/partialMatcherProject/tools/parmik/SRR12432009/RmYN02/cmp/blast/ParmikFnReads.txt /p/lava/parmik/partialMatcherProject/tools/parmik/SRR12432009/RmYN02/cmp/blast/kmerFrequLocFileAddress.txt  /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/RmYN02_queries /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads  30000 17000000
