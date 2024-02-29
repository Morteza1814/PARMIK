#!/bin/bash
#SBATCH --job-name="parmikfn"
#SBATCH --error="fn.err"
#SBATCH --output="fn.out"
#SBATCH --mem=250G
#SBATCH --time=960
#SBATCH --nodelist=puma02 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik /p/lava/parmik/partialMatcherProject/tools/parmik/alnResultsWithSW/SRR12432009/NC_045512/cmp/blast/ParmikFnReads.txt /p/lava/parmik/partialMatcherProject/tools/parmik/alnResultsWithSW/SRR12432009/NC_045512/cmp/blast/kmerFrequLocFileAddress.txt /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads 30000 30000000
