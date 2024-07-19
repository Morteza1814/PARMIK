#!/bin/bash
#SBATCH --job-name="index_K10T0"
#SBATCH --error="index_K10T0.err"
#SBATCH --output="index_K10T0.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado06 
#SBATCH --partition=cpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 0 -s 50 -c 150 -t 0 -k 10 -i 17000000 -j 30000 \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/SRR14381434_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381434/ \
-n /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/SRR12432009/index/kmerRanges/kmerRanges