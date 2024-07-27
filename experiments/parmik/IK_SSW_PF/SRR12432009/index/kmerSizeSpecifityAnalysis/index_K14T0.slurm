#!/bin/bash
#SBATCH --job-name="index_K14T0"
#SBATCH --error="index_K14T0.err"
#SBATCH --output="index_K14T0.out"
#SBATCH --mem=400000
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu
#SBATCH --nodelist=cortado10
#SBATCH --reservation=rgq5aw_112
/u/rgq5aw/GIT/PARMIK/src/parmik -a 0 -s 50 -c 150 -t 0 -k 14 -i 17000000 -j 30000 \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/kmerSizeSpecifityAnalysis/ \
-n /u/rgq5aw/GIT/PARMIK/experiments/parmik/IK_SSW_PF/SRR12432009/index/kmerRanges/kmerRanges