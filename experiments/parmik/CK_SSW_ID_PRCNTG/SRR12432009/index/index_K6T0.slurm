#!/bin/bash
#SBATCH --job-name="index_K6T0"
#SBATCH --error="index_K6T0.err"
#SBATCH --output="index_K6T0.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado06 
#SBATCH --partition=main
/u/rgq5aw/GIT/PARMIK/src/parmik -a 0 -s 50 -c 150 -e 2 -t 0 -k 6 -i 17000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out -n /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/SRR12432009/index/kmerRanges/kmerRanges