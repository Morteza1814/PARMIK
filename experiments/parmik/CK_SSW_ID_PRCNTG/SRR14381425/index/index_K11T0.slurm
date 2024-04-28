#!/bin/bash
#SBATCH --job-name="index_K11T0"
#SBATCH --error="index_K11T0.err"
#SBATCH --output="index_K11T0.out"
#SBATCH --mem=500G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado08 
#SBATCH --partition=main
#SBATCH --reservation=rgq5aw_79
/u/rgq5aw/GIT/PARMIK/src/parmik -a 0 -s 50 -c 150 -e 2 -t 0 -k 11 -i 27000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381425/SRR14381425_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381425/ -o /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/parmik/SRR12432009/NC_045512/ -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out -n /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/SRR12432009/index/kmerRanges/kmerRanges