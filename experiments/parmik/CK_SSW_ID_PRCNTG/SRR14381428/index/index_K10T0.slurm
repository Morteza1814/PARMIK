#!/bin/bash
#SBATCH --job-name="index_K10T0"
#SBATCH --error="index_K10T0.err"
#SBATCH --output="index_K10T0.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cheetah03 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 0 -s 50 -c 150 -e 2 -t 0 -k 10 -i 40000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381428/SRR14381428_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381428/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out -n /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/SRR12432009/index/kmerRanges/kmerRanges