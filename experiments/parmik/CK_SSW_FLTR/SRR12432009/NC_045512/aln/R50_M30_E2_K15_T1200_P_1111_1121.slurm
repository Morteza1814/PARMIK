#!/bin/bash
#SBATCH --job-name="R50_M30_E2_K15_T1200_P_1111_1121"
#SBATCH --error="R50_M30_E2_K15_T1200_P_1111_1121.err"
#SBATCH --output="R50_M30_E2_K15_T1200_P_1111_1121.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=adriatic05
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 1 -x -m 30 -s 50 -c 150 -e 2 -t 1200 -k 15 -z 14 -i 17000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_FLTR/SRR12432009/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_FLTR/PenaltySets/1111_1121 