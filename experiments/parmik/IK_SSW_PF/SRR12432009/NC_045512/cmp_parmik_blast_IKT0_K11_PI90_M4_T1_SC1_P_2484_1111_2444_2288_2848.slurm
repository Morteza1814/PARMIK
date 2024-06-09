#!/bin/bash
#SBATCH --job-name="cmp_parmik_blast_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848"
#SBATCH --error="cmp_parmik_blast_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.err"
#SBATCH --output="cmp_parmik_blast_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.out"
#SBATCH --mem=400G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado09 
#SBATCH --partition=cpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 2 -l blast -s 48 -c 150 -m 14 -t 0 -k 11 -d 90 -i 17000000 -j 30000 -x \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ \
-o /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512 \
-b /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/mts_50M/blastn_ws11_pi90_mts50M_results.out \
-p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848 \
-z /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848/parmikAln.txt