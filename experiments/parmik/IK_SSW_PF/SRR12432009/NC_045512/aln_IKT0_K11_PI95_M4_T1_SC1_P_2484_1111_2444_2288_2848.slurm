#!/bin/bash
#SBATCH --job-name="aln_IKT0_K11_PI95_M4_T1_SC1_P_2484_1111_2444_2288_2848"
#SBATCH --error="aln_IKT0_K11_PI95_M4_T1_SC1_P_2484_1111_2444_2288_2848.err"
#SBATCH --output="aln_IKT0_K11_PI95_M4_T1_SC1_P_2484_1111_2444_2288_2848.out"
#SBATCH --mem=400G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado04 
#SBATCH --partition=main
/u/rgq5aw/GIT/PARMIK/src/parmik -a 1 -s 48 -c 150 -m 14 -t 0 -k 11 -d 95 -i 17000000 -j 30000 -x \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ \
-o /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512 \
-b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out \
-p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848