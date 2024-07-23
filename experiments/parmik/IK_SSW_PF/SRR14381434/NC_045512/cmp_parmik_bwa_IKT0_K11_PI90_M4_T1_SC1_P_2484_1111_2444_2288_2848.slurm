#!/bin/bash
#SBATCH --job-name="cmp_parmik_bwa_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848"
#SBATCH --error="cmp_parmik_bwa_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.err"
#SBATCH --output="cmp_parmik_bwa_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.out"
#SBATCH --mem=450G
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 2 -l bwa -s 48 -c 150 -m 14 -t 0 -k 11 -d 90 -i 17000000 -j 30000 -x \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/SRR14381434_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381434/ \
-o /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR14381434/NC_045512 \
-b /p/lava/parmik/partialMatcherProject/tools/bwa/SRR14381434/NC_045512/NC_045512.2_SRR14381434.sam \
-p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848 \
-z /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR14381434/NC_045512/IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848/parmikAln.txt