#!/bin/bash
#SBATCH --job-name="fneval_IKT2000_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848"
#SBATCH --error="fneval_IKT2000_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.err"
#SBATCH --output="fneval_IKT2000_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.out"
#SBATCH --mem=400G
#SBATCH --time=08:00:00
#SBATCH --nodelist=cortado02 
#SBATCH --partition=main
/u/rgq5aw/GIT/PARMIK/src/parmik \
/p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT2000_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848/parmikAln.txt \
30000 \
17000000 \
/p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT2000_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848/parmikExpensiveKmerFNs.txt \
/p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ek_T2000_K11_r17000000 \
/p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads \
11 \
/p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT2000_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848