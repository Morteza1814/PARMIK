#!/bin/bash
#SBATCH --job-name="scEval_IKT0_K11_PI90_M4_T1_P_2484_1111_2444_2288_2848"
#SBATCH --error="scEval_IKT0_K11_PI90_M4_T1_P_2484_1111_2444_2288_2848.err"
#SBATCH --output="scEval_IKT0_K11_PI90_M4_T1_P_2484_1111_2444_2288_2848.out"
#SBATCH --mem=600G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=adriatic05 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik \
/p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848/parmikAln.txt \
/p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT0_K11_PI90_M4_T1_SC0_P_2484_1111_2444_2288_2848/parmikAln.txt \
30000 \
/p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT0_K11_PI90_M4_T1_SC0_P_2484_1111_2444_2288_2848 \
/p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries  