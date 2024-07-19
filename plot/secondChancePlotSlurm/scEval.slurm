#!/bin/bash
#SBATCH --job-name="scEval"
#SBATCH --error="scEval.err"
#SBATCH --output="scEval.out"
#SBATCH --mem=400G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado02 
#SBATCH --partition=cpu
module load python3
python3 ../polishingEval.py /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT0_K11_PI90_M4_T1_SC0_P_2484_1111_2444_2288_2848/SC_vs_NoSC.txt /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/IKT0_K11_PI90_M4_T1_SC0_P_2484_1111_2444_2288_2848/SC_vs_NoSC.png