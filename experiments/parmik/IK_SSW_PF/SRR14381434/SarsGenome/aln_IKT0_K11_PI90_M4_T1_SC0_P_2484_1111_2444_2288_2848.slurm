#!/bin/bash
#SBATCH --job-name="aln_IKT0_K11_PI90_M4_T1_SC0_P_2484_1111_2444_2288_2848"
#SBATCH --error="aln_IKT0_K11_PI90_M4_T1_SC0_P_2484_1111_2444_2288_2848.err"
#SBATCH --output="aln_IKT0_K11_PI90_M4_T1_SC0_P_2484_1111_2444_2288_2848.out"
#SBATCH --mem=400G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado08 
#SBATCH --partition=cpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 1 -s 48 -c 150 -m 14 -t 0 -k 11 -d 90 -i 17000000 -j 30000 -u -x \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/SRR14381434_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381434/ \
-o /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR14381434/SarsGenome \
-b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381434/SarsGenome/blastn_ws30_results.out \
-p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848