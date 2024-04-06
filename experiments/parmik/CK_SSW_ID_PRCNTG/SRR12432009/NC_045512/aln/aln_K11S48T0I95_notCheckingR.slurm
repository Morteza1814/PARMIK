#!/bin/bash
#SBATCH --job-name="aln_K11S48T0I95_notCheckingR"
#SBATCH --error="aln_K11S48T0I95_notCheckingR.err"
#SBATCH --output="aln_K11S48T0I95_notCheckingR.out"
#SBATCH --mem=240G
#SBATCH --time=960
#SBATCH --nodelist=titanx03 
#SBATCH --partition=gnolim
/u/rgq5aw/GIT/PARMIK/src/parmik -a 1 -s 48 -c 150 -t 0 -k 11 -d 95 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/NC_045512/notCheckingR -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848