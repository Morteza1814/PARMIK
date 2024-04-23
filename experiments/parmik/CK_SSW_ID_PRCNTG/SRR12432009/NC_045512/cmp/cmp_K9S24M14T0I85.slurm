#!/bin/bash
#SBATCH --job-name="cmp_K9S24M14T0I85"
#SBATCH --error="cmp_K9S24M14T0I85.err"
#SBATCH --output="cmp_K9S24M14T0I85.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=adriatic03 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 2 -l blast -s 24 -m 14 -c 150 -t 0 -k 9 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/NC_045512 -b /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws9_pi85_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848