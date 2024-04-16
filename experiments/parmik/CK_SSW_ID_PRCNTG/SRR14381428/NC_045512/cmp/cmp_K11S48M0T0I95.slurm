#!/bin/bash
#SBATCH --job-name="cmp_K11S48M0T0I95"
#SBATCH --error="cmp_K11S48M0T0I95.err"
#SBATCH --output="cmp_K11S48M0T0I95.out"
#SBATCH --mem=240G
#SBATCH --time=960
#SBATCH --nodelist=titanx02 
#SBATCH --partition=gnolim
/u/rgq5aw/GIT/PARMIK/src/parmik -a 2 -l blast -s 48 -c 150 -t 0 -k 11 -d 95 -i 40000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381428/SRR14381428_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381428/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR14381428/NC_045512 -b /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381428/NC_045512/blastn_ws11_pi95_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848