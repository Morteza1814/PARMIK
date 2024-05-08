#!/bin/bash
#SBATCH --job-name="blcmp_parmik_K11S48M13T0I85"
#SBATCH --error="blcmp_parmik_K11S48M13T0I85.err"
#SBATCH --output="blcmp_parmik_K11S48M13T0I85.out"
#SBATCH --mem=500G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado09
#SBATCH --partition=main
#SBATCH --reservation=rgq5aw_79
/u/rgq5aw/GIT/PARMIK/src/parmik -a 4 -l parmik -s 48 -c 150 -m 13 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/NC_045512 -b /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/NC_045512/aln/noPolishing/pmAln_R48_PI85_L13_M3_E2_K11_T0_P_2484_1111_2444_2288_2848.txt -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848 -z /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR12432009/NC_045512