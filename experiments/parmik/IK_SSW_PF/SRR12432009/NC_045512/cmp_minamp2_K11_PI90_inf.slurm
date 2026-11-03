#!/bin/bash
#SBATCH --job-name="cmp_minamp2_K11_PI90_inf"
#SBATCH --error="cmp_minamp2_K11_PI90_inf.err"
#SBATCH --output="cmp_minamp2_K11_PI90_inf.out"
#SBATCH --mem=500000
#SBATCH --time=4-00:00:00
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
stdbuf -oL -eL /u/rgq5aw/GIT/PARMIK/src/parmik -a 4 -l minimap2 -s 48 -c 150 -m 14 -t 0 -k 11 -d 90 -i 17000000 -j 30000 -x --comparison-threads 32 --comparison-summary-only \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ \
-o /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/NC_045512/minimap_temp \
-b /p/lava/parmik/partialMatcherProject/tools/minimap2/SRR12432009/NC_045512/minimap2_run2_ikt_inf.sam \
-p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848 \
-z /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR12432009/NC_045512
