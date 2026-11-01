#!/bin/bash
#SBATCH --job-name="cmp_bowtie_K11_PI95"
#SBATCH --error="cmp_bowtie_K11_PI95.err"
#SBATCH --output="cmp_bowtie_K11_PI95.out"
#SBATCH --mem=500000
#SBATCH --time=4-00:00:00
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 4 -l bowtie -s 48 -c 150 -m 14 -t 0 -k 11 -d 95 -i 17000000 -j 30000 -x \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ \
-o /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR12432009/SarsGenome \
-b /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR12432009/SarsGenome/alignment_L11_a_local \
-p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848 \
-z /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR12432009/SarsGenome