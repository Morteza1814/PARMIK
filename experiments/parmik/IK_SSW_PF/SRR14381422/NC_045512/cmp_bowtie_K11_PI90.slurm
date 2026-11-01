#!/bin/bash
#SBATCH --job-name="cmp_bowtie_K11_PI90"
#SBATCH --error="cmp_bowtie_K11_PI90.err"
#SBATCH --output="cmp_bowtie_K11_PI90.out"
#SBATCH --mem=500000
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 4 -l bowtie -s 48 -c 150 -m 14 -t 0 -k 11 -d 90 -i 27000000 -j 30000 -x \
-r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads \
-q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ \
-o /p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF/SRR14381422/NC_045512 \
-b /p/lava/parmik/partialMatcherProject/tools/bowtie2/SRR14381422/NC_045512/alignment_L11_a \
-p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848 \
-z /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/NC_045512