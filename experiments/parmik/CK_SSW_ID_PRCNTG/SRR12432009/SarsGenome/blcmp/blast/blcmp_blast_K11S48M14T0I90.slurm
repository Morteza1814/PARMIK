#!/bin/bash
#SBATCH --job-name="blcmp_blast_K11S48M14T0I90"
#SBATCH --error="blcmp_blast_K11S48M14T0I90.err"
#SBATCH --output="blcmp_blast_K11S48M14T0I90.out"
#SBATCH --mem=500G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cheetah02 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 4 -l blast -s 48 -c 150 -m 14 -t 0 -k 11 -d 90 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/SarsGenome -b /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/SarsGenome/blastn_ws11_pi90_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848 -z /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR12432009/SarsGenome