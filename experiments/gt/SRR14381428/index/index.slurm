#!/bin/bash
#SBATCH --job-name="parmikindex"
#SBATCH --error="index.err"
#SBATCH --output="index.out"
#SBATCH --mem=800G
#SBATCH --time=960
#SBATCH --nodelist=adriatic01 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 0 -m 30 -s 50 -c 150 -e 2 -t 0 -k 30 -z 29 -i 40000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381428/SRR14381428_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/gt/SRR14381428/index/ -o /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/parmik/SRR12432009/NC_045512/ -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out