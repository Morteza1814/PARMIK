#!/bin/bash
#SBATCH --job-name="parmikindex"
#SBATCH --error="index.err"
#SBATCH --output="index.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado03
#SBATCH --partition=main
/u/rgq5aw/GIT/PARMIK/src/parmik -a 0 -m 30 -s 50 -c 150 -e 2 -t 500 -k 15 -z 14 -i 17000000 -j 30000 -r /u/rgq5aw/GIT/PARMIK/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /u/rgq5aw/GIT/PARMIK/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/parmik/SRR12432009/index/ -o /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/parmik/SRR12432009/NC_045512/ -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out