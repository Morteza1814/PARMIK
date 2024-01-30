#!/bin/bash
#SBATCH --job-name="parmikaln"
#SBATCH --error="aln.err"
#SBATCH --output="aln.out"
#SBATCH --mem=600G
#SBATCH --time=960
#SBATCH --nodelist=adriatic06
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 1 -x -m 30 -s 50 -c 150 -e 2 -t 0 -k 30 -z 29 -i 27000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/SRR14381434_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/gt/SRR14381434/index -o /p/lava/parmik/partialMatcherProject/tools/gt/SRR14381434/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out