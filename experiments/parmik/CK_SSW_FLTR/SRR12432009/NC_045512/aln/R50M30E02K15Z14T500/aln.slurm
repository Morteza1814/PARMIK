#!/bin/bash
#SBATCH --job-name="parmikaln"
#SBATCH --error="aln.err"
#SBATCH --output="aln.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=adriatic04 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 1 -x -m 30 -s 50 -c 150 -e 2 -t 500 -k 15 -z 14 -i 17000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/SRR12432009/index/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/alnResultsWithSW/SRR12432009/NC_045512/aln/R50M30E02K15Z14T500 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out