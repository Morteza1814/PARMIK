#!/bin/bash
#SBATCH --job-name="parmikcmp"
#SBATCH --error="cmp.err"
#SBATCH --output="cmp.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado02 
#SBATCH --partition=main
#SBATCH --reservation=rgq5aw_33
/u/rgq5aw/GIT/PARMIK/src/parmik -a 2 -l blast -x -m 30 -s 50 -c 150 -e 2 -t 500 -k 15 -z 14 -i 27000000 -j 30000 -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381425/SRR14381425_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/SRR14381425/index/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/SRR14381425/NC_045512 -b /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381425/NC_045512/blastn_ws30_results.out
