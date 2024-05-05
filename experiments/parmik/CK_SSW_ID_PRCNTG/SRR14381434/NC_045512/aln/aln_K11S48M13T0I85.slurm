#!/bin/bash
#SBATCH --job-name="aln_K11S48M13T0I85"
#SBATCH --error="aln_K11S48M13T0I85.err"
#SBATCH --output="aln_K11S48M13T0I85.out"
#SBATCH --mem=500G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado03 
#SBATCH --partition=main
/u/rgq5aw/GIT/PARMIK/src/parmik -a 1 -s 48 -c 150 -m 13 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/SRR14381434_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381434/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR14381434/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381434/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848