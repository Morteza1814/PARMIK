#!/bin/bash
#SBATCH --job-name="cmp_K10S36M16T0I90"
#SBATCH --error="cmp_K10S36M16T0I90.err"
#SBATCH --output="cmp_K10S36M16T0I90.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=slurm1 
#SBATCH --partition=main
#SBATCH --reservation=rgq5aw_69
/u/rgq5aw/GIT/PARMIK/src/parmik -a 2 -l blast -s 36 -m 16 -c 150 -t 0 -k 10 -d 90 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/parmik/CK_SSW_ID_PRCNTG/SRR12432009/NC_045512 -b /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws10_pi90_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848