#!/bin/bash
#SBATCH --job-name="aln_K11I85_NC_045512_Q_2000-2999"
#SBATCH --error="aln_K11I85_NC_045512_Q_2000-2999.err"
#SBATCH --output="aln_K11I85_NC_045512_Q_2000-2999.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --threads-per-core=2
#SBATCH --mem=60G
#SBATCH --time=4-00:00:00
#SBATCH --partition=nolim
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 27000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381425/SRR14381425_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_2000-2999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381425/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381425/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381425/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848