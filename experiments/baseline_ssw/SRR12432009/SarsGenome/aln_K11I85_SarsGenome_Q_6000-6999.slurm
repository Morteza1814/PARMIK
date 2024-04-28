#!/bin/bash
#SBATCH --job-name="aln_K11I85_SarsGenome_Q_6000-6999"
#SBATCH --error="aln_K11I85_SarsGenome_Q_6000-6999.err"
#SBATCH --output="aln_K11I85_SarsGenome_Q_6000-6999.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --threads-per-core=2
#SBATCH --mem=400G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cortado05 
#SBATCH --partition=main
#SBATCH --reservation=rgq5aw_79
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_chunks/SarsGenome_queries_6000-6999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR12432009/SarsGenome -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR12432009/SarsGenome/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848