#!/bin/bash
#SBATCH --job-name="aln_K11I85_SarsGenome_Q_0-999"
#SBATCH --error="aln_K11I85_SarsGenome_Q_0-999.err"
#SBATCH --output="aln_K11I85_SarsGenome_Q_0-999.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --threads-per-core=2
#SBATCH --mem=400G
#SBATCH --time=4-00:00:00
#SBATCH --nodelist=cheetah03 
#SBATCH --partition=gpu
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381434/SRR14381434_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_chunks/SarsGenome_queries_0-999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381434/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381434/SarsGenome -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381434/SarsGenome/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848