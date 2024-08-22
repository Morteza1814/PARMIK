#!/bin/bash
#SBATCH --job-name="aln_K11I85_SarsGenome_Q_0-5000"
#SBATCH --error="aln_K11I85_SarsGenome_Q_0-5000.err"
#SBATCH --output="aln_K11I85_SarsGenome_Q_0-5000.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --threads-per-core=2
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu
# echo "Starting the job at $(date)"
# /u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_chunks/SarsGenome_queries_0-999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/SarsGenome -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/SarsGenome/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_chunks/SarsGenome_queries_1000-1999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/SarsGenome -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/SarsGenome/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_chunks/SarsGenome_queries_2000-2999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/SarsGenome -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/SarsGenome/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_chunks/SarsGenome_queries_3000-3999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/SarsGenome -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/SarsGenome/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/SarsGenome_chunks/SarsGenome_queries_4000-4999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/SarsGenome -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/SarsGenome/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848