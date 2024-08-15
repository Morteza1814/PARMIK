#!/bin/bash
#SBATCH --job-name="aln_K11I85_NC_045512_Q_25000-30000"
#SBATCH --error="aln_K11I85_NC_045512_Q_25000-30000.err"
#SBATCH --output="aln_K11I85_NC_045512_Q_25000-30000.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --threads-per-core=2
#SBATCH --mem=60G
#SBATCH --time=4-00:00:00
#SBATCH --partition=cpu
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_25000-25999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_26000-26999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_27000-27999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_28000-28999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381422/SRR14381422_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_29000-29999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381422/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381422/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381422/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
