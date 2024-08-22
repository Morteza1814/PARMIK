#!/bin/bash
#SBATCH --job-name="aln_K11I85_NC_045512_SRR14381418_Q_15000-20000"
#SBATCH --error="aln_K11I85_NC_045512_SRR14381418_Q_15000-20000.err"
#SBATCH --output="aln_K11I85_NC_045512_SRR14381418_Q_15000-20000.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --threads-per-core=2
#SBATCH --mem=60G
#SBATCH --time=4-00:00:00
#SBATCH --partition=gpu
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381418/SRR14381418_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_15000-15999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381418/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381418/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381418/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381418/SRR14381418_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_16000-16999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381418/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381418/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381418/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381418/SRR14381418_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_17000-17999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381418/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381418/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381418/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381418/SRR14381418_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_18000-18999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381418/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381418/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381418/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
echo "Starting the job at $(date)"
/u/rgq5aw/GIT/PARMIK/src/parmik -a 3 -s 48 -c 150 -t 0 -k 11 -d 85 -i 17000000 -j 30000 -x -r /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381418/SRR14381418_reads -q /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_chunks/NC_045512_queries_19000-19999 -f /p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR14381418/ -o /p/lava/parmik/partialMatcherProject/tools/baseline_ssw/SRR14381418/NC_045512 -b /u/rgq5aw/GIT/PARMIK/partialMatcherProject/tools/blast/SRR14381418/NC_045512/blastn_ws30_results.out -p /u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_ID_PRCNTG/PenaltySets/2484_1111_2444_2288_2848
