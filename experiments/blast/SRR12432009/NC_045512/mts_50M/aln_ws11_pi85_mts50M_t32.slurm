#!/bin/bash
#SBATCH --job-name="aln_ws11_pi85_mts50M_t32"
#SBATCH --error="aln_ws11_pi85_mts50M_t32.err"
#SBATCH --output="aln_ws11_pi85_mts50M_t32.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=2
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=adriatic05
#SBATCH --partition=gpu
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -word_size 11 -perc_identity 85 -max_target_seqs 50000000 -num_threads 32 \
-db /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/index/blastIndex \
-query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-out /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/mts_50M/blastn_ws11_pi85_mts50M_t32_results.out \
-outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps evalue"