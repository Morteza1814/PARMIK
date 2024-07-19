#!/bin/bash
#SBATCH --job-name="aln_ws11_pi90"
#SBATCH --error="aln_ws11_pi90.err"
#SBATCH --output="aln_ws11_pi90.out"
#SBATCH --mem=100G
#SBATCH --time=960
#SBATCH --partition=cpu

/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -word_size 11 -perc_identity 90 -max_target_seqs 50000000 \
-db /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381434/index/blastIndex \
-query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries \
-out /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381434/NC_045512/mts_50M/blastn_ws11_pi90_mts50M_results.out \
-outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps evalue"