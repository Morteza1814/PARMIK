#!/bin/bash
#SBATCH --job-name="aln_ws11_pi85_mts10000"
#SBATCH --error="aln_ws11_pi85_mts10000.err"
#SBATCH --output="aln_ws11_pi85_mts10000.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado07
#SBATCH --partition=main
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -db /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/index/blastIndex -query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -word_size 11 -perc_identity 85 -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws11_pi85_mts10000_results.out -max_target_seqs 10000 -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps"