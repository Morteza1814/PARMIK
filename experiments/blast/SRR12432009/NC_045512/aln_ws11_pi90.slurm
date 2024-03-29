#!/bin/bash
#SBATCH --job-name="aln_ws11_pi90"
#SBATCH --error="aln_ws11_pi90.err"
#SBATCH --output="aln_ws11_pi90.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado06
#SBATCH --partition=main
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -db /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/index/blastIndex -query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -word_size 11 -perc_identity 90 -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws11_pi90_results.out -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps"