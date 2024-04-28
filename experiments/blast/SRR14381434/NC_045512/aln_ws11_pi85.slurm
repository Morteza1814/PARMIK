#!/bin/bash
#SBATCH --job-name="aln_ws11_pi85"
#SBATCH --error="aln_ws11_pi85.err"
#SBATCH --output="aln_ws11_pi85.out"
#SBATCH --mem=240G
#SBATCH --time=2-00:00:00
#SBATCH --nodelist=cortado04
#SBATCH --partition=main
#SBATCH --reservation=rgq5aw_79
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -db /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381434/index/blastIndex -query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -word_size 11 -perc_identity 85 -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381434/NC_045512/blastn_ws11_pi85_results.out -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps"