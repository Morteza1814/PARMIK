#!/bin/bash
#SBATCH --job-name="aln_ws9_pi85"
#SBATCH --error="aln_ws9_pi85.err"
#SBATCH --output="aln_ws9_pi85.out"
#SBATCH --mem=400G
#SBATCH --time=2-00:00:00
#SBATCH --nodelist=cheetah02
#SBATCH --partition=gpu
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -db /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/index/blastIndex -query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -word_size 9 -perc_identity 85 -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws9_pi85_results.out -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps"