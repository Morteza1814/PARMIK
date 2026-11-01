#!/bin/bash
#SBATCH --job-name="SRR12432009_Ebola"
#SBATCH --error="SRR12432009_Ebola.err"
#SBATCH --output="SRR12432009_Ebola.out"
#SBATCH --mem=80G
#SBATCH --time=1-00:00:00
#SBATCH --partition=gpu
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -db /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/index/blastIndex -query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/ebola_queries -word_size 11 -perc_identity 90 -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/Ebola/blastn_ws11_pi90_results.out -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps"
