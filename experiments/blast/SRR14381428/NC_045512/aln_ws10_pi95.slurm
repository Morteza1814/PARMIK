#!/bin/bash
#SBATCH --job-name="aln_ws10_pi95"
#SBATCH --error="aln_ws10_pi95.err"
#SBATCH --output="aln_ws10_pi95.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=adriatic04
#SBATCH --partition=gpu
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -db /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381428/index/blastIndex -query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -word_size 10 -perc_identity 95 -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381428/NC_045512/blastn_ws10_pi95_results.out -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq pident gaps"