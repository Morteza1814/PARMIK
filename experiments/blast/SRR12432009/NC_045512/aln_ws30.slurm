#!/bin/bash
#SBATCH --job-name="aln_ws30"
#SBATCH --error="aln_ws30.err"
#SBATCH --output="aln_ws30.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado01
#SBATCH --partition=main
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/blastn -db /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/index/blastIndex -query /p/lava/parmik/partialMatcherProject/datasets/queryFiles/covidQueryFiles/NC_045512_queries -word_size 30 -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/NC_045512/blastn_ws30_results.out -outfmt "6 qseqid sseqid qstart qend sstart send length mismatch sstrand qseq sseq"