#!/bin/bash
#SBATCH --job-name="blastIndex"
#SBATCH --error="blastIndex.err"
#SBATCH --output="blastIndex.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado01
#SBATCH --partition=main
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/makeblastdb -in /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR12432009/SRR12432009_reads -dbtype nucl -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR12432009/index/blastIndex