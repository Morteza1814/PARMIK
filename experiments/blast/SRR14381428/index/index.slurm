#!/bin/bash
#SBATCH --job-name="blastIndex"
#SBATCH --error="blastIndex.err"
#SBATCH --output="blastIndex.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado03
#SBATCH --partition=main
/p/lava/parmik/blast/ncbi-blast-2.15.0+/bin/makeblastdb -in /p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/SRR14381428/SRR14381428_reads -dbtype nucl -out /p/lava/parmik/partialMatcherProject/tools/blast/SRR14381428/index/blastIndex