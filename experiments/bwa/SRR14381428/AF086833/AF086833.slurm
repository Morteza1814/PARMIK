#!/bin/bash
#SBATCH --job-name="AF086833"
#SBATCH --error="AF086833.err"
#SBATCH --output="AF086833.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado03
#SBATCH --partition=main
/u/rgq5aw/GIT/GenomePartialMatching/bwa/bwa mem /bigtemp/rgq5aw/partialMatcherProject/tools/bwa/Index/SRR14381428/SRR14381428_unique /bigtemp/rgq5aw/partialMatcherProject/datasets/queryFiles/ebolaGenome/query_150_AF086833.fasta > /bigtemp/rgq5aw/partialMatcherProject/tools/bwa/SRR14381428/AF086833.sam