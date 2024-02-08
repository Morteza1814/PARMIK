#!/bin/bash
#SBATCH --job-name="RmYN01"
#SBATCH --error="RmYN01.err"
#SBATCH --output="RmYN01.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado03
#SBATCH --partition=main
/u/rgq5aw/GIT/GenomePartialMatching/bwa/bwa mem /bigtemp/rgq5aw/partialMatcherProject/tools/bwa/Index/SRR14381428/SRR14381428_unique /bigtemp/rgq5aw/partialMatcherProject/datasets/queryFiles/covidGenomes/sampledFromCovid/query_150_RmYN01 > /bigtemp/rgq5aw/partialMatcherProject/tools/bwa/SRR14381428/RmYN01.sam