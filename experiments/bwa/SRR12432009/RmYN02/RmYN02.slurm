#!/bin/bash
#SBATCH --job-name="RmYN02"
#SBATCH --error="RmYN02.err"
#SBATCH --output="RmYN02.out"
#SBATCH --mem=400G
#SBATCH --time=960
#SBATCH --nodelist=cortado05
#SBATCH --partition=main
/u/rgq5aw/GIT/GenomePartialMatching/bwa/bwa mem /bigtemp/rgq5aw/partialMatcherProject/datasets/metagenomicDB/SRR12432009/bwaIndex/SRR12432009_unique /bigtemp/rgq5aw/partialMatcherProject/datasets/queryFiles/covidGenomes/sampledFromCovid/query_150_RmYN02 > /bigtemp/rgq5aw/partialMatcherProject/tools/bwa/RmYN02_SRR12432009.sam