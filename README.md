# PARMIK
**PA**rtial **R**ead **M**atching with **I**nexpensive **K**-mers

## Directory Structure ##
- `src`: Contains source code for the project
-  `scripts`: Contains scripts
## Prerequisites ##
Before you begin, ensure you have the following installed on your system, (section):
- Ubuntu: All testing has been done on Ubuntu 22.04+ Operating System.
- GCC: The GNU Compiler Collection, specifically g++9 which supports C++11 or later.
- Make: The build utility to automate the compilation.
- OpenMP: Support for parallel programming in C++.
- Python3 for running the scripts

## How to Compile ##
To compile, use:
```bash
make
```
To clean up all compiled files:
```bash
make clean
```
## Download Datasets ##
To download datasets, we used _SRA Toolkit (v3.0.7)_. Here is the command we used to download a metagenomic dataset (SRR12432009):
```bash
sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR12432009 -p --fasta --outdir <outputDir>
 ```
Replace `<outputDir>` with the path to your desired output directory.
## PARMIK parameters:
Below are the PARMIK's parameters in alphabetical order:
- `-a`, `--mode`: PARMIK mode (*required*)
  - PARMIK operation mode. It can get these values:
    - `PARMIK_MODE_INDEX (0)`
    - `PARMIK_MODE_ALIGN (1)`
    - `PARMIK_MODE_COMPARE (2)`
    - `PARMIK_MODE_BASELINE (3)`
    - `PARMIK_MODE_CMP_BASELINE (4)`
- `-b`, `--toolFileAddress`: Other Tool Alignment File Address (*required for compare mode*)
  - The address of the output of the other tool (BLAST, BWA, etc)
- `-c`, `--contigSize`: Contig Size (*default = 150*)
  - Length of the contigs 
- `-d`, `--percentageIdentity`: Percentage Identity (*default = 90%*)
  - Minimum Percentage of Identity in the alignment
- `-e`, `--editDistance`: Max Edit Distance (i/d/s) (*default = 2*)
  - Maximum edit distance (including Substitutions and InDels) allowed in the alignment
- `-f`, `--ikiAddress`: Inexpensive K-mer Index Address (*required*)
  - The path to the Inexpensive K-mer Index (IKI)
- `-h`, `--help`: Help
- `-i`, `--readCount`: Number of Metagenomic Reads (*default = 1*)
  - Number of reads in the Metagenomic dataset 
- `-j`, `--queryCount`: Number of Queries (*default = 1*)
  - Number of queries in the Query dataset
- `-k`, `--kmerLen`: K-mer Length (*default = 16*)
  - Length of the K-mer
- `-l`, `--otherTool`: The Other Tool Name (*required for compare mode*)
  - Name of other tool (bwa, blast, etc)
- `-m`, `--minExactMatchLen`: Minimum Exact Match Length (*default = 0*)
  - Min length of exact match required for alignment
  - M = (minExactMatchLen - K + 1)
- `-n`, `--kmerRangesFileAddress`: K-mer Ranges File Address
  - K-mer ranges file address required for calculating the inexpensive k-mer threshold
- `-o`, `--outputDir`: Output Directory (*required for all modes except index mode*)
  - Directory to dump the alignment results
- `-p`, `--penaltyFileAddress`: Penalty File Address
  - The penalty score sets used for the alignment step
- `-q`, `--query`: Query File Address (*required*)
  - Path to the query dataset file
- `-r`, `--read`: Metageomic Read Data Base Address (*required*)
  - Path to the read metagenomic dataset file
- `-s`, `--regionSize`: Region Size (*default = 48*)
  - Minimum size of the alignment
- `-t`, `--cheapKmerThreshold`: Cheap (Inexpensive) Kmer Threshold (*required*)
  - -t 0: includes all k-mers in the IKI (Inexpensive K-mer Index).
- `-u`, `--isSecondChanceOff`: Turn Second Chance Off
  - Turn off the second chance
  - This is a flag. When included, it disables the second chance mechanism (sets the flag to true).
- `-v`, `--verboseLog`: Verbose Logging (*default = false*)
- `-w`, `--numThreads`: Number of Threads (*default = 1*)
- `-x`, `--isIndexOffline`: Is the read index offline
  - This is a flag. When included, it enables the write/read the IKI to/from storage.
  - If not included, PARMIK creates and use the IKI on the fly
- `-z`, `--baselineBaseAddress`: BaseLine file base address  (*required for compare mode*)
  - Base address of the baseline alignment outputs

### Help ###
To display help for general usage:
```bash
./parmik --help
```

