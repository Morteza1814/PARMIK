# PARMIK
**PA**rtial **R**ead **M**atching with **I**nexpensive **K**-mers

### PARMIK Arguments:

- `-a`, `--mode`: PARMIK mode
  - **Values**: 
    - `PARMIK_MODE_INDEX (0)`
    - `PARMIK_MODE_ALIGN (1)`
    - `PARMIK_MODE_COMPARE (2)`
    - `PARMIK_MODE_BASELINE (3)`
    - `PARMIK_MODE_CMP_BASELINE (4)`
- `-b`, `--toolFileAddress`: Other Tool Alignment File Address
  - The address of the output of the other tool (BLAST, BWA, etc)
- `-c`, `--contigSize`: Contig Size
  - Length of the contigs (default = 150)
- `-d`, `--percentageIdentity`: Percentage Identity
- `-e`, `--editDistance`: Max Edit Distance (i/d/s)
- `-f`, `--ikiAddress`: Inexpensive K-mer Index Address
- `-h`, `--help`: Help
- `-i`, `--readCount`: Number of Metageomic Reads
- `-j`, `--queryCount`: Number of Queries
- `-k`, `--kmerLen`: K-mer Length
- `-l`, `--otherTool`: The Other Tool (bwa, blast, etc)
- `-m`, `--minExactMatchLen`: Minimum Exact Match Length
- `-n`, `--kmerRangesFileAddress`: K-mer Ranges File Address
- `-o`, `--outputDir`: Output Directory
- `-p`, `--penaltyFileAddress`: Penalty File Address
- `-q`, `--query`: Query File Address
- `-r`, `--read`: Metageomic Read Data Base Address
- `-s`, `--regionSize`: Region Size
- `-t`, `--cheapKmerThreshold`: Cheap Kmer Threshold
- `-u`, `--isSecondChanceOff`: Turn Second Chance Off
- `-v`, `--verboseLog`: Verbose Logging
- `-w`, `--numThreads`: Number of Threads
- `-x`, `--isIndexOffline`: Is the read index offline
- `-z`, `--baselineBaseAddress`: BaseLine file base address

