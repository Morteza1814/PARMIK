#ifndef CONFIGS_H
#define CONFIGS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <iomanip>
// default value of config params
#define KMER_SZ 16
#define CONTIG_SZ 150
#define REGION_SZ 48
// #define MIN_EXACT_MATCH_LEN 32
#define NUMBER_OF_ALLOWED_EDIT_DISTANCES 2
#define NUMBER_OF_READS 1
#define NUMBER_OF_QUERIES 1

using namespace std;
// config params
typedef struct config
{
    string readDatabaseAddress;
    string queryFileAddress;
    string outputDir;
    string readFileName;
    string queryFileName;
    string otherToolOutputFileAddress;
    string offlineIndexAddress;
    string penaltyFileAddress;
    string kmerRangesFileAddress;
    string baselineBaseAddress;
    unsigned int readsCount;
    unsigned int queryCount;
    unsigned int kmerLength;
    unsigned int cheapKmerThreshold;
    unsigned int minExactMatchLen;
    unsigned int regionSize;
    unsigned int contigSize;
    unsigned int parmikMode;
    unsigned int numThreads;
    double identityPercentage = 0.9;
    double inDelPenalty = 2;
    double subPenalty = 1;
    // unsigned int numberOfKmers;
    unsigned int editDistance;
    bool isIndexOffline = false;
    bool isVerboseLog = false;
    bool noOutputFileDump = false;
    string otherTool;
    bool isSecondChanceOff = false;
} Config;

#endif