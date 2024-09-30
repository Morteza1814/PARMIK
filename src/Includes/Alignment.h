#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#define DEBUG_MODE  0
#define EXE_MODE    1

uint64_t gAlignmentFoundWithNoPolish = 0;
uint64_t gAlignmentFoundWithPolish = 0;
uint64_t gAlignmentDumped = 0;
uint64_t gQueriesHaveAtLeastOneAlignment = 0;
uint64_t gAlignmentFoundWithPolishLargerThanBest = 0;

typedef struct Penalty
{
    int matchPenalty = 1;
    int mismatchPenalty = 1;
    int gapOpenPenalty = 1;
    int gapExtendPenalty = 1;
} Penalty;

typedef struct Alignment
{
    int readID = -1;            //read ID of the alignment
    int queryID = -1;            //query ID of the alignment
    string read;                //the region in read evaluated for the partial match
    string query;               //the region in query evaluated for the partial match
    string alignedRead;         //the region in read aligned for the partial match
    string alignedQuery;        //the region in query aligned for the partial match
    string editDistanceTypes;   //edit distance types of the region of the partial match
    string mismatchPositions;   // positions of mismatches in BWA
    uint8_t editDistance = 0;           //number of edit distances detected in the region
    uint8_t partialMatchSize = 0;       //partial match region size
    vector<uint16_t> editLocations;  //the positions of the edit distances in the partial match region
    string cigar;               //CIGAR string of the alignment
    uint8_t substitutions = 0;       //number of substitutions
    uint8_t inDels = 0;       //number of InDel
    uint8_t matches = 0;       //number of matched bp
    uint8_t readRegionStartPos = 0;
    uint8_t readRegionEndPos = 0;
    uint8_t queryRegionStartPos = 0;
    uint8_t queryRegionEndPos = 0;
    uint8_t flag = 0;                  // determines the strand for now
    uint8_t score = 0;
    uint8_t criteriaCode = 0;     // criteria code for the alignment
} Alignment;

typedef struct Alignmen_Baseline
{
    int readID = -1;            //read ID of the alignment
    int queryID = -1;            //query ID of the alignment
    string cigar;               //CIGAR string of the alignment
    uint8_t substitutions = 0;       //number of substitutions
    uint8_t inDels = 0;       //number of InDel
    uint8_t matches = 0;       //number of matched bp
    uint8_t flag = 0;                  // determines the strand for now
} Alignmen_Baseline;

#endif