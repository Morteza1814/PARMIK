#ifndef ALIGNER_H
#define ALIGNER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include "Utils.h"

using namespace std;

class Aligner {
public:
  typedef struct Alignment
  {
      int readID = -1;            //read ID of the alignment
      int queryID = -1;            //query ID of the alignment
      string read;                //the region in read evaluated for the partial match
      string query;               //the region in query evaluated for the partial match
      string alignedRead;         //the region in read aligned for the partial match
      string alignedQuery;        //the region in query aligned for the partial match
      string editDistanceTypes;   //edit distance types of the region of the partial match
      unsigned int editDistance = 0;           //number of edit distances detected in the region
      int partialMatchSize = 0;       //partial match region size
      vector<unsigned int> editPositions;  //the positions of the edit distances in the partial match region
      string cigar;               //CIGAR string of the alignment
      uint32_t numberOfSub = 0;       //number of substitutions
      uint32_t numberOfInDel = 0;       //number of InDel
      uint32_t numberOfMatches = 0;       //number of matched bp
      uint16_t readRegionStartPos = 0;
      uint16_t readRegionEndPos = 0;
      uint16_t queryRegionStartPos = 0;
      uint16_t queryRegionEndPos = 0;
      uint16_t flag = 0;                  // determines the strand for now
      uint16_t score = 0;
  } Alignment;

  void smithWatermanAligner(Alignment &aln)
  {
    int32_t maskLen = strlen(aln.query.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner(1,1,1,1);
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the query to the ref
    aligner.Align(aln.query.c_str(), aln.read.c_str(), aln.read.size(), filter, &alignment, maskLen); 

    aln.readRegionStartPos = alignment.ref_begin;
    aln.readRegionEndPos = alignment.ref_end;
    aln.queryRegionStartPos = alignment.query_begin;
    aln.queryRegionEndPos = alignment.query_end;
    aln.cigar = alignment.cigar_string;
    aln.editDistance = alignment.mismatches;
    aln.score = alignment.sw_score;
  }

};

#endif