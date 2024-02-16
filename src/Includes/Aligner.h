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
#include "sw/ssw_cpp.h"

using namespace std;

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
      uint16_t substitutions = 0;       //number of substitutions
      uint16_t inDels = 0;       //number of InDel
      uint16_t matches = 0;       //number of matched bp
      uint16_t readRegionStartPos = 0;
      uint16_t readRegionEndPos = 0;
      uint16_t queryRegionStartPos = 0;
      uint16_t queryRegionEndPos = 0;
      uint16_t flag = 0;                  // determines the strand for now
      uint16_t score = 0;
  } Alignment;

class Aligner {
public:

  void parseCigar(const std::string& cigar, uint16_t& matches, uint16_t& substitutions, uint16_t& inDels, std::vector<unsigned int>& editLocations) {
    substitutions = inDels = 0;
    int currentPos = 0; // Current position in the read
    int len = cigar.length();

    for (int i = 0; i < len; ++i) {
        // Parse the numeric part of CIGAR operation
        int num = 0;
        while (i < len && isdigit(cigar[i])) {
            num = num * 10 + (cigar[i] - '0');
            ++i;
        }

        // Extract the CIGAR operation
        char op = cigar[i];

        // Perform actions based on CIGAR operation
        switch (op) {
            case '=':
                // Match or mismatch (substitution)
                currentPos += num;
                matches += num;
                break;
            case 'X':
                // Substitution
                substitutions += num;
                for (int j = 0; j < num; ++j) {
                    editLocations.push_back(currentPos + j);
                }
                currentPos += num;
                break;
            case 'I':
            case 'D':
                // Insertion
                inDels += num;
                for (int j = 0; j < num; ++j) {
                    editLocations.push_back(currentPos + j);
                }
                currentPos += num;
                break;
            case 'S':
                // Soft clipping
                currentPos += num;
                break;
            default:
                // Unsupported CIGAR operation
                std::cerr << "Unsupported CIGAR operation: " << op << std::endl;
                break;
        }
    }
}

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
    parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, aln.editPositions);
  }

};

#endif