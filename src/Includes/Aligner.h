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
    vector<unsigned int> editLocations;  //the positions of the edit distances in the partial match region
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
private:
    size_t regionSize;
    uint32_t allowedEditDistance;
    uint32_t contigSize;
public:
    Aligner(size_t R, uint32_t a, uint32_t c) : regionSize(R), allowedEditDistance(a), contigSize(c){}

    void findMemAndExtend(Alignment & aln){
        uint16_t memSize = aln.editLocations[0];
        int start = -1, end = 0;
        int memStart = -1, memEnd = -1;
        //find MEM
        uint16_t memSize = 0, maxMemSize = 0;
        for (int i = 0; i < aln.editLocations.size(); ++i) {
            if (i == 0) {
                end = i;
                memSize = aln.editLocations[i];
            } else {
                start = i - 1;
                end = i;
                memSize = aln.editLocations[i] - aln.editLocations[i-1] - 1;
            }
            if (memSize > maxMemSize) {
                maxMemSize = memSize;
                memStart = start;
                memEnd = end;
            }
        }
        //(aln.queryRegionEndPos - aln.queryRegionStartPos) to make the beginning of the region 0
        if (aln.queryRegionEndPos - aln.queryRegionStartPos - aln.editLocations[aln.editLocations.size()-1] > maxMemSize) {
            maxMemSize = aln.queryRegionEndPos - aln.queryRegionStartPos - aln.editLocations[aln.editLocations.size()-1];
            memStart = aln.editLocations.size() - 1;
            memEnd = -1;
        }
        //extend to the left as far as possible
        start = memStart;
        end = memEnd;
        uint16_t ed = 0;
        while (start >= 0 && ed <= allowedEditDistance)
        {
            start--;
            ed++;
        }
        if (ed >= allowedEditDistance) 
        {
            ed = allowedEditDistance;
            if (start < 0) 
            {
                start = -1;
            }
        } else 
        {
            if (start < 0)
            {
                start = -1;
            }
        }
        //extend to the right as far as possible
        if (ed < allowedEditDistance)
        {
            while (end < aln.editLocations.size() && ed < allowedEditDistance)
            {
                end++;
                ed++;
            }
            if (ed >= allowedEditDistance)
            {
                ed = allowedEditDistance;
                if (end >= aln.editLocations.size())
                {
                    end = -1;
                }
            }
            else
            {
                end = -1;
            }
        }
        // calculate the largest aligment
        uint16_t maxAlnSize = 0, uint16_t alnSize = 0;
        if (end == -1 && start == -1) {
            alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos;
        } else if (end == -1 && start > -1) {
            alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos - aln.editLocations[start];
        } else if (end > -1 && start == -1) {
            alnSize = aln.editLocations[end];
        } else {
            alnSize = aln.editLocations[end] - aln.editLocations[start] - 1;
        }
        maxAlnSize = alnSize;
        int maxAlnStart = start, maxAlnEnd = end;
        start++;
        end++;
        while (start < memStart && end <= aln.editLocations.size())
        {
            if (end >= aln.editLocations.size())
            {
                end = -1;
            }
            if (end == -1 && start == -1) {
            alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos;
            } else if (end == -1 && start > -1) {
                alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos - aln.editLocations[start];
            } else if (end > -1 && start == -1) {
                alnSize = aln.editLocations[end];
            } else {
                alnSize = aln.editLocations[end] - aln.editLocations[start] - 1;
            }
            if (alnSize > maxAlnSize)
            {
                maxAlnSize = alnSize;
                maxAlnStart = start;
                maxAlnEnd = end;
            }
            start++;
            end++;
        }
        //find the actual location of max  alignment start and end
        uint16_t maxAlnStartPos = 0, maxAlnEndPos = 0;
        if (maxAlnStart < 0){
            maxAlnStartPos = 0;
        } else {
            maxAlnStartPos = aln.editLocations[maxAlnStart] + 1;
        }
        if (maxAlnEnd == -1) {
            maxAlnEndPos = aln.queryRegionEndPos - aln.queryRegionStartPos;
        } else {
            maxAlnEndPos = aln.editLocations[maxAlnEnd] - 1;
        }
        // set the region start and end in R and Q
        aln.readRegionStartPos = maxAlnStartPos + aln.readRegionStartPos;
        aln.readRegionEndPos = maxAlnEndPos + aln.readRegionStartPos;
        aln.queryRegionStartPos = maxAlnStartPos + aln.queryRegionStartPos;
        aln.queryRegionEndPos = maxAlnEndPos + aln.queryRegionStartPos;
        //change the cigar based on the new region
        std::string newCigar = "";            
    }  

    Alignment align(uint16_t matchPen, uint16_t subPen, uint16_t gapoPen, uint16_t gapextPen)
    {
        Alignment aln;
        //prepare the aln
        smithWatermanAligner(aln, matchPen, subPen, gapoPen, gapextPen);
        if (aln.editDistance <= allowedEditDistance) {
            //if (checkAlnBasedOnCriteria(aln)) 
            // return aln;
            //else return NuLL
        } else {
            //find the mem
            findMemAndExtend(aln);

        }

    }

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
                    // currentPos += num;
                    break;
                default:
                    // Unsupported CIGAR operation
                    std::cerr << "Unsupported CIGAR operation: " << op << std::endl;
                    break;
            }
        }
    }

    void smithWatermanAligner(Alignment &aln, uint16_t matchPen, uint16_t subPen, uint16_t gapoPen, uint16_t gapextPen)
    {
        int32_t maskLen = strlen(aln.query.c_str())/2;
        maskLen = maskLen < 15 ? 15 : maskLen;

        // Declares a default Aligner
        StripedSmithWaterman::Aligner aligner(matchPen, subPen, gapoPen, gapextPen);
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
        parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, aln.editLocations);
    }

};

#endif