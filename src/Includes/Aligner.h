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
    size_t editDistance = 0;           //number of edit distances detected in the region
    size_t partialMatchSize = 0;       //partial match region size
    vector<uint16_t> editLocations;  //the positions of the edit distances in the partial match region
    string cigar;               //CIGAR string of the alignment
    size_t substitutions = 0;       //number of substitutions
    size_t inDels = 0;       //number of InDel
    size_t matches = 0;       //number of matched bp
    size_t readRegionStartPos = 0;
    size_t readRegionEndPos = 0;
    size_t queryRegionStartPos = 0;
    size_t queryRegionEndPos = 0;
    size_t flag = 0;                  // determines the strand for now
    size_t score = 0;
    size_t criteriaCode = 0;     // criteria code for the alignment
} Alignment;

template <typename contigIndT>
class Aligner {
private:
    size_t regionSize;
    size_t allowedEditDistance;
    size_t contigSize;
    size_t minExactMatchLength;
public:
    Aligner(size_t R, size_t a, size_t c, size_t m) : regionSize(R), allowedEditDistance(a), contigSize(c), minExactMatchLength(m) {}

    string convertCigarToStr(const string& cigar) {
        stringstream simplifiedCigar;
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
                    simplifiedCigar << string(num, '=');
                    break;
                case 'X':
                    // Match or mismatch (substitution)
                    simplifiedCigar << string(num, 'X');
                    break;
                case 'I':
                    // Insertion
                    simplifiedCigar << string(num, 'I');
                    break;
                case 'D':
                    // Deletion
                    simplifiedCigar << string(num, 'D');
                    break;
                case 'S':
                    // Soft clipping
                    simplifiedCigar << string(num, 'S');
                    break;
                default:
                    // Handle other CIGAR operations if needed
                    break;
            }
        }

        return simplifiedCigar.str();
    }

    string convertStrToCigar(const string cigarStr) {
        stringstream cigar;
        char prev = cigarStr[0];
        uint16_t num = 0;
        for (size_t i = 1; i < cigarStr.length(); ++i) {
            char cur = cigarStr[i];
            num++;
            if (prev != cur)
            {
                cigar << num << prev;
                num = 0;
                prev = cur;
            }
        }
        num++;
        cigar << num << prev;
        return cigar.str();
    }

    void findMemAndExtend(Alignment &aln){
        int start = -1, end = -1;
        int memStart = -1, memEnd = -1;
        string cigarStr = convertCigarToStr(aln.cigar);
        // cout << "CIGAR str: " << cigarStr << endl;
        // insertions adds to the read positions and deletions adds to the query positions
        size_t alignmentOffset = (cigarStr.size() > contigSize) ? cigarStr.size() - contigSize : 0; // if the alignment size goes over 150bp
        //find MEM
        uint16_t memSize = 0, maxMemSize = 0;
        for (size_t i = 0; i < aln.editLocations.size(); ++i) {
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
        if (aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[aln.editLocations.size()-1] > maxMemSize) {
            maxMemSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[aln.editLocations.size()-1];
            memStart = aln.editLocations.size() - 1;
            memEnd = -1;
        }
        // cout << "MEM start: " << memStart << " end: " << memEnd << " maxMemSize: " << maxMemSize << endl;
        //extend to the left as far as possible
        start = memStart;
        end = memEnd;
        uint16_t ed = 0;
        while (start >= 0 && ed < allowedEditDistance)
        {
            start--;
            ed++;
        }
        if (ed >= allowedEditDistance) 
            ed = allowedEditDistance;

        if (start < 0) 
        {
            start = -1;
        }
        //extend to the right as far as possible
        if (ed < allowedEditDistance && end > -1)
        {
            while ((uint32_t) end < aln.editLocations.size() && ed < allowedEditDistance)
            {
                end++;
                ed++;
            }
            if (ed >= allowedEditDistance)
            {
                ed = allowedEditDistance;
                if ((uint32_t) end >= aln.editLocations.size())
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
        uint16_t maxAlnSize = 0, alnSize = 0;
        if (end == -1 && start == -1) {
            alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset;
        } else if (end == -1 && start > -1) {
            alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[start];
        } else if (end > -1 && start == -1) {
            alnSize = aln.editLocations[end];
        } else {
            alnSize = aln.editLocations[end] - aln.editLocations[start] - 1;
        }
        maxAlnSize = alnSize;
        int maxAlnStart = start, maxAlnEnd = end;
        while (start < memStart && (uint32_t) end <= aln.editLocations.size() && end > -1)
        {
            start++;
            end++;
            if ((uint32_t) end >= aln.editLocations.size())
            {
                end = -1;
            }
            if (end == -1 && start == -1) {
                alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset;
            } else if (end == -1 && start > -1) {
                alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[start];
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
            
        }
        //find the actual location of max  alignment start and end
        uint16_t maxAlnStartPos = 0, maxAlnEndPos = 0;
        if (maxAlnStart < 0){
            maxAlnStartPos = 0;
        } else {
            maxAlnStartPos = aln.editLocations[maxAlnStart] + 1;
            // for (int i = 0; i <= maxAlnStart; ++i) {
            //     aln.editLocations.erase(aln.editLocations.begin() + i);
            // }
        }
        if (maxAlnEnd == -1) {
            maxAlnEndPos = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset;
        } else {
            maxAlnEndPos = aln.editLocations[maxAlnEnd] - 1;
            // for (size_t i = maxAlnEnd; i < aln.editLocations.size(); ++i) {
            //     aln.editLocations.erase(aln.editLocations.begin() + i);
            // }
        }
        //Update aln
        // set the region start and end in R and Q
        size_t queryClips = aln.queryRegionStartPos;
        // cout << "queryClips: " << queryClips << endl;
        size_t insertionsBeforeStart = 0, deletionsBeforeStart = 0;
        for (size_t i = queryClips; i < queryClips + maxAlnStartPos; i++) {
            if (cigarStr[i] == 'I') {
                insertionsBeforeStart++;
            } else if (cigarStr[i] == 'D') {
                deletionsBeforeStart++;
            }
        }
        size_t insertionsBeforeEnd = 0, deletionsBeforeEnd = 0;
        for (size_t i = queryClips; i < queryClips + maxAlnEndPos; i++) {
            if (cigarStr[i] == 'I') {
                insertionsBeforeEnd++;
            } else if (cigarStr[i] == 'D') {
                deletionsBeforeEnd++;
            }
        }
        // cout << "maxAlnStart: " << maxAlnStart << " maxAlnEnd: " << maxAlnEnd << endl;
        // cout << "maxAlnStartPos: " << maxAlnStartPos << " maxAlnEndPos: " << maxAlnEndPos << endl;
        aln.readRegionEndPos = maxAlnEndPos + aln.readRegionStartPos - insertionsBeforeEnd;
        aln.readRegionStartPos = maxAlnStartPos + aln.readRegionStartPos - insertionsBeforeStart;
        aln.queryRegionEndPos = maxAlnEndPos + aln.queryRegionStartPos - deletionsBeforeEnd;
        aln.queryRegionStartPos = maxAlnStartPos + aln.queryRegionStartPos - deletionsBeforeStart;
        // cout << "deletionsBeforeStart: " << deletionsBeforeStart << " deletionsBeforeEnd: " << deletionsBeforeEnd << endl;
        // cout << "insertionsBeforeStart: " << insertionsBeforeStart << " insertionsBeforeEnd: " << insertionsBeforeEnd << endl;
        // cout << "aln.queryRegionEndPos: " << aln.queryRegionEndPos << " aln.queryRegionStartPos: " << aln.queryRegionStartPos << endl;
        // cout << "aln.readRegionEndPos: " << aln.readRegionEndPos << " aln.readRegionStartPos: " << aln.readRegionStartPos << endl;
        aln.matches = 0;
        aln.substitutions = 0;
        aln.inDels = 0;
        aln.partialMatchSize = maxAlnEndPos - maxAlnStartPos + 1;
        //change the cigar based on the new region
        string newCigarStr = cigarStr.substr(maxAlnStartPos + queryClips, aln.partialMatchSize);
        // cout << "newCigarStr: " << newCigarStr << endl;
        aln.cigar = convertStrToCigar(newCigarStr);
        vector<uint16_t> edits;
        parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, edits);
        aln.editDistance = aln.substitutions + aln.inDels;
    }  

    Alignment alignDifferentPenaltyScores(string query, string read, bool isForwardStran, vector<Penalty> &penalties)
    {
        Alignment bestAlignment;
        if (penalties.size() == 0)
        {
            bestAlignment.read = read;
            bestAlignment.query = query;
            if (!isForwardStran) bestAlignment.flag = 16;
            align(bestAlignment, 1, 1, 1, 1);
            //check the alignment based on the criteria
            bool criteriaCheck = checkAlingmentCriteria(bestAlignment);
            if (criteriaCheck || (!criteriaCheck && (bestAlignment.criteriaCode <= 3))){
                return bestAlignment;
            }
        }
        for (auto penalty : penalties)
        {
            Alignment aln;
            aln.read = read;
            aln.query = query;
            if (!isForwardStran) aln.flag = 16;
            align(aln, penalty.matchPenalty, penalty.mismatchPenalty, penalty.gapOpenPenalty, penalty.gapExtendPenalty);
            //check the alignment based on the criteria
            bool criteriaCheck = checkAlingmentCriteria(aln);
            if (criteriaCheck || (!criteriaCheck && (aln.criteriaCode <= 3))){
                 if (aln.partialMatchSize > bestAlignment.partialMatchSize || (aln.partialMatchSize == bestAlignment.partialMatchSize && aln.editDistance < bestAlignment.editDistance))
                    bestAlignment = aln;
            }
        }
        return bestAlignment;
    }

    void align(Alignment &aln,uint16_t matchPen, uint16_t subPen, uint16_t gapoPen, uint16_t gapextPen)
    {
        //do the alignment using smith waterman
        smithWatermanAligner(aln, matchPen, subPen, gapoPen, gapextPen);
        // cout << "---------- before -------" << endl;
        // cout << "Score: " << aln.score << endl;
        // cout << "cigar: " << aln.cigar << endl;
        // cout << "read start pos: " << aln.readRegionStartPos << endl;
        // cout << "read end pos: " << aln.readRegionEndPos << endl;
        // cout << "query start pos: " << aln.queryRegionStartPos << endl;
        // cout << "query end pos: " << aln.queryRegionEndPos << endl;
        // cout << "substitutions: " << aln.substitutions << endl;
        // cout << "inDels: " << aln.inDels << endl;
        // cout << "matches: " << aln.matches << endl;
        // cout << "editDistance: " << aln.editDistance << endl;
        // cout << "edit pos: ";
        // for(auto it = aln.editLocations.begin(); it!= aln.editLocations.end(); it++){
        //     cout << *it << " ";
        // }
        if (aln.editDistance > allowedEditDistance) {
            // exclude the additional edit distance
            findMemAndExtend(aln);
        }  
        // cout << "---------- after -------" << endl;
        // cout << "Score: " << aln.score << endl;
        // cout << "cigar: " << aln.cigar << endl;
        // cout << "read start pos: " << aln.readRegionStartPos << endl;
        // cout << "read end pos: " << aln.readRegionEndPos << endl;
        // cout << "query start pos: " << aln.queryRegionStartPos << endl;
        // cout << "query end pos: " << aln.queryRegionEndPos << endl;
        // cout << "substitutions: " << aln.substitutions << endl;
        // cout << "inDels: " << aln.inDels << endl;
        // cout << "matches: " << aln.matches << endl;
        // cout << "editDistance: " << aln.editDistance << endl;
        // cout << "edit pos: ";
        // for(auto it = aln.editLocations.begin(); it!= aln.editLocations.end(); it++){
        //     cout << *it << " ";
        // }
    }

     void dumpSam(ofstream &oSam, Alignment l)
    {
        oSam << l.queryID << '\t' << l.flag << '\t' << l.readID << '\t' << l.readRegionStartPos << '\t'
                << "*" << '\t' << l.cigar << '\t' << "*" << '\t' << "*" << '\t' << "*" << '\t' 
                << l.read << '\t' << "*" << '\t' << "NM:i:" + to_string(l.substitutions) << '\t' << "CC" << l.criteriaCode << '\n';
    }

    bool checkAlingmentCriteria(Alignment &l)
    {
        /*criteriaCode
        0000 -> accepted
        0x01 -> front region dismissed because region's starting position higher than 0 (each missed bp is 1 edit distance)
        0x02 -> back region dismissed because region's starting position lower than contigSize - regionSize (each missed bp is 1 edit distance)
        0x04 -> region starts from 100 + allowedEditDistance
        0x08 -> Either alignment len or edit distance does not match the criteria
        */
        //check the alignment len and edit distance
        if((l.editDistance > allowedEditDistance) || (((uint32_t)l.partialMatchSize < (uint32_t)regionSize)))
        {
            // cout << "0- edit distance[" << l.editDistance << "] > allowed edit distance or partial match size [" << l.partialMatchSize <<"] < region size" << endl;
            //TODO: did not check whether there is a min exact match region in the partial match region
            // if(l.queryRegionStartPos + minExactMatchLength <= regionSize || l.queryRegionEndPos - minExactMatchLength >= contigSize - regionSize)
            l.criteriaCode |= 0x08;
            return false;
        }
        bool frontRegionDismissed = false, backRegionDismissed = false;
        //check the starting position of the alignment for the front region
        if(allowedEditDistance >= l.queryRegionStartPos){
            auto editsAllwedInFrontRegion = allowedEditDistance - l.queryRegionStartPos ;
            if(l.substitutions + l.inDels > editsAllwedInFrontRegion)
            {
                frontRegionDismissed = true;
                // cout << "1- front region dismissed" << endl;
            }
        } else {
            frontRegionDismissed = true;
            // cout << "2- front region dismissed" << endl;
        }
        if(frontRegionDismissed)
            l.criteriaCode |= 0x01;
        //check the starting position of the alignment for the back region
        if(l.queryRegionStartPos + l.partialMatchSize >= contigSize - allowedEditDistance){ //query start position + alignment len should be larger than conig size - allowed edit
            auto editsAllwedInBackRegion = l.queryRegionStartPos + l.partialMatchSize - (contigSize - allowedEditDistance);
            if(l.substitutions + l.inDels > editsAllwedInBackRegion)
            {
                backRegionDismissed = true;
                // cout << "3- back region dismissed" << endl;
            }
        }else {
            backRegionDismissed = true;
            // cout << "4- back region dismissed" << endl;
        }
        if(backRegionDismissed)
            l.criteriaCode |= 0x02;
        if(frontRegionDismissed && backRegionDismissed)
            return false;
        if(l.queryRegionStartPos > (contigSize - regionSize + allowedEditDistance)){ //first bp of back region starts from 100 + allowed edit distance
            // cout << "5- first bp of back region starts from 100" << endl;
            l.criteriaCode |= 0x04;
            return false;
        }
        // cout << "6- not meet criteria" << endl;
        return true;
    }

    map<contigIndT, string> readContigsFromMap(tsl::robin_map<uint32_t, string>& reads,  set<contigIndT>& contigIdSet)
    {
        contigIndT setElement;
        unsigned int setInd = 0;
        map<contigIndT, string> readContigs;
        while (true)
        {        
            if (setInd >= 0 && setInd < contigIdSet.size())
            {
                auto it = contigIdSet.begin();
                advance(it, setInd);
                setElement = *it;
                setInd++;
                readContigs[setElement] = reads[setElement];
            }
            else
            {
                break;
            }
        }
        return readContigs;
    }

    void findPartiaMatches(tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, IndexContainer<contigIndT, contigIndT>& frontMinThCheapSeedReads, IndexContainer<contigIndT, contigIndT>& backMinThCheapSeedReads, contigIndT queryCount, bool isForwardStrand, string parmikAlignments, vector<Penalty> penalties)
    {
        ofstream pAln(parmikAlignments, ios::app);
        cout << "Starting alignment for all queries [" << (isForwardStrand ? ("fwd"):("rev")) << "]..." << endl;
        for (size_t i = 0; i < queryCount; i++)
        {
            if(i % 1000 == 0)
                cout << i << " queries processed..." << endl;
            auto itq = queries.find(i);
            if (itq == queries.end())
                continue;
            string query = itq->second;
            if (query.find('n') != string::npos || query.find('N') != string::npos)
                continue;
            // read the candidate reads of cheap k-mer filter of front
            auto frontReadSet = frontMinThCheapSeedReads.get(i);
            map<contigIndT, string> frontCandidateReads = readContigsFromMap(reads, frontReadSet);
            tsl::robin_map <uint32_t, Alignment> alignments;
            for (auto it = frontCandidateReads.begin(); it != frontCandidateReads.end(); it++)
            {
                // Alignment aln;
                // aln.query = query;
                // aln.read = it->second;
                // align(aln, 1, 1, 2, 1);
                Alignment aln = alignDifferentPenaltyScores(query, it->second, isForwardStrand, penalties);
                aln.queryID = i;
                aln.readID = it->first;
                //check the alignment based on the criteria
                // if (!isForwardStrand) aln.flag = 16;
                // bool criteriaCheck = checkAlingmentCriteria(aln);
                // if (criteriaCheck || (!criteriaCheck && (aln.criteriaCode <= 3)))
                if (aln.partialMatchSize > 0){
                    alignments.insert(make_pair(it->first, aln));
                }
            }
            auto backReadSet = backMinThCheapSeedReads.get(i);
            map<contigIndT, string> backCandidateReads = readContigsFromMap(reads, backReadSet);
            for (auto it = backCandidateReads.begin(); it != backCandidateReads.end(); it++)
            {
                // Alignment aln;
                // aln.query = query;
                // aln.read = it->second;
                // align(aln, 1, 1, 2, 1);
                Alignment aln = alignDifferentPenaltyScores(query, it->second, isForwardStrand, penalties);
                aln.queryID = i;
                aln.readID = it->first;
                //check the alignment based on the criteria
                // if (!isForwardStrand) aln.flag = 16;
                // bool criteriaCheck = checkAlingmentCriteria(aln);
                // if (criteriaCheck || (!criteriaCheck && (aln.criteriaCode <= 3))){
                if (aln.partialMatchSize > 0){
                    auto ita = alignments.find(it->first);
                    if (ita == alignments.end() || (ita != alignments.end() && ita->second.partialMatchSize < aln.partialMatchSize)){
                        alignments.insert(make_pair(it->first, aln));
                    }
                }
            }
            //dump the alignments
            size_t matchesAccepted = 0, matchesPassedWithStartingPosOverED = 0;
            for (auto it = alignments.begin(); it!= alignments.end(); it++)
            {
                if (it->second.criteriaCode <= 3) // 00, 01, 02, 03
                {
                    dumpSam(pAln, it->second);
                    matchesAccepted++;
                }
                if (it->second.criteriaCode > 0 && it->second.criteriaCode <= 3 ) // 00
                {
                    matchesPassedWithStartingPosOverED++;
                }

            }
            cout << "queryID: " << i << ", " << (isForwardStrand ? ("fwd"):("rev")) <<", total matches: " << alignments.size() << i << ", matches accepted: " << matchesAccepted << ", matches passed with starting pos over ED: " << matchesPassedWithStartingPosOverED << endl;
        }
    }

    void parseCigar(const string& cigar, size_t& matches, size_t& substitutions, size_t& inDels, vector<uint16_t> &editLocations) {
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
                    cerr << "Unsupported CIGAR operation: " << op << endl;
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
        // cout << "aln.query_begin: " << alignment.query_begin << ", aln.query_end: " << alignment.query_end << ", aln.ref_begin: " << alignment.ref_begin << ", aln.ref_end: " << alignment.ref_end << ", cigar_string: " << alignment.cigar_string << endl;
        parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, aln.editLocations);
        aln.partialMatchSize = aln.matches + aln.substitutions + aln.inDels;
        // for(auto it = aln.editLocations.begin(); it!= aln.editLocations.end(); it++){
        //     cout << *it << " ";
        // }
    }

};

#endif