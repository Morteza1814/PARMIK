#ifndef ALIGNER_H
#define ALIGNER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include <chrono>
#include "Utils.h"
#include "sw/ssw_cpp.h"
#include "PostFilter.h"
#include "Alignment.h"
#include <omp.h>

using namespace std;

template <typename contigIndT>
class Aligner {
private:
    uint16_t regionSize;
    uint16_t allowedEditDistance;
    uint16_t contigSize;
    uint16_t k_;
    uint16_t minNumExactMatchKmer;
    double percentageIdentity;
    bool isSecondChanceOff;
    int num_threads = 1;
public:
    Aligner(uint16_t R, uint16_t a, uint16_t c, uint16_t k, uint16_t m, double i, bool sc, int th) : regionSize(R), allowedEditDistance(a), contigSize(c), k_(k), minNumExactMatchKmer(m), percentageIdentity(i), isSecondChanceOff(sc), num_threads(th) {}

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

    string convertStrToCigar(const string cigarStr, uint32_t queryS, uint32_t queryE) {
        stringstream cigar;
        char prev = cigarStr[0];
        uint16_t num = 0;
        // cout << "queryS: " << queryS << ", queryE: " << queryE << endl;
        if (prev != 'S' && prev != 'H' && queryS > 0) {
            cigar << queryS << 'S';
        }
        for (size_t i = 1; i < cigarStr.length(); ++i) {
            char cur = cigarStr[i];
            num++;
            if (prev != cur)
            {
                if (prev == 'S' || prev == 'H') {
                    num += queryS;
                }
                cigar << num << prev;
                num = 0;
                prev = cur;
            }
        }
        num++;
        int remainedClips = contigSize - queryE - 1;
        if (remainedClips > 0 && (prev == 'S' || prev == 'H')) {
            num += remainedClips;
            remainedClips = 0;
        }
        cigar << num << prev;
        if (remainedClips > 0) {
            cigar << remainedClips << 'S';
        }
        return cigar.str();
    }

    // void findMemAndExtend(Alignment &aln){
    //     int start = -1, end = -1;
    //     int memStart = -1, memEnd = -1;
    //     string cigarStr = convertCigarToStr(aln.cigar);
    //     // cout << "CIGAR str: " << cigarStr << endl;
    //     // insertions adds to the read positions and deletions adds to the query positions
    //     size_t alignmentOffset = (cigarStr.size() > contigSize) ? cigarStr.size() - contigSize : 0; // if the alignment size goes over 150bp
    //     //find MEM
    //     uint16_t memSize = 0, maxMemSize = 0;
    //     for (size_t i = 0; i < aln.editLocations.size(); ++i) {
    //         if (i == 0) {
    //             end = i;
    //             memSize = aln.editLocations[i];
    //         } else {
    //             start = i - 1;
    //             end = i;
    //             memSize = aln.editLocations[i] - aln.editLocations[i-1] - 1;
    //         }
    //         if (memSize > maxMemSize) {
    //             maxMemSize = memSize;
    //             memStart = start;
    //             memEnd = end;
    //         }
    //     }
    //     //(aln.queryRegionEndPos - aln.queryRegionStartPos) to make the beginning of the region 0
    //     if (aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[aln.editLocations.size()-1] > maxMemSize) {
    //         maxMemSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[aln.editLocations.size()-1];
    //         memStart = aln.editLocations.size() - 1;
    //         memEnd = -1;
    //     }
    //     // cout << "MEM start: " << memStart << " end: " << memEnd << " maxMemSize: " << maxMemSize << endl;
    //     //extend to the left as far as possible
    //     start = memStart;
    //     end = memEnd;
    //     uint16_t ed = 0;
    //     while (start >= 0 && ed < allowedEditDistance)
    //     {
    //         start--;
    //         ed++;
    //     }
    //     if (ed >= allowedEditDistance) 
    //         ed = allowedEditDistance;

    //     if (start < 0) 
    //     {
    //         start = -1;
    //     }
    //     //extend to the right as far as possible
    //     if (ed < allowedEditDistance && end > -1)
    //     {
    //         while ((uint32_t) end < aln.editLocations.size() && ed < allowedEditDistance)
    //         {
    //             end++;
    //             ed++;
    //         }
    //         if (ed >= allowedEditDistance)
    //         {
    //             ed = allowedEditDistance;
    //             if ((uint32_t) end >= aln.editLocations.size())
    //             {
    //                 end = -1;
    //             }
    //         }
    //         else
    //         {
    //             end = -1;
    //         }
    //     }
    //     // calculate the largest aligment
    //     uint16_t maxAlnSize = 0, alnSize = 0;
    //     if (end == -1 && start == -1) {
    //         alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset;
    //     } else if (end == -1 && start > -1) {
    //         alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[start];
    //     } else if (end > -1 && start == -1) {
    //         alnSize = aln.editLocations[end];
    //     } else {
    //         alnSize = aln.editLocations[end] - aln.editLocations[start] - 1;
    //     }
    //     maxAlnSize = alnSize;
    //     int maxAlnStart = start, maxAlnEnd = end;
    //     while (start < memStart && (uint32_t) end <= aln.editLocations.size() && end > -1)
    //     {
    //         start++;
    //         end++;
    //         if ((uint32_t) end >= aln.editLocations.size())
    //         {
    //             end = -1;
    //         }
    //         if (end == -1 && start == -1) {
    //             alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset;
    //         } else if (end == -1 && start > -1) {
    //             alnSize = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset - aln.editLocations[start];
    //         } else if (end > -1 && start == -1) {
    //             alnSize = aln.editLocations[end];
    //         } else {
    //             alnSize = aln.editLocations[end] - aln.editLocations[start] - 1;
    //         }
    //         if (alnSize > maxAlnSize)
    //         {
    //             maxAlnSize = alnSize;
    //             maxAlnStart = start;
    //             maxAlnEnd = end;
    //         }
            
    //     }
    //     //find the actual location of max  alignment start and end
    //     uint16_t maxAlnStartPos = 0, maxAlnEndPos = 0;
    //     if (maxAlnStart < 0){
    //         maxAlnStartPos = 0;
    //     } else {
    //         maxAlnStartPos = aln.editLocations[maxAlnStart] + 1;
    //         // for (int i = 0; i <= maxAlnStart; ++i) {
    //         //     aln.editLocations.erase(aln.editLocations.begin() + i);
    //         // }
    //     }
    //     if (maxAlnEnd == -1) {
    //         maxAlnEndPos = aln.queryRegionEndPos - aln.queryRegionStartPos + alignmentOffset;
    //     } else {
    //         maxAlnEndPos = aln.editLocations[maxAlnEnd] - 1;
    //         // for (size_t i = maxAlnEnd; i < aln.editLocations.size(); ++i) {
    //         //     aln.editLocations.erase(aln.editLocations.begin() + i);
    //         // }
    //     }
    //     //Update aln
    //     // set the region start and end in R and Q
    //     size_t queryClips = aln.queryRegionStartPos;
    //     // cout << "queryClips: " << queryClips << endl;
    //     size_t insertionsBeforeStart = 0, deletionsBeforeStart = 0;
    //     for (size_t i = queryClips; i < queryClips + maxAlnStartPos; i++) {
    //         if (cigarStr[i] == 'I') {
    //             insertionsBeforeStart++;
    //         } else if (cigarStr[i] == 'D') {
    //             deletionsBeforeStart++;
    //         }
    //     }
    //     size_t insertionsBeforeEnd = 0, deletionsBeforeEnd = 0;
    //     for (size_t i = queryClips; i < queryClips + maxAlnEndPos; i++) {
    //         if (cigarStr[i] == 'I') {
    //             insertionsBeforeEnd++;
    //         } else if (cigarStr[i] == 'D') {
    //             deletionsBeforeEnd++;
    //         }
    //     }
    //     // cout << "maxAlnStart: " << maxAlnStart << " maxAlnEnd: " << maxAlnEnd << endl;
    //     // cout << "maxAlnStartPos: " << maxAlnStartPos << " maxAlnEndPos: " << maxAlnEndPos << endl;
    //     aln.readRegionEndPos = maxAlnEndPos + aln.readRegionStartPos - insertionsBeforeEnd;
    //     aln.readRegionStartPos = maxAlnStartPos + aln.readRegionStartPos - insertionsBeforeStart;
    //     aln.queryRegionEndPos = maxAlnEndPos + aln.queryRegionStartPos - deletionsBeforeEnd;
    //     aln.queryRegionStartPos = maxAlnStartPos + aln.queryRegionStartPos - deletionsBeforeStart;
    //     // cout << "deletionsBeforeStart: " << deletionsBeforeStart << " deletionsBeforeEnd: " << deletionsBeforeEnd << endl;
    //     // cout << "insertionsBeforeStart: " << insertionsBeforeStart << " insertionsBeforeEnd: " << insertionsBeforeEnd << endl;
    //     // cout << "aln.queryRegionEndPos: " << aln.queryRegionEndPos << " aln.queryRegionStartPos: " << aln.queryRegionStartPos << endl;
    //     // cout << "aln.readRegionEndPos: " << aln.readRegionEndPos << " aln.readRegionStartPos: " << aln.readRegionStartPos << endl;
    //     aln.matches = 0;
    //     aln.substitutions = 0;
    //     aln.inDels = 0;
    //     aln.partialMatchSize = maxAlnEndPos - maxAlnStartPos + 1;
    //     //change the cigar based on the new region
    //     string newCigarStr = cigarStr.substr(maxAlnStartPos + queryClips, aln.partialMatchSize);
    //     // cout << "newCigarStr: " << newCigarStr << endl;
    //     aln.cigar = convertStrToCigar(newCigarStr, aln.queryRegionStartPos, aln.queryRegionEndPos);
    //     vector<uint16_t> edits;
    //     Utilities<uint32_t> util;
    //     util.parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, edits);
    //     aln.editDistance = aln.substitutions + aln.inDels;
    // }  

    Alignment alignDifferentPenaltyScores(string query, string read, uint32_t queryID, uint32_t readID, bool isForwardStran, vector<Penalty> &penalties)
    {
        Alignment bestAlignment;
        PostFilter pf(regionSize, k_, contigSize, minNumExactMatchKmer, percentageIdentity, isSecondChanceOff);
        if (penalties.size() == 0)
        {
            bestAlignment.read = read;
            bestAlignment.query = query;
            bestAlignment.readID = readID;
            bestAlignment.queryID = queryID;
            if (!isForwardStran) bestAlignment.flag = 16;
            align(bestAlignment, 1, 1, 1, 1);
            //check the alignment based on the criteria
            bool criteriaCheck = pf.checkAndUpdateBasedOnAlingmentCriteria(bestAlignment);
            if (criteriaCheck ){
                return bestAlignment;
            }
        }
        vector<string> repeatedCigars;
        for (auto penalty : penalties)
        {
            if(DEBUG_MODE) 
            {
                cout << "<<<<<<<<<<<<Match: " << penalty.matchPenalty << " Mismatch: " << penalty.mismatchPenalty << " GapOpen: " << penalty.gapOpenPenalty << " GapExtend: " << penalty.gapExtendPenalty << ">>>>>>>>>>>>>" << endl;
            }
            Alignment aln;
            aln.readID = readID;
            aln.queryID = queryID;
            aln.read = read;
            aln.query = query;
            if (!isForwardStran) aln.flag = 16;
            // cout << "penalties: " << penalty.matchPenalty <<  penalty.mismatchPenalty <<  penalty.gapOpenPenalty << penalty.gapExtendPenalty << endl;
            align(aln, penalty.matchPenalty, penalty.mismatchPenalty, penalty.gapOpenPenalty, penalty.gapExtendPenalty);
            if(find(repeatedCigars.begin(), repeatedCigars.end(), aln.cigar) != repeatedCigars.end())
                continue;
            else
                repeatedCigars.push_back(aln.cigar);
            //check the alignment based on the criteria
            bool criteriaCheck = pf.checkAndUpdateBasedOnAlingmentCriteria(aln);
            if (criteriaCheck){
                gAlignmentFoundWithPolish++;
                if (aln.partialMatchSize > bestAlignment.partialMatchSize || (aln.partialMatchSize == bestAlignment.partialMatchSize && aln.editDistance < bestAlignment.editDistance))
                {
                    gAlignmentFoundWithPolishLargerThanBest++;
                    bestAlignment = aln;
                }
            }
        }
        return bestAlignment;
    }

    void align(Alignment &aln,uint16_t matchPen, uint16_t subPen, uint16_t gapoPen, uint16_t gapextPen)
    {
        //do the alignment using smith waterman
        smithWatermanAligner(aln, matchPen, subPen, gapoPen, gapextPen);
        // cout << "penalties: " << matchPen <<  subPen <<  gapoPen << gapextPen << endl;
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
        // cerr << "queryID: " << l.queryID << ", readID: "<< l.readID << "readRegionStartPos: " << l.readRegionStartPos << endl;
        assert((l.readRegionStartPos >= 0 && l.readRegionStartPos < contigSize) && "wrong readRegionStartPos");
        oSam << l.queryID << '\t' << l.flag << '\t' << l.readID << '\t' << l.readRegionStartPos << '\t'
                << "*" << '\t' << l.cigar << '\t' << "*" << '\t' << "*" << '\t' << "*" << '\t' 
                << l.read << '\t' << "*" << '\t' << "NM:i:" + to_string(l.substitutions) << '\t' << "CC" << l.criteriaCode << '\n';
    }

    map<contigIndT, string> readContigsFromMap(tsl::robin_map<uint32_t, string>& reads,  unordered_set<contigIndT>& contigIdSet)
    {
        // contigIndT setElement;
        // unsigned int setInd = 0;
        map<contigIndT, string> readContigs;
        for (auto it = contigIdSet.begin(); it!= contigIdSet.end(); it++)
        {
            auto itt = reads.find(*it);
            if (itt == reads.end())
                assert(false && "contig not found in reads");
            readContigs[*it] = itt->second;
        }
        // while (true)
        // {        
        //     if (setInd >= 0 && setInd < contigIdSet.size())
        //     {
        //         auto it = contigIdSet.begin();
        //         advance(it, setInd);
        //         setElement = *it;
        //         setInd++;
        //         readContigs[setElement] = reads[setElement];
        //     }
        //     else
        //     {
        //         break;
        //     }
        // }
        return readContigs;
    }

    void findPartiaMatches(tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, Container<contigIndT, contigIndT>& minThCheapSeedReads, contigIndT queryCount, bool isForwardStrand, string parmikAlignments, vector<Penalty> penalties)
    {
        ofstream pAln;
        if (isForwardStrand) {
            pAln.open(parmikAlignments);
        } else {
            pAln.open(parmikAlignments, std::ios::app);
        }
        multiset<uint32_t> matchesPerQuery;
        cout << "Starting alignment for all queries [" << (isForwardStrand ? ("fwd"):("rev")) << "]..." << endl;
        omp_set_num_threads(num_threads);
        #pragma omp parallel for
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
            auto readSet = minThCheapSeedReads.get(i);
            map<contigIndT, string> candidateReads = readContigsFromMap(reads, readSet);
            tsl::robin_map <uint32_t, Alignment> alignments;
            bool haveAlign = false;
            auto start = chrono::high_resolution_clock::now();
            for (auto it = candidateReads.begin(); it != candidateReads.end(); it++)
            {
                Alignment aln = alignDifferentPenaltyScores(query, it->second, i, it->first, isForwardStrand, penalties);
                if (aln.partialMatchSize > 0){
                    if (num_threads == 1) gAlignmentDumped++;
                    haveAlign = true;
                    alignments.insert(make_pair(it->first, aln));
                }
            }
            if(num_threads == 1 && haveAlign) gQueriesHaveAtLeastOneAlignment++;
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            if (num_threads == 1) 
                cout << "query [" << i << "]: # of matches found " << alignments.size() << " and it took: " << static_cast<double>(duration.count()) / 1'000'000.0 << "ms" << endl;
            //dump the alignments
            uint32_t matchesAccepted = 0;
            #pragma omp critical
            {
                for (auto it = alignments.begin(); it!= alignments.end(); it++)
                {
                    dumpSam(pAln, it->second);
                    matchesAccepted++;
                }
                matchesPerQuery.insert(alignments.size());
            }
            // cout << "queryID: " << i << ", " << (isForwardStrand ? ("fwd"):("rev")) <<", total matches: " << alignments.size() << i << ", matches accepted: " << matchesAccepted << ", matches passed with starting pos over ED: " << matchesPassedWithStartingPosOverED << endl;
        }
        Utilities<uint32_t> util; 
        tuple<uint32_t, uint32_t, uint32_t> matchesPerQueryTuple = util.calculateStatistics(matchesPerQuery);
        printf("PARMIK's matches per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(matchesPerQueryTuple), get<1>(matchesPerQueryTuple), get<2>(matchesPerQueryTuple));
    }

    void smithWatermanAligner(Alignment &aln, uint16_t matchPen, uint16_t subPen, uint16_t gapoPen, uint16_t gapextPen)
    {
        Utilities<uint32_t> util;
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
        util.parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, aln.editLocations);
        aln.partialMatchSize = aln.matches + aln.substitutions + aln.inDels;
        // for(auto it = aln.editLocations.begin(); it!= aln.editLocations.end(); it++){
        //     cout << *it << " ";
        // }
    }

};

#endif