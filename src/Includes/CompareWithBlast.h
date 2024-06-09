#ifndef COMAPREWITHBLAST_H
#define COMAPREWITHBLAST_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include "BlastReader.h"
#include "Aligner.h"
#include "Alignment.h"
#include <vector>
#include "Utils.h"

#define CHECK_EXACT_MATCH_CRITERION_ 0

using namespace std;

class CompareWithBlast {
private:
    double percentageIdentity;
public: 

    CompareWithBlast(double pi) : percentageIdentity(pi) {}

    void determineEditsLocationsAndType(string read, string query, vector<int>& editPos) {
        // Check if the sequences have the same length
        if (read.length() != query.length()) {
            // cout << "q : " << query << endl;
            // cout << "r : " << read << endl;
            cerr << "Error: Sequences have different lengths." << endl;
            return;
        }
        // Iterate through the sequences and compare corresponding characters
        for (size_t i = 0; i < read.length(); ++i) {
            if (read[i] != query[i]) {
                editPos.push_back(i);
            } 
        }
 
    }

    string getCigarStr(uint32_t queryS, string queryAligned, string readAligned, uint32_t contigSize)
    {
        stringstream cigar;
        for(uint32_t i = 0; i < queryS; i++) {
            cigar << 'S';
        }
        for(uint32_t i = 0; i < queryAligned.size(); i++){
            if(queryAligned[i] == readAligned[i])
                cigar << '=';
            else {
                if(queryAligned[i] == '-')
                    cigar << 'D';
                else if(readAligned[i] == '-')
                    cigar << 'I';
                else
                    cigar << 'X';
            }
        }
        int remainedClips = contigSize - (queryS + queryAligned.size());
        if(remainedClips > 0){
            for(int i = 0; i < remainedClips; i++) {
                cigar << 'S';
            }
        }
        return cigar.str();
    }

    string convertCigarToStr(const string& cigar, bool skipClips = false) {
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
                    if(!skipClips) simplifiedCigar << string(num, 'S');
                    break;
                default:
                    // Handle other CIGAR operations if needed
                    break;
            }
        }

        return simplifiedCigar.str();
    }

    bool hasConsecutiveMatches(const std::string& s, int k) {
        int count = 0;
        for (char c : s) {
            if (c == '=') {
                count++;
                if (count >= k) {
                    return true;
                }
            } else {
                count = 0; // Reset count if the current character is not '='
            }
        }
        return false;
    }

    uint32_t getMatchesCount(string cigarStr) {
        uint32_t cnt = 0;
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == '=') {
                cnt++;
            }
        }
        return cnt;
    }

    string trimClips(string cigarStr){
        string trimmedCigarStr = "";
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'S' || cigarStr[i] == 'H') {
                continue;
            }
            trimmedCigarStr += cigarStr[i];
        }
        return trimmedCigarStr;
    }

    bool checkIdentityPercentange(string cigarStr, bool isTrimClips = false) {
        if (isTrimClips) cigarStr = trimClips(cigarStr);
        uint32_t matches = getMatchesCount(cigarStr);
        double identity =  (double) matches / (double) cigarStr.size();
        if(DEBUG_MODE) cout << "matches: " << matches << ", cigarStr.size : " << cigarStr.size() << ", identity: " << identity << endl;
        if(identity < percentageIdentity)
            return false;
        return true;
    }

    void comparePmWithBlast(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, const uint32_t queryCount, const string& outputAddress)
    {
        Utilities<uint16_t> utils;
        //read the blast file
        BlastReader blastReader(cfg.otherToolOutputFileAddress);
        vector<Alignment> blastAlignments, parmikAlignments;
        blastReader.parseFile(queryCount, blastAlignments);
        //read the parmik file
        SamReader parmikSam(cfg.baselineBaseAddress);
        parmikSam.parseFile(queryCount, parmikAlignments, false);
        ofstream outputFile(outputAddress);
        if(!outputFile.is_open()){
            cout << "error openin file: " << outputAddress << endl;
            return;
        }

        uint32_t blastFN = 0, parmikFN = 0, totalTN = 0, numberOfQueryContainN = 0;
        uint64_t blastOutperform = 0, parmikOutperform = 0, equalPerformance = 0;
        uint64_t sameReadBlastOutperform = 0, differentReadBlastOutperform = 0, sameReadPmOutperform = 0, differentReadPmOutperform = 0, sameReadEqual = 0, differentReadEqual = 0;
        multiset<uint16_t> sameReadBlastOutperformBps, differentReadBlastOutperformBps, sameReadPmOutperformBps, differentReadPmOutperformBps;
        for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            if(queryInd % 1000 == 0) {
                cout << "queries processed: " << queryInd << " / " << queryCount << endl;
            }
            string query = queries[queryInd];
            if(query.find('N') != string::npos || query.find('n') != string::npos)
            {
                cout << "query contains N!" << endl;
                numberOfQueryContainN++; 
                // if(REPORT_ALN_PER_Q) alnPerQ << queryInd << " 0 0"<< endl;
                continue;
            }
            vector<Alignment> query_blastAlignments;
            for (const Alignment& aln : blastAlignments) 
            {
                if ((uint32_t)aln.queryID == queryInd) {
                    if(CHECK_EXACT_MATCH_CRITERION_) {
                        string cigarStr = getCigarStr(aln.queryRegionStartPos, aln.alignedQuery, aln.alignedRead, cfg.contigSize);
                        if(hasConsecutiveMatches(cigarStr, cfg.kmerLength))
                            query_blastAlignments.push_back(aln);
                    } else {
                        query_blastAlignments.push_back(aln);
                    }
                }
            }
            size_t blastReadPerQuery = query_blastAlignments.size();

            vector<Alignment> query_parmikAlignments;
            for (const Alignment& aln : parmikAlignments) 
            {
                if ((uint32_t)aln.queryID == queryInd) {
                    if(CHECK_EXACT_MATCH_CRITERION_) {
                        string cigarStr = convertCigarToStr(aln.cigar, true);
                        if(hasConsecutiveMatches(cigarStr, cfg.kmerLength))
                            query_parmikAlignments.push_back(aln);
                    } else {
                        query_parmikAlignments.push_back(aln);
                    }
                }
            }
            size_t parmikReadPerQuery = query_parmikAlignments.size();

            // cout << "# of reads found by BLAST : " << blastReadPerQuery << endl;
            // cout << "# of reads found by PARMIK : " << parmikReadPerQuery << endl;
            vector<uint32_t> blastTPReadIds;
            if(blastReadPerQuery == 0 && parmikReadPerQuery == 0){
                //TN
                totalTN++;
                // cmp << "TN for BLAST" << endl;
            } else if(blastReadPerQuery == 0 && parmikReadPerQuery > 0) {
                //FN
                blastFN++;
            } else if(blastReadPerQuery > 0 && parmikReadPerQuery == 0) {
                //FN
                parmikFN++;
            } else { // both > 0
                //get the best PARMIK alignment
                Alignment bestAlnBlast;
                Alignment bestAlnPm;
                for (auto it = query_blastAlignments.begin(); it != query_blastAlignments.end(); it++) {
                    Alignment blastAlnn = (*it);
                    if (blastAlnn.partialMatchSize > bestAlnBlast.partialMatchSize) // only exact matches
                    {
                        bestAlnBlast = blastAlnn;
                    }
                    else if (blastAlnn.partialMatchSize == bestAlnBlast.partialMatchSize)
                    {
                        if (blastAlnn.inDels + blastAlnn.substitutions < bestAlnBlast.inDels + bestAlnBlast.substitutions) // InDel has the same wight as substitution
                        {
                            bestAlnBlast = blastAlnn;
                        }
                    }
                }
                Alignment parmikAlnForBlastSameReadID;
                for (auto it = query_parmikAlignments.begin(); it != query_parmikAlignments.end(); it++) {
                    Alignment pmAlnn = (*it);
                    if(bestAlnBlast.partialMatchSize > 0 && bestAlnBlast.readID == pmAlnn.readID) {
                        parmikAlnForBlastSameReadID = pmAlnn;
                    }
                    if (pmAlnn.matches + pmAlnn.inDels + pmAlnn.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                    {
                        bestAlnPm = pmAlnn;
                    }
                    else if (pmAlnn.matches + pmAlnn.inDels + pmAlnn.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions)
                    {
                        if (pmAlnn.inDels + pmAlnn.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) // InDel has the same wight as substitution
                        {
                            bestAlnPm = pmAlnn;
                        }
                    }
                }
                bool foundSameRead = false;
                if(bestAlnBlast.readID == bestAlnPm.readID) {
                    foundSameRead = true;
                }
                if (bestAlnBlast.partialMatchSize > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                {
                    blastOutperform++;
                    if (foundSameRead){
                        sameReadBlastOutperform++;
                        sameReadBlastOutperformBps.insert(bestAlnBlast.partialMatchSize - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                        outputFile << "sameReadBlastOutperform";
                    } else {
                        differentReadBlastOutperform++;
                        differentReadBlastOutperformBps.insert(bestAlnBlast.partialMatchSize - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                        outputFile << "parmikAlnForBlastSameReadID: " << parmikAlnForBlastSameReadID.cigar << ", M: " << parmikAlnForBlastSameReadID.matches << ", S: " << parmikAlnForBlastSameReadID.substitutions << ", InDels: " << parmikAlnForBlastSameReadID.inDels << endl;
                        outputFile << "differentReadBlastOutperform";
                    }
                    outputFile << " (larger aln_length) for [" << bestAlnBlast.queryID << ", " << bestAlnBlast.readID << "]:, blast cigar: " 
                    << getCigarStr(bestAlnBlast.queryRegionStartPos, bestAlnBlast.alignedQuery, bestAlnBlast.alignedRead, cfg.contigSize) << ", M: " << bestAlnBlast.partialMatchSize - bestAlnBlast.inDels - bestAlnBlast.substitutions << ", S: " << bestAlnBlast.substitutions << ", InDels: " << bestAlnBlast.inDels
                    << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                } else if (bestAlnBlast.partialMatchSize < bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                {
                    parmikOutperform++;
                    if (foundSameRead){
                        sameReadPmOutperform++;
                        sameReadPmOutperformBps.insert((bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) - bestAlnBlast.partialMatchSize);
                        outputFile << "sameReadPmOutperform";
                    } else {
                        differentReadPmOutperform++;
                        differentReadPmOutperformBps.insert((bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) - bestAlnBlast.partialMatchSize);
                        outputFile << "parmikAlnForBlastSameReadID: " << parmikAlnForBlastSameReadID.cigar << ", M: " << parmikAlnForBlastSameReadID.matches << ", S: " << parmikAlnForBlastSameReadID.substitutions << ", InDels: " << parmikAlnForBlastSameReadID.inDels << endl;
                        outputFile << "differentReadPmOutperform";
                    }
                    outputFile << " (larger aln_length) for [" << bestAlnBlast.queryID << ", " << bestAlnBlast.readID << "]:, blast cigar: " 
                    << getCigarStr(bestAlnBlast.queryRegionStartPos, bestAlnBlast.alignedQuery, bestAlnBlast.alignedRead, cfg.contigSize) << ", M: " << bestAlnBlast.partialMatchSize - bestAlnBlast.inDels - bestAlnBlast.substitutions << ", S: " << bestAlnBlast.substitutions << ", InDels: " << bestAlnBlast.inDels
                    << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                } else if (bestAlnBlast.partialMatchSize == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions)
                {
                    if (bestAlnBlast.inDels + bestAlnBlast.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) // InDel has the same wight as substitution
                    {
                        blastOutperform++;
                        if (foundSameRead){
                            sameReadBlastOutperform++;
                            sameReadBlastOutperformBps.insert(bestAlnBlast.partialMatchSize - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                            outputFile << "sameReadBlastOutperform";
                        } else {
                            differentReadBlastOutperform++;
                            differentReadBlastOutperformBps.insert(bestAlnBlast.partialMatchSize - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                            outputFile << "parmikAlnForBlastSameReadID: " << parmikAlnForBlastSameReadID.cigar << ", M: " << parmikAlnForBlastSameReadID.matches << ", S: " << parmikAlnForBlastSameReadID.substitutions << ", InDels: " << parmikAlnForBlastSameReadID.inDels << endl;
                            outputFile << "differentReadBlastOutperform";
                        }
                        outputFile << " (fewer edits) for [" << bestAlnBlast.queryID << ", " << bestAlnBlast.readID << "]:, blast cigar: " 
                        << getCigarStr(bestAlnBlast.queryRegionStartPos, bestAlnBlast.alignedQuery, bestAlnBlast.alignedRead, cfg.contigSize) << ", M: " << bestAlnBlast.partialMatchSize - bestAlnBlast.inDels - bestAlnBlast.substitutions << ", S: " << bestAlnBlast.substitutions << ", InDels: " << bestAlnBlast.inDels
                        << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                    } else if (bestAlnBlast.inDels + bestAlnBlast.substitutions > bestAlnPm.inDels + bestAlnPm.substitutions)
                    {
                        parmikOutperform++;
                        if (foundSameRead){
                            sameReadPmOutperform++;
                            sameReadPmOutperformBps.insert((bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) - bestAlnBlast.partialMatchSize);
                         outputFile << "sameReadPmOutperform";
                        } else {
                            differentReadPmOutperform++;
                            differentReadPmOutperformBps.insert((bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) - bestAlnBlast.partialMatchSize);
                            outputFile << "parmikAlnForBlastSameReadID: " << parmikAlnForBlastSameReadID.cigar << ", M: " << parmikAlnForBlastSameReadID.matches << ", S: " << parmikAlnForBlastSameReadID.substitutions << ", InDels: " << parmikAlnForBlastSameReadID.inDels << endl;
                            outputFile << "differentReadPmOutperform";
                        }
                        outputFile << " (fewer edits) for [" << bestAlnBlast.queryID << ", " << bestAlnBlast.readID << "]:, blast cigar: " 
                        << getCigarStr(bestAlnBlast.queryRegionStartPos, bestAlnBlast.alignedQuery, bestAlnBlast.alignedRead, cfg.contigSize) << ", M: " << bestAlnBlast.partialMatchSize - bestAlnBlast.inDels - bestAlnBlast.substitutions << ", S: " << bestAlnBlast.substitutions << ", InDels: " << bestAlnBlast.inDels
                        << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                    } else{
                        equalPerformance++;
                        if (foundSameRead){
                            sameReadEqual++;
                            outputFile << "sameReadEqual";
                        } else {
                            differentReadEqual++;
                            outputFile << "differentReadEqual";
                        }
                        outputFile << " for [" << bestAlnBlast.queryID << ", " << bestAlnBlast.readID << "]:, blast cigar: " 
                        << getCigarStr(bestAlnBlast.queryRegionStartPos, bestAlnBlast.alignedQuery, bestAlnBlast.alignedRead, cfg.contigSize) << ", M: " << bestAlnBlast.partialMatchSize - bestAlnBlast.inDels - bestAlnBlast.substitutions << ", S: " << bestAlnBlast.substitutions << ", InDels: " << bestAlnBlast.inDels
                        << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                    }
                }
            }
        }
        outputFile << "<<<<<<<<<<<<<<<<<<<<<<<<Final Results>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        outputFile << "Total Blast TN : " << totalTN << endl;
        outputFile << "Total Blast FN : " << blastFN << endl;
        outputFile << "Total Parmik FNN : " << parmikFN << endl;
        outputFile << "Total Blast Outperform : " << blastOutperform << endl;
        outputFile << "Total Parmik Outperform : " << parmikOutperform << endl;
        outputFile << "Total Blast Equal : " << equalPerformance << endl;
        outputFile << "Same Read Blast Outperform : " << sameReadBlastOutperform << endl;
        outputFile << "Different Read Blast Outperform : " << differentReadBlastOutperform << endl;
        outputFile << "Same Read Parmik Outperform : " << sameReadPmOutperform << endl;
        outputFile << "Different Read Parmik Outperform : " << differentReadPmOutperform << endl;
        outputFile << "Same Read Equal : " << sameReadEqual << endl;
        outputFile << "Different Read Equal : " << differentReadEqual << endl;
        pair<uint16_t, uint16_t> sameReadBlastOutperformBpsTuple = utils.calculateStatistics2(sameReadBlastOutperformBps);
        printf("No. of Bp (same read) BWA outperforms => [average: %d, median: %d]\n", get<0>(sameReadBlastOutperformBpsTuple), get<1>(sameReadBlastOutperformBpsTuple));
        pair<uint16_t, uint16_t> differentReadBlastOutperformBpsTuple = utils.calculateStatistics2(differentReadBlastOutperformBps);
        printf("No. of Bp (different read) BWA outperforms => [average: %d, median: %d]\n", get<0>(differentReadBlastOutperformBpsTuple), get<1>(differentReadBlastOutperformBpsTuple));
        pair<uint16_t, uint16_t> sameReadPmOutperformBpsTuple = utils.calculateStatistics2(sameReadPmOutperformBps);
        printf("No. of Bp (same read) PARMIK outperforms => [average: %d, median: %d]\n", get<0>(sameReadPmOutperformBpsTuple), get<1>(sameReadPmOutperformBpsTuple));
        pair<uint16_t, uint16_t> differentReadPmOutperformBpsTuple = utils.calculateStatistics2(differentReadPmOutperformBps);
        printf("No. of Bp (different read) PARMIK outperforms => [average: %d, median: %d]\n", get<0>(differentReadPmOutperformBpsTuple), get<1>(differentReadPmOutperformBpsTuple));
    }
};

#endif