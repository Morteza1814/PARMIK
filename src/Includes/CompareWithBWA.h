#ifndef COMAPREWITHBWA_H
#define COMAPREWITHBWA_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include "Aligner.h"
#include "Alignment.h"

#define CHECK_EXACT_MATCH_CRITERION__ 0

using namespace std;

class ComparatorWithBWA {
private:
    double percentageIdentity;
public: 

    ComparatorWithBWA(double pi) : percentageIdentity(pi) {}


    string convertCigarToStr(const string& cigar, bool skipClips = false) {
        //TODO: check if this is correct
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
                case 'M':
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
                case 'H':
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
            if (c == 'M') {
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
            if (cigarStr[i] == 'M') {
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

    bool checkIdentityPercentange(const Alignment &aln) {
        // if (isTrimClips) cigarStr = trimClips(cigarStr);
        // uint32_t matches = getMatchesCount(cigarStr);
        double identity =  (double) aln.matches / (double) aln.partialMatchSize;
        if(DEBUG_MODE) cout << "matches: " << aln.matches << ", cigarStr.size : " << aln.partialMatchSize << ", identity: " << identity << endl;
        if(identity < percentageIdentity)
            return false;
        return true;
    }

    void comparePmWithBwa(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, const uint32_t queryCount, const string& outputAddress)
    {
        vector<Alignment> bwaAlignments, parmikAlignments;
        //read the bwa file
        SamReader bwaSam(cfg.otherToolOutputFileAddress);
        bwaSam.parseFile(queryCount, bwaAlignments, false, true);
        //read the parmik file
        SamReader parmikSam(cfg.baselineBaseAddress);
        parmikSam.parseFile(queryCount, parmikAlignments, false, false);
        ofstream outputFile(outputAddress);
        if(!outputFile.is_open()){
            cout << "error openin file: " << outputAddress << endl;
            return;
        }

        uint32_t bwaFN = 0, parmikFN = 0, totalTN = 0, numberOfQueryContainN = 0;
        uint64_t bwaOutperform = 0, parmikOutperform = 0, equalPerformance = 0;
        uint64_t sameReadBwaOutperform = 0, differentReadBwaOutperform = 0, sameReadPmOutperform = 0, differentReadPmOutperform = 0, sameReadEqual = 0, differentReadEqual = 0;
        uint64_t bwaLowPI_FN = 0, bwaLowPIOutperformBestPm = 0;
        for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            vector<Alignment> bwaLowPI_Alignments;
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
            vector<Alignment> query_bwaAlignments;
            for (const Alignment& aln : bwaAlignments) 
            {
                if ((uint32_t)aln.queryID == queryInd) {
                    if (checkIdentityPercentange(aln)) {
                        if(CHECK_EXACT_MATCH_CRITERION__) {
                            string cigarStr = convertCigarToStr(aln.cigar, true);
                            if(hasConsecutiveMatches(cigarStr, cfg.kmerLength))
                                query_bwaAlignments.push_back(aln);
                        } else {
                            query_bwaAlignments.push_back(aln);
                        }
                    } else {
                        bwaLowPI_Alignments.push_back(aln);
                        bwaLowPI_FN++;
                    }
                }
            }
            size_t bwaReadPerQuery = query_bwaAlignments.size();

            vector<Alignment> query_parmikAlignments;
            for (const Alignment& aln : parmikAlignments) 
            {
                if ((uint32_t)aln.queryID == queryInd) {
                    if(CHECK_EXACT_MATCH_CRITERION__) {
                        string cigarStr = convertCigarToStr(aln.cigar, true);
                        if(hasConsecutiveMatches(cigarStr, cfg.kmerLength))
                            query_parmikAlignments.push_back(aln);
                    } else {
                        query_parmikAlignments.push_back(aln);
                    }
                }
            }
            size_t parmikReadPerQuery = query_parmikAlignments.size();

            // cout << "# of reads found by BWA : " << bwaReadPerQuery << endl;
            // cout << "# of reads found by PARMIK : " << parmikReadPerQuery << endl;
            bool bothHasTP = false, bwaHasLowPIFN = false;
            vector<uint32_t> bwaTPReadIds;
            if(bwaReadPerQuery == 0 && parmikReadPerQuery == 0){
                //TN
                totalTN++;
                // cmp << "TN for BWA" << endl;
            } else if(bwaReadPerQuery == 0 && parmikReadPerQuery > 0) {
                //FN
                bwaFN++;
                if(bwaLowPI_Alignments.size() > 0) {
                    bwaHasLowPIFN = true;
                }
            } else if(bwaReadPerQuery > 0 && parmikReadPerQuery == 0) {
                //FN
                parmikFN++;
            } else { // both > 0
                bothHasTP = true;
            }
            if(bothHasTP || bwaHasLowPIFN) {
                //get the best PARMIK alignment
                Alignment bestAlnBwa;
                Alignment bestAlnPm;
                for (auto it = query_parmikAlignments.begin(); it != query_parmikAlignments.end(); it++) {
                    Alignment pmAlnn = (*it);
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
                if(!bwaHasLowPIFN) {
                    for (auto it = query_bwaAlignments.begin(); it != query_bwaAlignments.end(); it++) {
                        Alignment bwaAlnn = (*it);
                        if (bwaAlnn.matches + bwaAlnn.inDels + bwaAlnn.substitutions > bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions) // only exact matches
                        {
                            bestAlnBwa = bwaAlnn;
                        }
                        else if (bwaAlnn.matches + bwaAlnn.inDels + bwaAlnn.substitutions == bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions)
                        {
                            if (bwaAlnn.inDels + bwaAlnn.substitutions < bestAlnBwa.inDels + bestAlnBwa.substitutions) // InDel has the same wight as substitution
                            {
                                bestAlnBwa = bwaAlnn;
                            }
                        }
                    }
                    //check if the best alignment id is the same
                    bool foundSameRead = false;
                    if(bestAlnBwa.readID == bestAlnPm.readID) {
                        foundSameRead = true;
                    }
                    if (bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                    {
                        outputFile << "Bwa ouperformed: larger aln_length for [" << bestAlnBwa.queryID << ", " << bestAlnBwa.readID << 
                        "]: Bwa aln_len: " <<   bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions << 
                        " - parmik aln_len: " <<  bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << endl;
                        bwaOutperform++;
                        if (foundSameRead){
                            outputFile << "sameReadBwaOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            sameReadBwaOutperform++;
                        } else {
                            outputFile << "differentReadBwaOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            differentReadBwaOutperform++;
                        }
                    } else if (bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions < bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                    {
                        if (foundSameRead){
                            outputFile << "sameReadPmOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            sameReadPmOutperform++;
                        } else {
                            outputFile << "differentReadPmOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            differentReadPmOutperform++;
                        }
                        parmikOutperform++;
                    } else if (bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions)
                    {
                        if (bestAlnBwa.inDels + bestAlnBwa.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) // InDel has the same wight as substitution
                        {
                            outputFile << "bwa ouperformed: smaller edits for [" << bestAlnBwa.queryID << ", " << bestAlnBwa.readID << 
                            "]: Bwa edit_len: " <<   bestAlnBwa.inDels + bestAlnBwa.substitutions << 
                            " - parmik edit_len: " <<  bestAlnPm.inDels + bestAlnPm.substitutions << endl;
                            bwaOutperform++;
                            if (foundSameRead){
                                sameReadBwaOutperform++;
                                outputFile << "sameReadBwaOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            } else {
                                differentReadBwaOutperform++;
                                outputFile << "differentReadBwaOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            }
                        } else if (bestAlnBwa.inDels + bestAlnBwa.substitutions > bestAlnPm.inDels + bestAlnPm.substitutions)
                        {
                            if (foundSameRead){
                            outputFile << "sameReadPmOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            sameReadPmOutperform++;
                            } else {
                                outputFile << "differentReadPmOutperform, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                                differentReadPmOutperform++;
                            }
                            parmikOutperform++;
                        } else{
                            if (foundSameRead){
                                sameReadEqual++;
                                outputFile << "sameReadEqual, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            } else {
                                differentReadEqual++;
                                outputFile << "differentReadEqual, bwa cigar: " << bestAlnBwa.cigar << " - parmik cigar: " << bestAlnPm.cigar << endl;
                            }
                            equalPerformance++;
                        }
                    }
                }
                // check whether the bwa Low PI FN are larger than its best alignment
                for (auto it = bwaLowPI_Alignments.begin(); it!= bwaLowPI_Alignments.end(); it++) {
                    Alignment bwaLowPIAln = (*it);
                    if (bwaLowPIAln.matches + bwaLowPIAln.inDels + bwaLowPIAln.substitutions >= bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                    {
                        outputFile << "bwaLowPIAln is longer than bestAlnPm [" << bwaLowPIAln.cigar << ", " << bestAlnPm.cigar << endl;
                        bwaLowPIOutperformBestPm++;
                    }
                    else if (bwaLowPIAln.matches + bwaLowPIAln.inDels + bwaLowPIAln.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions)
                    {
                        if (bwaLowPIAln.inDels + bwaLowPIAln.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) // InDel has the same wight as substitution
                        {
                            outputFile << "bwaLowPIAln is longer than bestAlnPm [" << bwaLowPIAln.cigar << ", " << bestAlnPm.cigar << endl;
                            bwaLowPIOutperformBestPm++;
                        }
                    }
                }
            }
        }
        outputFile << "Total Bwa TN : " << totalTN << endl;
        outputFile << "Total Bwa FN : " << bwaFN << endl;
        outputFile << "Total Parmik FNN : " << parmikFN << endl;
        outputFile << "Total Bwa Outperform : " << bwaOutperform << endl;
        outputFile << "Total Parmik Outperform : " << parmikOutperform << endl;
        outputFile << "Total Bwa Equal : " << equalPerformance << endl;
        outputFile << "Same Read Bwa Outperform : " << sameReadBwaOutperform << endl;
        outputFile << "Different Read Bwa Outperform : " << differentReadBwaOutperform << endl;
        outputFile << "Same Read Parmik Outperform : " << sameReadPmOutperform << endl;
        outputFile << "Different Read Parmik Outperform : " << differentReadPmOutperform << endl;
        outputFile << "Same Read Equal : " << sameReadEqual << endl;
        outputFile << "Different Read Equal : " << differentReadEqual << endl;
        outputFile << "bwaLowPI_FN : " << bwaLowPI_FN << endl;
        outputFile << "bwaLowPI_OutperformBestPm : " << bwaLowPIOutperformBestPm << endl;
        outputFile.close();
    }
    
};

#endif