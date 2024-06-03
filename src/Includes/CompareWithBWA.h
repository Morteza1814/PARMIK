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

    void comparePmWithBwa(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, const uint32_t queryCount)
    {
        vector<Alignment> bwaAlignments, parmikAlignments;
        //read the bwa file
        SamReader bwaSam(cfg.otherToolOutputFileAddress);
        bwaSam.parseFile(queryCount, bwaAlignments, false, true);
        //read the parmik file
        SamReader parmikSam(cfg.baselineBaseAddress);
        parmikSam.parseFile(queryCount, parmikAlignments, false, false);

        uint32_t bwaFN = 0, parmikFN = 0, totalTN = 0, numberOfQueryContainN = 0;
        uint64_t bwaOutperform = 0, parmikOutperform = 0, equalPerformance = 0;
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
            vector<uint32_t> bwaTPReadIds;
            if(bwaReadPerQuery == 0 && parmikReadPerQuery == 0){
                //TN
                totalTN++;
                // cmp << "TN for BWA" << endl;
            } else if(bwaReadPerQuery == 0 && parmikReadPerQuery > 0) {
                //FN
                bwaFN++;
            } else if(bwaReadPerQuery > 0 && parmikReadPerQuery == 0) {
                //FN
                parmikFN++;
            } else { // both > 0
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
                if (bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                {
                    bwaOutperform++;
                } else if (bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions < bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                {
                    parmikOutperform++;
                } else if (bestAlnBwa.matches + bestAlnBwa.inDels + bestAlnBwa.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions)
                {
                    if (bestAlnBwa.inDels + bestAlnBwa.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) // InDel has the same wight as substitution
                    {
                        bwaOutperform++;
                    } else if (bestAlnBwa.inDels + bestAlnBwa.substitutions > bestAlnPm.inDels + bestAlnPm.substitutions)
                    {
                        parmikOutperform++;
                    } else{
                        equalPerformance++;
                    }
                }
            }
        }
        cout << "Total Bwa TN : " << totalTN << endl;
        cout << "Total Bwa FN : " << bwaFN << endl;
        cout << "Total Parmik FNN : " << parmikFN << endl;
        cout << "Total Bwa Outperform : " << bwaOutperform << endl;
        cout << "Total Parmik Outperform : " << parmikOutperform << endl;
        cout << "Total Bwa Equal : " << equalPerformance << endl;
    }
    
};

#endif