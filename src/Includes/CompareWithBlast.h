#ifndef COMAPREWITHBLAST_H
#define COMAPREWITHBLAST_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include "BlastReader.h"
#include "Aligner.h"
#include "PostFilter.h"
#include "Alignment.h"
#include <vector>

#define REPORT_BEST_ALN 0
#define REPORT_PARMIK_FN 0
#define REPORT_ALN_PER_Q 0

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

    void comparePmWithBlast(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, const uint32_t queryCount)
    {
        //read the blast file
        BlastReader blastReader(cfg.otherToolOutputFileAddress);
        vector<Alignment> blastAlignments, parmikAlignments;
        blastReader.parseFile(queryCount, blastAlignments);
        //read the parmik file
        SamReader parmikSam(cfg.baselineBaseAddress);
        parmikSam.parseFile(queryCount, parmikAlignments, false);

        uint32_t blastFN = 0, parmikFN = 0, totalTN = 0, numberOfQueryContainN = 0;
        uint64_t blastOutperform = 0, parmikOutperform = 0, equalPerformance = 0;
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
                    query_blastAlignments.push_back(aln);
                }
            }
            size_t blastReadPerQuery = query_blastAlignments.size();

            vector<Alignment> query_parmikAlignments;
            for (const Alignment& aln : parmikAlignments) 
            {
                if ((uint32_t)aln.queryID == queryInd) {
                    query_parmikAlignments.push_back(aln);
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
                for (auto it = query_blastAlignments.begin(); it != query_blastAlignments.end(); it++) {
                    Alignment blastAlnn = (*it);
                    if (blastAlnn.matches + blastAlnn.inDels + blastAlnn.substitutions > bestAlnBlast.matches + bestAlnBlast.inDels + bestAlnBlast.substitutions) // only exact matches
                    {
                        bestAlnBlast = blastAlnn;
                    }
                    else if (blastAlnn.matches + blastAlnn.inDels + blastAlnn.substitutions == bestAlnBlast.matches + bestAlnBlast.inDels + bestAlnBlast.substitutions)
                    {
                        if (blastAlnn.inDels + blastAlnn.substitutions < bestAlnBlast.inDels + bestAlnBlast.substitutions) // InDel has the same wight as substitution
                        {
                            bestAlnBlast = blastAlnn;
                        }
                    }
                }
                if (bestAlnBlast.matches + bestAlnBlast.inDels + bestAlnBlast.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                {
                    blastOutperform++;
                } else if (bestAlnBlast.matches + bestAlnBlast.inDels + bestAlnBlast.substitutions < bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) // only exact matches
                {
                    parmikOutperform++;
                } else if (bestAlnBlast.matches + bestAlnBlast.inDels + bestAlnBlast.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions)
                {
                    if (bestAlnBlast.inDels + bestAlnBlast.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) // InDel has the same wight as substitution
                    {
                        blastOutperform++;
                    } else if (bestAlnBlast.inDels + bestAlnBlast.substitutions > bestAlnPm.inDels + bestAlnPm.substitutions)
                    {
                        parmikOutperform++;
                    } else{
                        equalPerformance++;
                    }
                }
            }
        }
        cout << "Total Blast TN : " << totalTN << endl;
        cout << "Total Blast FN : " << blastFN << endl;
        cout << "Total Parmik FNN : " << parmikFN << endl;
        cout << "Total Blast Outperform : " << blastOutperform << endl;
        cout << "Total Parmik Outperform : " << parmikOutperform << endl;
        cout << "Total Blast Equal : " << equalPerformance << endl;
    }
};

#endif