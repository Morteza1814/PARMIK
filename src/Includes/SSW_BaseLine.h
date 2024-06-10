#ifndef SSW_BASELINE_H
#define SSW_BASELINE_H

#include <omp.h>
#include "Alignment.h"

class SSW_BaseLine {
private:
    uint16_t regionSize; // min region size to be considered for the alignement
    double percentageIdentity;
public:
    SSW_BaseLine(uint16_t R, double i) : regionSize(R), percentageIdentity(i) {}

    uint32_t getMatchesCount(string cigarStr) {
        uint32_t cnt = 0;
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == '=') {
                cnt++;
            }
        }
        return cnt;
    }

    bool checkIdentityPercentange(string cigarStr){
        uint32_t matches = getMatchesCount(cigarStr);
        double identity =  (double) matches / (double) cigarStr.size();
        if(DEBUG_MODE) cout << "matches: " << matches << ", cigarStr.size : " << cigarStr.size() << ", identity: " << identity << endl;
        if(identity < percentageIdentity)
            return false;
        return true;
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

    Alignment alignDifferentPenaltyScores(string query, string read, uint32_t queryID, uint32_t readID, bool isForwardStran, vector<Penalty> &penalties)
    {
        Alignment bestAlignment;
        if (penalties.size() == 0)
        {
            Alignment aln;
            aln.read = read;
            aln.query = query;
            aln.readID = readID;
            aln.queryID = queryID;
            if (!isForwardStran) aln.flag = 16;
            smithWatermanAligner(aln, 1, 1, 1, 1);
            string cigarStr = convertCigarToStr(aln.cigar);
            cigarStr = trimClips(cigarStr);
            if(hasConsecutiveMatches(cigarStr, regionSize) && checkIdentityPercentange(cigarStr)) {
                if(DEBUG_MODE) cout << "cigar len was smaller than the R or K at the beginning\n";
                return aln;
            } else {
                //return an empty alignment
                return bestAlignment;
            }
        }
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
            smithWatermanAligner(aln, penalty.matchPenalty, penalty.mismatchPenalty, penalty.gapOpenPenalty, penalty.gapExtendPenalty);
            //check the alignment based on the criteria
            string cigarStr = convertCigarToStr(aln.cigar);
            cigarStr = trimClips(cigarStr);
            if(hasConsecutiveMatches(cigarStr, regionSize) && checkIdentityPercentange(cigarStr)) {
                if(DEBUG_MODE) cout << "cigar len was smaller than the R or K at the beginning\n";
                if (aln.partialMatchSize > bestAlignment.partialMatchSize || (aln.partialMatchSize == bestAlignment.partialMatchSize && aln.editDistance < bestAlignment.editDistance))
                    bestAlignment = aln;
            }
        }
        return bestAlignment;
    }

    void findPartiaMatches(tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, uint32_t queryCount, bool isForwardStrand, string parmikAlignments, vector<Penalty> penalties, uint32_t queryBaseIndex)
    {
        ofstream pAln;
        if (isForwardStrand) {
            pAln.open(parmikAlignments);
        } else {
            pAln.open(parmikAlignments, std::ios::app);
        }
        multiset<uint32_t> matchesPerQuery;
        cout << "Starting alignment for all queries [" << (isForwardStrand ? ("fwd"):("rev")) << "]..." << endl;
        cout << "queryBaseIndex: " << queryBaseIndex << endl;
        for (size_t i = 0; i < queryCount; i++)
        {
            // if(i % 100 == 0)
                // cout << i << " queries processed..." << endl;
            uint32_t qID = i + queryBaseIndex; 
            auto itq = queries.find(qID);
            if (itq == queries.end())
                continue;
            string query = itq->second;
            if (query.find('n') != string::npos || query.find('N') != string::npos)
                continue;
            // read the candidate reads of cheap k-mer filter of front
            tsl::robin_map <uint32_t, Alignment> alignments;
            auto start = chrono::high_resolution_clock::now();
            #pragma omp parallel for
            for (size_t j = 0; j < reads.size(); j++)
            {
                // if(j % 10000 == 0)
                //     cout << j << " reads processed in Thread: " << omp_get_thread_num() << endl;
                auto itr = reads.find(j);
                if (itr == reads.end())
                    continue;
                string read = itr->second;
                Alignment aln = alignDifferentPenaltyScores(query, read, qID, j, isForwardStrand, penalties);
                if (aln.partialMatchSize > 0){
                    #pragma omp critical
                    {
                        alignments.insert(make_pair(j, aln));
                    }
                }
            }
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            cout << "query [" << qID << "]: # of matches found " << alignments.size() << " and it took: " << static_cast<double>(duration.count()) / 1'000'000.0 << "ms" << endl;
            //dump the alignments
            for (auto it = alignments.begin(); it!= alignments.end(); it++)
            {
                dumpSam(pAln, it->second);
            }
            matchesPerQuery.insert(alignments.size());
            // cout << "queryID: " << qID << ", " << (isForwardStrand ? ("fwd"):("rev")) <<", total matches: " << alignments.size() << i << ", matches accepted: " << matchesAccepted << ", matches passed with starting pos over ED: " << matchesPassedWithStartingPosOverED << endl;
        }
        Utilities<uint32_t> util; 
        tuple<uint32_t, uint32_t, uint32_t> matchesPerQueryTuple = util.calculateStatistics(matchesPerQuery);
        printf("PARMIK's matches per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(matchesPerQueryTuple), get<1>(matchesPerQueryTuple), get<2>(matchesPerQueryTuple));
    }

    void dumpSam(ofstream &oSam, Alignment l)
    {
        // cerr << "queryID: " << l.queryID << ", readID: "<< l.readID << "readRegionStartPos: " << l.readRegionStartPos << endl;
        // assert((l.readRegionStartPos >= 0 && l.readRegionStartPos < contigSize) && "wrong readRegionStartPos");
        /*oSam << l.queryID << '\t' << l.flag << '\t' << l.readID << '\t' << l.readRegionStartPos << '\t'
                << "*" << '\t' << l.cigar << '\t' << "*" << '\t' << "*" << '\t' << "*" << '\t' 
                << l.read << '\t' << "*" << '\t' << "NM:i:" + to_string(l.substitutions) << '\n';*/
        oSam << l.queryID << '\t' << l.flag << '\t' << l.readID << '\t' << l.readRegionStartPos << '\t'
            << l.cigar << '\t' << "NM:i:" + to_string(l.substitutions) << '\n';
    }
};

#endif