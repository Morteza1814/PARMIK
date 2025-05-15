#ifndef COMAPREWITHBASELINE_H
#define COMAPREWITHBASELINE_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include "BlastReader.h"
#include "Aligner.h"
#include "Alignment.h"
#include "SamReader.h"
#include "Container.h"

#define REPORT_BEST_ALN 0
#define REPORT_baseLine_FN 0
#define REPORT_ALN_PER_Q 0
#define CHECK_EXACT_MATCH_CRITERION 1

using namespace std;


class CompareWithBaseLine {
private:
    double percentageIdentity;
    //Alignment sizes
    //FN
    map<uint16_t, uint64_t> tool2FN_aln_sz;
    map<uint16_t, uint64_t> baselineFN_aln_sz;
    //TP
    map<uint16_t, uint64_t> tool2TP_aln_sz;
    map<uint16_t, uint64_t> baselineTP_aln_sz;
    //best aln
    map<uint16_t, uint64_t> tool2Best_aln_sz;
    map<uint16_t, uint64_t> baselineBest_aln_sz;
    
    //differences
    //TP
    map<uint16_t, uint64_t> tool2TPOutperform_aln_sz_difference;
    map<uint16_t, uint64_t> baselineTPOutperform_aln_sz_difference;
    //equal aln_sz but diferent edits
    map<uint16_t, uint64_t> tool2TPOutperform_edit_difference;
    map<uint16_t, uint64_t> baselineTPOutperform_edit_difference;
    //best aln
    map<uint16_t, uint64_t> tool2BestOutperform_aln_sz_difference;
    map<uint16_t, uint64_t> baselineBestOutperform_aln_sz_difference;
    //equal aln_sz but diferent edits
    map<uint16_t, uint64_t> tool2BestOutperform_edit_difference;
    map<uint16_t, uint64_t> baselineBestOutperform_edit_difference;
public: 

    CompareWithBaseLine(double pi) : percentageIdentity(pi) {}

    void addTool2TPOutperform(uint16_t difference) {
        if (tool2TPOutperform_aln_sz_difference.find(difference) == tool2TPOutperform_aln_sz_difference.end()) {
            tool2TPOutperform_aln_sz_difference[difference] = 1;
        } else {
            tool2TPOutperform_aln_sz_difference[difference]++;
        }
    }

    void addBaselineTPOutperform(uint16_t difference) {
        if (baselineTPOutperform_aln_sz_difference.find(difference) == baselineTPOutperform_aln_sz_difference.end()) {
            baselineTPOutperform_aln_sz_difference[difference] = 1;
        } else {
            baselineTPOutperform_aln_sz_difference[difference]++;
        }
    }

    void addTool2BestOutperform(uint16_t difference) {
        if (tool2BestOutperform_aln_sz_difference.find(difference) == tool2BestOutperform_aln_sz_difference.end()) {
            tool2BestOutperform_aln_sz_difference[difference] = 1;
        } else {
            tool2BestOutperform_aln_sz_difference[difference]++;
        }
    }

    void addBaselineBestOutperform(uint16_t difference) {
        if (baselineBestOutperform_aln_sz_difference.find(difference) == baselineBestOutperform_aln_sz_difference.end()) {
            baselineBestOutperform_aln_sz_difference[difference] = 1;
        } else {
            baselineBestOutperform_aln_sz_difference[difference]++;
        }
    }

    void addTool2TPOutperform_edit(uint16_t difference) {
        if (tool2TPOutperform_edit_difference.find(difference) == tool2TPOutperform_edit_difference.end()) {
            tool2TPOutperform_edit_difference[difference] = 1;
        } else {
            tool2TPOutperform_edit_difference[difference]++;
        }
    }

    void addBaselineTPOutperform_edit(uint16_t difference) {
        if (baselineTPOutperform_edit_difference.find(difference) == baselineTPOutperform_edit_difference.end()) {
            baselineTPOutperform_edit_difference[difference] = 1;
        } else {
            baselineTPOutperform_edit_difference[difference]++;
        }
    }

    void addTool2BestOutperform_edit(uint16_t difference) {
        if (tool2BestOutperform_edit_difference.find(difference) == tool2BestOutperform_edit_difference.end()) {
            tool2BestOutperform_edit_difference[difference] = 1;
        } else {
            tool2BestOutperform_edit_difference[difference]++;
        }
    }

    void addBaselineBestOutperform_edit(uint16_t difference) {
        if (baselineBestOutperform_edit_difference.find(difference) == baselineBestOutperform_edit_difference.end()) {
            baselineBestOutperform_edit_difference[difference] = 1;
        } else {
            baselineBestOutperform_edit_difference[difference]++;
        }
    }

    void addTool2TPAlnSz(uint16_t aln_sz) {
        if (tool2TP_aln_sz.find(aln_sz) == tool2TP_aln_sz.end()) {
            tool2TP_aln_sz[aln_sz] = 1;
        } else {
            tool2TP_aln_sz[aln_sz]++;
        }
    }

    void addBaselineTPAlnSz(uint16_t aln_sz) {
        if (baselineTP_aln_sz.find(aln_sz) == baselineTP_aln_sz.end()) {
            baselineTP_aln_sz[aln_sz] = 1;
        } else {
            baselineTP_aln_sz[aln_sz]++;
        }
    }

    void removeBaselineTPAlnSz(uint16_t aln_sz) {
        if (baselineTP_aln_sz.find(aln_sz) == baselineTP_aln_sz.end()) {
            //do nothing
        } else {
            if(baselineTP_aln_sz[aln_sz] > 0) baselineTP_aln_sz[aln_sz]--;
        }
    }

    void addTool2BestAlnSz(uint16_t aln_sz) {
        if (tool2Best_aln_sz.find(aln_sz) == tool2Best_aln_sz.end()) {
            tool2Best_aln_sz[aln_sz] = 1;
        } else {
            tool2Best_aln_sz[aln_sz]++;
        }
    }

    void addBaselineBestAlnSz(uint16_t aln_sz) {
        if (baselineBest_aln_sz.find(aln_sz) == baselineBest_aln_sz.end()) {
            baselineBest_aln_sz[aln_sz] = 1;
        } else {
            baselineBest_aln_sz[aln_sz]++;
        }
    }

    void addTool2FNAlnSz(uint16_t aln_sz) {
        if (tool2FN_aln_sz.find(aln_sz) == tool2FN_aln_sz.end()) {
            tool2FN_aln_sz[aln_sz] = 1;
        } else {
            tool2FN_aln_sz[aln_sz]++;
        }
    }

    void addBaselineFNAlnSz(uint16_t aln_sz) {
        if (baselineFN_aln_sz.find(aln_sz) == baselineFN_aln_sz.end()) {
            baselineFN_aln_sz[aln_sz] = 1;
        } else {
            baselineFN_aln_sz[aln_sz]++;
        }
    }

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

    string blastGetCigarStr(uint32_t queryS, string queryAligned, string readAligned, uint32_t contigSize, bool skipClips = false)
    {
        stringstream cigar;
        if(!skipClips) {
            for(uint32_t i = 0; i < queryS; i++) {
                cigar << 'S';
            }
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
        if(remainedClips > 0 && !skipClips){
            for(int i = 0; i < remainedClips; i++) {
                cigar << 'S';
            }
        }
        return cigar.str();
    }

    void readChunkOfBaseLine(vector<Alignmen_Baseline> &baseLineSamAlignments, const string alignmentFileAddressPrefix, uint32_t queryRangeBase, uint32_t queryRangeSize)
    {
        baseLineSamAlignments.clear();
        string alignmentFileAddress = alignmentFileAddressPrefix + "_" + to_string(queryRangeBase) + "-" + to_string(queryRangeBase + queryRangeSize - 1) + ".txt";
        // cout << "alignmentFileAddress: " << alignmentFileAddress << endl;
        uint32_t queryCount = queryRangeBase + queryRangeSize - 1;
        SamReader baseLineSam(alignmentFileAddress);
        baseLineSam.parseFileBaseline(queryCount, baseLineSamAlignments, true);
        cout << "Sam: Read " << baseLineSamAlignments.size() << " alignments." << endl;
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

    string getCigarStr(Alignment aln, Config cfg, string toolName){
        if(toolName == "BLAST" || toolName == "blast") {
            return blastGetCigarStr(aln.queryRegionStartPos, aln.alignedQuery, aln.alignedRead, cfg.contigSize, true);
        } else {
            return convertCigarToStr(aln.cigar, true);
        }
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

    bool checkIdentityPercentangeWithAln(Alignment aln) {
        double identity =  (double) aln.matches / (double) aln.partialMatchSize;
        if(identity < percentageIdentity)
            return false;
        return true;
    }

    int longestRunOfMatches(const std::string& cigarStr) {
        int maxRun = 0;
        int currentRun = 0;

        for (char c : cigarStr) {
            if (c == '=') {
                currentRun++;
                maxRun = std::max(maxRun, currentRun);
            } else {
                currentRun = 0; // Reset the current run if the character is not '='
            }
        }

        return maxRun;
    }

    //baseLine is always baseline
    void compareWithBaseLine(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress, 
    const uint32_t queryCount, const string alnReportAddressBase, string baseLineFilePrefixAddress)
    {
        const string tool2Name = cfg.otherTool;
        cout << "Comparing Baseline and " + tool2Name << endl;
        ofstream cmp(comparisonResultsFileAddress);
        cout << "cmp file address: " << comparisonResultsFileAddress << endl;
        if(!cmp.is_open())
            cout << "cmp file not opened" << endl;
        ofstream alnPerQ;
        ofstream baseLineFn;
        ofstream fnAlnSz;
        ofstream tpAlnCmp, bestAlnCmp, tpEdCmp, bestEdCmp;
        ofstream tpAlnSz, bestAlnSz;
        // if(REPORT_ALN_PER_Q) alnPerQ.open(alnPerQueryFileAddress);
        // if(REPORT_baseLine_FN) baseLineFn.open(baseLineFnReadsFileAddress);
        // if(REPORT_BEST_ALN) bestAlnCmp.open(bestAlnCmpFileAddress);
        bestAlnCmp.open(alnReportAddressBase + "CmpAlnSz_Best.txt");
        tpAlnCmp.open(alnReportAddressBase + "CmpAlnSz_tp.txt");
        bestEdCmp.open(alnReportAddressBase + "CmpEd_Best.txt");
        tpEdCmp.open(alnReportAddressBase + "CmpEd_tp.txt");
        tpAlnSz.open(alnReportAddressBase + "AlnSz_tp.txt");
        bestAlnSz.open(alnReportAddressBase + "AlnSz_Best.txt");
        fnAlnSz.open(alnReportAddressBase + "AlnSz_fn.txt");
        //read the alignment files
        vector<Alignment> tool2Alignments;
        bool isBwotie = (tool2Name == "bowtie" || tool2Name == "BOWTIE") ? true:false;
        bool isParmik = (tool2Name == "parmik" || tool2Name == "PARMIK") ? true:false;

        if(tool2Name == "BLAST" || tool2Name == "blast") {
            BlastReader blastReader(cfg.otherToolOutputFileAddress);
            blastReader.parseFile(queryCount, tool2Alignments);
        } else if(isParmik || isBwotie) {
            SamReader parmikSam(cfg.otherToolOutputFileAddress);
            parmikSam.parseFile(queryCount, tool2Alignments, false, isBwotie);//bowtie's output should be parsed like bwa
        }
        IndexContainer<uint32_t, Alignment> queryTool2Alignment;
        for (const Alignment& aln : tool2Alignments)
        {
            queryTool2Alignment.put(aln.queryID, aln);
        }
        tool2Alignments.clear();


        uint32_t numberOfQueryContainN = 0,
        queriestool2FoundMatch = 0, queriesblfoundMatch = 0,
        baseLineTotalNumberOfReadIDs = 0, tool2TotalNumberOfReadIDs = 0,
        numnberOfBaseLine_FN = 0, tool2FN = 0, tool2FN_noCriteriaFittedMatches = 0, // places where tool2 found no matches or none of the matches conform with the criteria
        totalTool2TN = 0,
        totalbaseLineTN = 0,
        totaltool2TP = 0, 
        totaltool2FP = 0,
        tool2TP_nobaseLineMatches_allQ = 0, tool2TP_tool2Outperfomed_allQ = 0, tool2TP_baseLineOutperfomed_allQ = 0, tool2TP_tool2EqualbaseLine_allQ = 0,
        baseLineTP_notool2Matches_allQ = 0,
        tool2FP_lowAlnLen_allQ = 0, tool2FP_lowPercentageIdentity_allQ = 0,
        tool2Best_tool2Outperfomed_allQ = 0, tool2Best_baseLineOutperfomed_allQ = 0, tool2Best_tool2EqualbaseLine_allQ = 0;
        //all query level parameters sets (if total is needed, just add all the elements in the set)
        multiset<uint32_t> baseLineReadPerQuerySet, tool2ReadPerQuerySet;
        multiset<uint32_t> tool2TPPerQuery, tool2FPPerQuery;
        multiset<uint32_t> tool2TP_nobaseLineMatches_alnLen_allQ, tool2TP_tool2Outperfomed_alnLen_allQ, tool2TP_baseLineOutperfomed_alnLen_allQ, tool2TP_tool2EqualbaseLine_alnLen_allQ; 
        multiset<uint32_t> tool2TP_nobaseLineMatches_ed_allQ, tool2TP_tool2Outperfomed_ed_allQ, tool2TP_baseLineOutperfomed_ed_allQ, tool2TP_tool2EqualbaseLine_ed_allQ; 
        multiset<uint32_t> baseLineTP_notool2Matches_alnLen_allQ, baseLineTP_notool2Matches_ed_allQ;
        multiset<uint32_t> baseLineTP_tool2Outperfomed_alnLen_allQ, baseLineTP_baseLineOutperfomed_alnLen_allQ, baseLineTP_tool2EqualbaseLine_alnLen_allQ; 
        multiset<uint32_t> baseLineTP_tool2Outperfomed_ed_allQ, baseLineTP_baseLineOutperfomed_ed_allQ, baseLineTP_tool2EqualbaseLine_ed_allQ; 
        // best alignment comparisons
        multiset<uint32_t> tool2Best_tool2Outperfomed_alnLen_allQ, tool2Best_baseLineOutperfomed_alnLen_allQ, tool2Best_tool2EqualbaseLine_alnLen_allQ; 
        multiset<uint32_t> tool2Best_tool2Outperfomed_ed_allQ, tool2Best_baseLineOutperfomed_ed_allQ, tool2Best_tool2EqualbaseLine_ed_allQ; 
        multiset<uint32_t> baseLineBest_tool2Outperfomed_alnLen_allQ, baseLineBest_baseLineOutperfomed_alnLen_allQ, baseLineBest_tool2EqualbaseLine_alnLen_allQ; 
        multiset<uint32_t> baseLineBest_tool2Outperfomed_ed_allQ, baseLineBest_baseLineOutperfomed_ed_allQ, baseLineBest_tool2EqualbaseLine_ed_allQ; 
        multiset<uint32_t> tool2FP_lowAlnLen_alnLen_allQ, tool2FP_lowAlnLen_matchsize_allQ, tool2FP_lowerExactMatchKmers_alnLen_allQ, tool2FP_lowPercentageIdentity_alnLen_allQ;
        multiset<uint32_t> tool2FP_lowAlnLen_ed_allQ, tool2FP_lowerExactMatchKmers_ed_allQ, tool2FP_lowPercentageIdentity_ed_allQ;
        multiset<uint32_t> tool2FN_baseLine_alnLen_allQ, tool2FN_baseLine_ed_allQ;// baseLine alignment characteristics for all queries when tool2 did not find a match
        multiset<uint32_t> tool2FN_noCriteria_baseLine_alnLen_allQ, tool2FN_noCriteria_tool1_ed_allQ;
        multiset<uint32_t> tool2FN_noCriteria_baseLine_ed_allQ, baseLineFN_tool2_alnLen_allQ, baseLineFN_tool2_ed_allQ;// tool2 alignment characteristics for all queries when baseLine did not find a match
        vector<Alignmen_Baseline> baseLineSamAlignments;
        IndexContainer<uint32_t, Alignmen_Baseline> queryBaselineAlignment;
        for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            uint32_t baseLineBestAlnSize = 0;
            string query = queries[queryInd];
            if(query.find('N') != string::npos || query.find('n') != string::npos)
            {
                cmp << "query contains N!" << endl;
                numberOfQueryContainN++; 
                if(REPORT_ALN_PER_Q) alnPerQ << queryInd << " 0 0"<< endl;
                continue;
            }
            if(queryInd % 1000 == 0) {
                queryBaselineAlignment.clear();
                cout << "queries processed: " << queryInd << " / " << queryCount << endl;
                readChunkOfBaseLine(baseLineSamAlignments, baseLineFilePrefixAddress, queryInd, 1000);
                for (const Alignmen_Baseline& aln : baseLineSamAlignments)
                {
                    queryBaselineAlignment.put(aln.queryID, aln);
                }
                baseLineSamAlignments.clear();
            }
            //baseline alignments for query
            tsl::robin_map<uint32_t, Alignmen_Baseline> baselineReadID_Aln;
            Alignmen_Baseline baseLineBestAln;

            auto blQueryRrange = queryBaselineAlignment.getRange(queryInd);
            for (auto it = blQueryRrange.first; it != blQueryRrange.second; it++){
                Alignmen_Baseline aln = it->second;
                string cigarStr = convertCigarToStr(aln.cigar, true);
                // only add alignments whose pi > PI
                if (checkIdentityPercentange(cigarStr)) {
                    auto readIt = baselineReadID_Aln.find(aln.readID);
                    if(readIt == baselineReadID_Aln.end()) {
                        baselineReadID_Aln[aln.readID] = aln;
                        addBaselineTPAlnSz(aln.matches + aln.inDels + aln.substitutions);
                    } else {
                        auto readAln = readIt->second;
                        if ((readAln.matches + readAln.inDels + readAln.substitutions < aln.matches + aln.inDels + aln.substitutions) || 
                        ((readAln.matches + readAln.inDels + readAln.substitutions == aln.matches + aln.inDels + aln.substitutions) &&
                        (readAln.inDels + readAln.substitutions > aln.inDels + aln.substitutions))) {
                            baselineReadID_Aln[aln.readID] = aln;
                            removeBaselineTPAlnSz(readAln.matches + readAln.inDels + readAln.substitutions);
                            addBaselineTPAlnSz(aln.matches + aln.inDels + aln.substitutions);
                        }
                    }
                    if (aln.matches + aln.inDels + aln.substitutions > baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions) // only exact matches
                    {
                        baseLineBestAln = aln;
                    }
                    else if (aln.matches + aln.inDels + aln.substitutions == baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions)
                    {
                        if (aln.inDels + aln.substitutions < baseLineBestAln.inDels + baseLineBestAln.substitutions) // InDel has the same wight as substitution
                        {
                            baseLineBestAln = aln;
                        }
                    }
                }   
            }

            size_t baseLineReadPerQuery = baselineReadID_Aln.size();
            baseLineTotalNumberOfReadIDs += baseLineReadPerQuery;
            baseLineReadPerQuerySet.insert(baseLineReadPerQuery);
            if(baseLineReadPerQuery > 0){
                queriesblfoundMatch++;
                baseLineBestAlnSize = baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions;
                addBaselineBestAlnSz(baseLineBestAlnSize);
            }
            //query level parameters
            uint32_t tool2TP = 0, tool2FP = 0;
            uint32_t tool2TP_nobaseLineMatches = 0, tool2TP_tool2Outperfomed = 0, tool2TP_baseLineOutperfomed = 0, tool2TP_tool2EqualbaseLine = 0;
            uint32_t baseLineTP_notool2Matches = 0;
            uint32_t tool2Best_tool2Outperfomed = 0, tool2Best_baseLineOutperfomed = 0, tool2Best_tool2EqualbaseLine = 0;
            uint32_t tool2FP_lowAlnLen = 0, tool2FP_lowPercentageIdentity = 0;
            cmp << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            cmp << "Q : " << query << ", queryInd: " << queryInd << endl;
            Alignment bestAlntool2;

            tsl::robin_map<uint32_t, Alignment> query_t2Alignments;
            auto tool2QueryRrange = queryTool2Alignment.getRange(queryInd);
            for (auto it = tool2QueryRrange.first; it != tool2QueryRrange.second; it++){
                Alignment aln = it->second;
                auto readIt = query_t2Alignments.find(aln.readID);
                if(readIt == query_t2Alignments.end()) {
                    query_t2Alignments[aln.readID] = aln;
                } else {
                    Alignment readAln = readIt->second;
                    if ((readAln.matches + readAln.inDels + readAln.substitutions < aln.matches + aln.inDels + aln.substitutions) || 
                    ((readAln.matches + readAln.inDels + readAln.substitutions == aln.matches + aln.inDels + aln.substitutions) &&
                    (readAln.inDels + readAln.substitutions > aln.inDels + aln.substitutions))) {
                        query_t2Alignments[aln.readID] = aln;
                    }
                }
            }
            size_t tool2ReadPerQuery = query_t2Alignments.size();
            tool2TotalNumberOfReadIDs += tool2ReadPerQuery;
            tool2ReadPerQuerySet.insert(tool2ReadPerQuery);
            cmp << "# of reads found by " + tool2Name + " : " << tool2ReadPerQuery << endl;
            cmp << "# of reads found by Baseline : " << baseLineReadPerQuery << endl;
            
            vector<uint32_t> tool2TPReadIds;
            if(tool2ReadPerQuery == 0 && baseLineReadPerQuery == 0){
                //TN
                totalTool2TN++;
                cmp << "TN for " + tool2Name << endl;
            } else if(tool2ReadPerQuery == 0 && baseLineReadPerQuery > 0) {
                //tool2 FN
                tool2FN++;
                tool2FN_baseLine_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                tool2FN_baseLine_ed_allQ.insert(baseLineBestAln.inDels + baseLineBestAln.substitutions);
                cmp << "FN for " + tool2Name + ", Baseline alnlen : " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", ed : " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", readID: " <<  baseLineBestAln.readID << endl;
            } else {    // if both tool2 and baseLine have alignments or only tool2 has alignments
                for (auto it = query_t2Alignments.begin(); it != query_t2Alignments.end(); it++) 
                {
                    Alignment aln = it->second;
                    // string tool2R = reads[aln.readID]; 
                    // if(tool2R.find('N') != string::npos || tool2R.find('n') != string::npos)
                    // {
                    //     numberOfTool2ReadsContainingN++;
                    //     continue;
                    // }
                    bool isFP = false;
                    string cigarStr = getCigarStr(aln, cfg, tool2Name);
                    if(CHECK_EXACT_MATCH_CRITERION && !hasConsecutiveMatches(cigarStr, cfg.kmerLength)) {
                        // isFP = true; // this is FP in terms of not finding exact match, but in total it is TP
                        aln.criteriaCode = 0x10;
                        tool2FP_lowAlnLen++;
                        uint16_t ltm = longestRunOfMatches(cigarStr);
                        tool2FP_lowAlnLen_alnLen_allQ.insert(aln.partialMatchSize);
                        tool2FP_lowAlnLen_ed_allQ.insert(aln.substitutions + aln.inDels);
                        tool2FP_lowAlnLen_matchsize_allQ.insert(ltm);
                        cmp << "FP for " + tool2Name + ", due to low aln len, alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", match size: " << ltm << endl;
                    } else if ((isBwotie && !checkIdentityPercentangeWithAln(aln)) || !checkIdentityPercentange(cigarStr)){
                        isFP = true;
                        aln.criteriaCode = 0x80;
                        tool2FP_lowPercentageIdentity++;
                        tool2FP_lowPercentageIdentity_alnLen_allQ.insert(aln.partialMatchSize);
                        tool2FP_lowPercentageIdentity_ed_allQ.insert(aln.substitutions + aln.inDels);
                        cmp << "FP for " + tool2Name + ", due to low PercentageIdentity, alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << endl;
                    }
        
                    if (isFP){
                        tool2FP++;
                    }
                    else {
                        tool2TP++;
                        tool2TPReadIds.push_back(aln.readID);
                        addTool2TPAlnSz(aln.partialMatchSize);
                        bool blfound = false;
                        auto baselineReadIt = baselineReadID_Aln.find(aln.readID);
                        Alignmen_Baseline baselineAln;
                        if(baselineReadIt != baselineReadID_Aln.end())
                        {
                            baselineAln = baselineReadIt->second;
                            blfound = true;
                        }
                        if (!blfound)
                        {
                            tool2TP_nobaseLineMatches++;
                            tool2TP_nobaseLineMatches_alnLen_allQ.insert(aln.partialMatchSize);
                            tool2TP_nobaseLineMatches_ed_allQ.insert(aln.substitutions + aln.inDels);
                            cmp << "TP for " + tool2Name + ", no baseLine match, alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << endl;
                            addBaselineFNAlnSz(bestAlntool2.partialMatchSize);
                            if(REPORT_baseLine_FN) baseLineFn << queryInd << "\t" << aln.readID << "\t" << aln.flag << endl;
                        } else {
                            if (baselineAln.matches + baselineAln.inDels + baselineAln.substitutions > aln.partialMatchSize){ //baseLine outperformed
                                tool2TP_baseLineOutperfomed++;
                                tool2TP_baseLineOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                tool2TP_baseLineOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                baseLineTP_baseLineOutperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                baseLineTP_baseLineOutperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                                addBaselineTPOutperform(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions - aln.partialMatchSize);
                            } else if (baselineAln.matches + baselineAln.inDels + baselineAln.substitutions < aln.partialMatchSize){//tool2 outperformed
                                tool2TP_tool2Outperfomed++;
                                tool2TP_tool2Outperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                tool2TP_tool2Outperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                baseLineTP_tool2Outperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                baseLineTP_tool2Outperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                                cmp << "TP for " + tool2Name + ", " + tool2Name + " outperformed in terms of alnlen, " + tool2Name + " alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", baseLine alnlen: " << baselineAln.matches + baselineAln.inDels + baselineAln.substitutions << ", ed: " << baselineAln.substitutions + baselineAln.inDels << ", readID: " << baselineAln.readID << endl;
                                addTool2TPOutperform(aln.partialMatchSize - baselineAln.matches - baselineAln.inDels - baselineAln.substitutions);
                            } else {
                                if (baselineAln.substitutions + baselineAln.inDels < aln.substitutions + aln.inDels)//baseLine outperformed
                                {
                                    tool2TP_baseLineOutperfomed++;
                                    tool2TP_baseLineOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                    tool2TP_baseLineOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                    baseLineTP_baseLineOutperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                    baseLineTP_baseLineOutperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                                    cmp << "TP for " + tool2Name + ", baseLine outperformed in terms of ed, " + tool2Name + " alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", baseLine alnlen: " << baselineAln.matches + baselineAln.inDels << ", ed: " << baselineAln.substitutions + baselineAln.inDels << ", readID: " << baselineAln.readID << endl;
                                    addBaselineTPOutperform_edit(baselineAln.substitutions + baselineAln.inDels - aln.substitutions - aln.inDels);
                                } else if (baselineAln.substitutions + baselineAln.inDels > aln.substitutions + aln.inDels)//tool2 outperformed
                                {
                                    tool2TP_tool2Outperfomed++;
                                    tool2TP_tool2Outperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                    tool2TP_tool2Outperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                    baseLineTP_tool2Outperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                    baseLineTP_tool2Outperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                                    cmp << "TP for " + tool2Name + ", " + tool2Name + " outperformed in terms of ed, " + tool2Name + " alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", baseLine alnlen: " << baselineAln.matches + baselineAln.inDels << ", ed: " << baselineAln.substitutions + baselineAln.inDels << ", readID: " << baselineAln.readID << endl;
                                    addTool2TPOutperform_edit(aln.substitutions + aln.inDels - baselineAln.substitutions - baselineAln.inDels);
                                } // tool2 and baseLine performed equally
                                else{
                                    tool2TP_tool2EqualbaseLine++;
                                    tool2TP_tool2EqualbaseLine_alnLen_allQ.insert(aln.partialMatchSize);
                                    tool2TP_tool2EqualbaseLine_ed_allQ.insert(aln.substitutions + aln.inDels);
                                    baseLineTP_tool2EqualbaseLine_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                    baseLineTP_tool2EqualbaseLine_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                                }
                            }
                        }
                        //find the best alignment for tool2
                        if (aln.partialMatchSize  > bestAlntool2.partialMatchSize)
                        {
                            bestAlntool2 = aln;
                        } else if (aln.partialMatchSize == bestAlntool2.partialMatchSize)
                        {
                            if (aln.substitutions + aln.inDels < bestAlntool2.substitutions + bestAlntool2.inDels)
                            {
                                bestAlntool2 = aln;
                            }
                        }
                    }
                }
                addTool2BestAlnSz(bestAlntool2.partialMatchSize);
                for (auto itt = baselineReadID_Aln.begin(); itt != baselineReadID_Aln.end(); itt++) 
                {
                    Alignmen_Baseline aaln = itt->second;
                    bool tool2found = false;
                    for (const auto& element : tool2TPReadIds)
                    {
                        if(element == (uint32_t)aaln.readID)
                        {
                            tool2found = true;
                            break;
                        }
                    }
                    if(!tool2found)
                    {
                        baseLineTP_notool2Matches++;
                        baseLineTP_notool2Matches_alnLen_allQ.insert(aaln.matches + aaln.inDels + aaln.substitutions);
                        baseLineTP_notool2Matches_ed_allQ.insert(aaln.substitutions + aaln.inDels);
                        addTool2FNAlnSz(aaln.matches + aaln.inDels + aaln.substitutions);
                    }
                }
                baseLineTP_notool2Matches_allQ += baseLineTP_notool2Matches;
                cmp << "baseLine TP that " + tool2Name + " did not find: " << baseLineTP_notool2Matches << endl;
                cmp << tool2Name + " TP: " << tool2TP << "\t" << tool2Name + " FP: " << tool2FP << endl;
                tool2TPPerQuery.insert(tool2TP);
                totaltool2TP += tool2TP;
                tool2FPPerQuery.insert(tool2FP);
                totaltool2FP += tool2FP;
                if(tool2TP > 0){ 
                    queriestool2FoundMatch++;
                }
                if(baseLineReadPerQuery > 0 && tool2ReadPerQuery > 0 && tool2TP == 0){//baseLine TP
                    tool2FN_noCriteriaFittedMatches++;
                    tool2FN_noCriteria_baseLine_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                    tool2FN_noCriteria_baseLine_ed_allQ.insert(baseLineBestAln.inDels + baseLineBestAln.substitutions);
                    cmp << tool2Name + "FN_noCriteriaFitted and baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine readID: " << baseLineBestAln.readID << endl;
                }
                if(baseLineReadPerQuery == 0 && tool2TP > 0){//baseLine FN
                    baseLineFN_tool2_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                    baseLineFN_tool2_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                    cmp << "baseLine FN, and " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels << ", " + tool2Name + " readID: " << bestAlntool2.readID << endl;
                    numnberOfBaseLine_FN++;
                    addBaselineFNAlnSz(bestAlntool2.partialMatchSize);
                } else if(baseLineReadPerQuery > 0 && tool2TP > 0){//compare the best alignmnets for this queyry
                    if (baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions > bestAlntool2.partialMatchSize){ //baseLine outperformed
                        tool2Best_baseLineOutperfomed++;
                        tool2Best_baseLineOutperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                        tool2Best_baseLineOutperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                        baseLineBest_baseLineOutperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                        baseLineBest_baseLineOutperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                        cmp << "Best: baseLine outperformed " + tool2Name + ", baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions <<  ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels <<  ", " + tool2Name + " readID: " << bestAlntool2.readID << endl;
                        addBaselineBestOutperform((baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions) - (bestAlntool2.partialMatchSize));
                    } else if (baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions < bestAlntool2.partialMatchSize){//tool2 outperformed
                        tool2Best_tool2Outperfomed++;
                        tool2Best_tool2Outperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                        tool2Best_tool2Outperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                        baseLineBest_tool2Outperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                        baseLineBest_tool2Outperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                        cmp << "Best: " + tool2Name + " outperformed baseLine in terms of alnlen, baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions <<  ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels << ", " + tool2Name + " readID: " << bestAlntool2.readID << endl;
                        addTool2BestOutperform((bestAlntool2.partialMatchSize) - (baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions));
                    } else {
                        if (baseLineBestAln.substitutions + baseLineBestAln.inDels < bestAlntool2.substitutions + bestAlntool2.inDels)//baseLine outperformed
                        {
                            tool2Best_baseLineOutperfomed++;
                            tool2Best_baseLineOutperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                            tool2Best_baseLineOutperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                            baseLineBest_baseLineOutperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                            baseLineBest_baseLineOutperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                            cmp << "Best: baseLine outperformed " + tool2Name + ", baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels <<  ", " + tool2Name + " readID: " << bestAlntool2.readID << endl;
                            addBaselineBestOutperform_edit((bestAlntool2.substitutions + bestAlntool2.inDels) - (baseLineBestAln.substitutions + baseLineBestAln.inDels));
                        } else if (baseLineBestAln.substitutions + baseLineBestAln.inDels > bestAlntool2.substitutions + bestAlntool2.inDels)//tool2 outperformed
                        {
                            tool2Best_tool2Outperfomed++;
                            tool2Best_tool2Outperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                            tool2Best_tool2Outperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                            baseLineBest_tool2Outperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                            baseLineBest_tool2Outperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                            cmp << "Best: " + tool2Name + " outperformed baseLine in terms of ed, baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions <<  ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels <<  ", " + tool2Name + " readID: " << bestAlntool2.readID << endl;
                            addTool2BestOutperform_edit((baseLineBestAln.substitutions + baseLineBestAln.inDels) - (bestAlntool2.substitutions + bestAlntool2.inDels));
                        } // tool2 and baseLine performed equally
                        else{
                            tool2Best_tool2EqualbaseLine++;
                            tool2Best_tool2EqualbaseLine_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                            tool2Best_tool2EqualbaseLine_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                            baseLineBest_tool2EqualbaseLine_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                            baseLineBest_tool2EqualbaseLine_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                            cmp << "Best: " + tool2Name + " & baseLine perfromed equally, baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels << ", " + tool2Name + " readID: " << bestAlntool2.readID << endl;
                        }
                    }
                }
            }
            tool2TP_nobaseLineMatches_allQ += tool2TP_nobaseLineMatches; 
            tool2TP_tool2Outperfomed_allQ += tool2TP_tool2Outperfomed; 
            tool2TP_baseLineOutperfomed_allQ += tool2TP_baseLineOutperfomed;
            tool2TP_tool2EqualbaseLine_allQ += tool2TP_tool2EqualbaseLine;
            cmp << tool2Name + "TP_nobaseLineMatches: " << tool2TP_nobaseLineMatches << ", " + tool2Name + "TP_tool2Outperfomed: " << tool2TP_tool2Outperfomed << ", " + tool2Name + "TP_baseLineOutperfomed: " << tool2TP_baseLineOutperfomed << ", " + tool2Name + "TP_tool2EqualbaseLine: " << tool2TP_tool2EqualbaseLine << endl;
            tool2FP_lowAlnLen_allQ += tool2FP_lowAlnLen; 
            tool2FP_lowPercentageIdentity_allQ += tool2FP_lowPercentageIdentity;
            cmp  << tool2Name + "FP_lowAlnLen: " << tool2FP_lowAlnLen << ", " + tool2Name << "FP_lowPercentageIdentity: " << tool2FP_lowPercentageIdentity << endl;
            tool2Best_baseLineOutperfomed_allQ += tool2Best_baseLineOutperfomed;
            tool2Best_tool2Outperfomed_allQ += tool2Best_tool2Outperfomed;
            tool2Best_tool2EqualbaseLine_allQ += tool2Best_tool2EqualbaseLine;
            cmp << tool2Name + "Best_baseLineOutperfomed: " << tool2Best_baseLineOutperfomed << ", " + tool2Name + "Best_tool2Outperfomed: " << tool2Best_tool2Outperfomed << ", " + tool2Name + "Best_tool2EqualbaseLine: " << tool2Best_tool2EqualbaseLine << endl;
            // best_aln_sz.push_back(make_pair(baseLineBestAlnSize, tool2BestAlnSize)); 
            //baseLine and tool2 alignments per query
            if(REPORT_ALN_PER_Q) alnPerQ << queryInd << " " << baselineReadID_Aln.size() << " " << tool2TP << endl;
        }
        Utilities<uint32_t> util;   
        cmp << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Overall Comparison Results>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cmp << left << setw(80) << "# Of queries : " + to_string(queryCount) << endl;
        cmp << left << setw(80) << "# Of queries that baseLine found match : " << queriesblfoundMatch << endl;
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " found match (TP) : " << queriestool2FoundMatch << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries that Baseline didn't find TN : " << totalbaseLineTN << endl;
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find TN : " << totalTool2TN << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find FN (Total) : " << tool2FN + tool2FN_noCriteriaFittedMatches << endl;
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find any match FN : " << tool2FN << endl;
        pair<uint32_t, uint32_t> avgTool2FN_baseLine_alnLen_allQ = util.calculateStatistics2(tool2FN_baseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgTool2FN_baseLine_ed_allQ = util.calculateStatistics2(tool2FN_baseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) Baseline's alignment length where " + tool2Name + " FN: " << avgTool2FN_baseLine_alnLen_allQ.first << ", " << avgTool2FN_baseLine_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) Baseline's edit distance where " + tool2Name + " FN: " << avgTool2FN_baseLine_ed_allQ.first << ", " << avgTool2FN_baseLine_ed_allQ.second << endl;
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find based on criteria FN : " << tool2FN_noCriteriaFittedMatches << endl;
        pair<uint32_t, uint32_t> avgtool2FN_noCriteria_baseLine_alnLen_allQ = util.calculateStatistics2(tool2FN_noCriteria_baseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2FN_noCriteria_baseLine_ed_allQ = util.calculateStatistics2(tool2FN_noCriteria_baseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where " + tool2Name + " FN based on criteria: " << avgtool2FN_noCriteria_baseLine_alnLen_allQ.first << ", " << avgtool2FN_noCriteria_baseLine_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where " + tool2Name + " FN based on criteria: " << avgtool2FN_noCriteria_baseLine_ed_allQ.first << ", " << avgtool2FN_noCriteria_baseLine_ed_allQ.second << endl;
       
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries that baseLine didn't find FN : " << numnberOfBaseLine_FN << endl;
        pair<uint32_t, uint32_t> avgBaselineFN_tool2_alnLen_allQ = util.calculateStatistics2(baseLineFN_tool2_alnLen_allQ);
        pair<uint32_t, uint32_t> avgBaselineFN_tool2_ed_allQ = util.calculateStatistics2(baseLineFN_tool2_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where baseLine FN: " << avgBaselineFN_tool2_alnLen_allQ.first << ", " << avgBaselineFN_tool2_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where baseLine FN: " << avgBaselineFN_tool2_ed_allQ.first << ", " << avgBaselineFN_tool2_ed_allQ.second << endl;
       
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of readIDs found by Baseline (total): " << baseLineTotalNumberOfReadIDs << endl;
        cmp << left << setw(80) << "# Of readIDs found by " + tool2Name + " (total): " << tool2TotalNumberOfReadIDs << endl;
        // cmp << left << setw(80) << "# Of readIDs found by " + tool2Name + " containing N (total): " << numberOfTool2ReadsContainingN << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of total " + tool2Name + " (TP): " << totaltool2TP << endl;
        pair<uint32_t, uint32_t> avgtool2TPPerQuery = util.calculateStatistics2(tool2TPPerQuery);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s TP per query: " << avgtool2TPPerQuery.first << ", " << avgtool2TPPerQuery.second << endl;

        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where baseLine didn't find match: " << tool2TP_nobaseLineMatches_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2TP_nobaseLineMatches_alnLen_allQ = util.calculateStatistics2(tool2TP_nobaseLineMatches_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2TP_nobaseLineMatches_ed_allQ = util.calculateStatistics2(tool2TP_nobaseLineMatches_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where baseLine didn't find match: " << avgtool2TP_nobaseLineMatches_alnLen_allQ.first << ", " << avgtool2TP_nobaseLineMatches_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where baseLine didn't find match: " << avgtool2TP_nobaseLineMatches_ed_allQ.first << ", " << avgtool2TP_nobaseLineMatches_ed_allQ.second << endl;
        
        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where baseLine outperformed: " << tool2TP_baseLineOutperfomed_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2TP_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(tool2TP_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2TP_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(tool2TP_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where baseLine outperformed: " << avgtool2TP_baseLineOutperfomed_alnLen_allQ.first << ", " << avgtool2TP_baseLineOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where baseLine outperformed: " << avgtool2TP_baseLineOutperfomed_ed_allQ.first << ", " << avgtool2TP_baseLineOutperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgbaseLineTP_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(baseLineTP_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineTP_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(baseLineTP_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where baseLine outperformed: " << avgbaseLineTP_baseLineOutperfomed_alnLen_allQ.first << ", " << avgbaseLineTP_baseLineOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where baseLine outperformed: " << avgbaseLineTP_baseLineOutperfomed_ed_allQ.first << ", " << avgbaseLineTP_baseLineOutperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where " + tool2Name + " outperformed: " << tool2TP_tool2Outperfomed_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2TP_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(tool2TP_tool2Outperfomed_alnLen_allQ);
    
        pair<uint32_t, uint32_t> avgtool2TP_tool2Outperfomed_ed_allQ = util.calculateStatistics2(tool2TP_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where " + tool2Name + " outperformed: " << avgtool2TP_tool2Outperfomed_alnLen_allQ.first << ", " << avgtool2TP_tool2Outperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where " + tool2Name + " outperformed: " << avgtool2TP_tool2Outperfomed_ed_allQ.first << ", " << avgtool2TP_tool2Outperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgbaseLineTP_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(baseLineTP_tool2Outperfomed_alnLen_allQ);
 
        pair<uint32_t, uint32_t> avgbaseLineTP_tool2Outperfomed_ed_allQ = util.calculateStatistics2(baseLineTP_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where " + tool2Name + " outperformed: " << avgbaseLineTP_tool2Outperfomed_alnLen_allQ.first << ", " << avgbaseLineTP_tool2Outperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where " + tool2Name + " outperformed: " << avgbaseLineTP_tool2Outperfomed_ed_allQ.first << ", " << avgbaseLineTP_tool2Outperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where they perfromed equaly: " << tool2TP_tool2EqualbaseLine_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2TP_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(tool2TP_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2TP_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(tool2TP_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where they perfromed equaly: " << avgtool2TP_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgtool2TP_tool2EqualbaseLine_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where they perfromed equaly: " << avgtool2TP_tool2EqualbaseLine_ed_allQ.first << ", " << avgtool2TP_tool2EqualbaseLine_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgbaseLineTP_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(baseLineTP_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineTP_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(baseLineTP_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where they perfromed equaly: " << avgbaseLineTP_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgbaseLineTP_tool2EqualbaseLine_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where they perfromed equaly: " << avgbaseLineTP_tool2EqualbaseLine_ed_allQ.first << ", " << avgbaseLineTP_tool2EqualbaseLine_ed_allQ.second << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of " + tool2Name + " best where baseLine outperformed: " << tool2Best_baseLineOutperfomed_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2Best_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(tool2Best_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2Best_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(tool2Best_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best alignment length where baseLine outperformed: " << avgtool2Best_baseLineOutperfomed_alnLen_allQ.first << ", " << avgtool2Best_baseLineOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best edit distance where baseLine outperformed: " << avgtool2Best_baseLineOutperfomed_ed_allQ.first << ", " << avgtool2Best_baseLineOutperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgbaseLineBest_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(baseLineBest_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineBest_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(baseLineBest_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's best alignment length where baseLine outperformed: " << avgbaseLineBest_baseLineOutperfomed_alnLen_allQ.first << ", " << avgbaseLineBest_baseLineOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's best edit distance where baseLine outperformed: " << avgbaseLineBest_baseLineOutperfomed_ed_allQ.first << ", " << avgbaseLineBest_baseLineOutperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of " + tool2Name + " best where " + tool2Name + " outperformed: " << tool2Best_tool2Outperfomed_allQ  << endl;
        pair<uint32_t, uint32_t> avgtool2Best_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(tool2Best_tool2Outperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2Best_tool2Outperfomed_ed_allQ = util.calculateStatistics2(tool2Best_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best alignment length where " + tool2Name + " outperformed: " << avgtool2Best_tool2Outperfomed_alnLen_allQ.first << ", " << avgtool2Best_tool2Outperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best edit distance where " + tool2Name + " outperformed: " << avgtool2Best_tool2Outperfomed_ed_allQ.first << ", " << avgtool2Best_tool2Outperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(baseLineBest_tool2Outperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2Outperfomed_ed_allQ = util.calculateStatistics2(baseLineBest_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's best alignment length where " + tool2Name + " outperformed: " << avgbaseLineBest_tool2Outperfomed_alnLen_allQ.first << ", " << avgbaseLineBest_tool2Outperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's best edit distance where " + tool2Name + " outperformed: " << avgbaseLineBest_tool2Outperfomed_ed_allQ.first << ", " << avgbaseLineBest_tool2Outperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of " + tool2Name + " best where they perfromed equaly: " << tool2Best_tool2EqualbaseLine_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2Best_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(tool2Best_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2Best_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(tool2Best_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best alignment length where they perfromed equaly: " << avgtool2Best_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgtool2Best_tool2EqualbaseLine_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best edit distance where they perfromed equaly: " << avgtool2Best_tool2EqualbaseLine_ed_allQ.first << ", " << avgtool2Best_tool2EqualbaseLine_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(baseLineBest_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(baseLineBest_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's best alignment length where they perfromed equaly: " << avgbaseLineBest_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgbaseLineBest_tool2EqualbaseLine_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's best edit distance where they perfromed equaly: " << avgbaseLineBest_tool2EqualbaseLine_ed_allQ.first << ", " << avgbaseLineBest_tool2EqualbaseLine_ed_allQ.second << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of total " + tool2Name + " (FP): " << totaltool2FP << endl;
        pair<uint32_t, uint32_t> avgtool2FPPerQuery = util.calculateStatistics2(tool2FPPerQuery);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s FP per query: " << avgtool2FPPerQuery.first << ", " << avgtool2FPPerQuery.second << endl;

        cmp << left << setw(80) << "# Of " + tool2Name + " (FP) that alignment length < R: " << tool2FP_lowAlnLen_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2FP_lowAlnLen_alnLen_allQ = util.calculateStatistics2(tool2FP_lowAlnLen_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2FP_lowAlnLen_ed_allQ = util.calculateStatistics2(tool2FP_lowAlnLen_ed_allQ);
        pair<uint32_t, uint32_t> avgtool2FP_lowAlnLen_matchsize_allQ = util.calculateStatistics2(tool2FP_lowAlnLen_matchsize_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where alignment length < R: " << avgtool2FP_lowAlnLen_alnLen_allQ.first << ", " << avgtool2FP_lowAlnLen_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where alignment length < R: " << avgtool2FP_lowAlnLen_ed_allQ.first << ", " << avgtool2FP_lowAlnLen_ed_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s match size where alignment length < R: " << avgtool2FP_lowAlnLen_matchsize_allQ.first << ", " << avgtool2FP_lowAlnLen_matchsize_allQ.second << endl;

        // cmp << left << setw(80) << "# Of " + tool2Name + " (FP) that # of kmers is lower than regionSize" << tool2FP_lowerExactMatchKmers_allQ << endl;
        // pair<uint32_t, uint32_t> avgtool2FP_lowerExactMatchKmers_alnLen_allQ = util.calculateStatistics2(tool2FP_lowerExactMatchKmers_alnLen_allQ);
        // pair<uint32_t, uint32_t> avgtool2FP_lowerExactMatchKmers_ed_allQ = util.calculateStatistics2(tool2FP_lowerExactMatchKmers_ed_allQ);
        // cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where # of kmers < MinExactNumberOfKmers" <<  avgtool2FP_lowerExactMatchKmers_alnLen_allQ.first << ", " << avgtool2FP_lowerExactMatchKmers_alnLen_allQ.second << endl;
        // cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where # of kmers < MinExactNumberOfKmers" << avgtool2FP_lowerExactMatchKmers_ed_allQ.first << ", " << avgtool2FP_lowerExactMatchKmers_ed_allQ.second << endl;
        
        cmp << left << setw(80) << "# Of " + tool2Name + " (FP) that percentage identity is lower than PI" << tool2FP_lowPercentageIdentity_allQ << endl;
        pair<uint32_t, uint32_t> avgtool2FP_lowPercentageIdentity_alnLen_allQ = util.calculateStatistics2(tool2FP_lowPercentageIdentity_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2FP_lowPercentageIdentity_ed_allQ = util.calculateStatistics2(tool2FP_lowPercentageIdentity_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where percentage identity is lower than PI" <<  avgtool2FP_lowPercentageIdentity_alnLen_allQ.first << ", " << avgtool2FP_lowPercentageIdentity_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where percentage identity is lower than PI" << avgtool2FP_lowPercentageIdentity_ed_allQ.first << ", " << avgtool2FP_lowPercentageIdentity_ed_allQ.second << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of baseLine (TP) where " + tool2Name + " didn't find match: " << baseLineTP_notool2Matches_allQ << endl;
        pair<uint32_t, uint32_t> avgbaseLineTP_notool2Matches_alnLen_allQ = util.calculateStatistics2(baseLineTP_notool2Matches_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineTP_notool2Matches_ed_allQ = util.calculateStatistics2(baseLineTP_notool2Matches_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where " + tool2Name + " didn't find match: " << avgbaseLineTP_notool2Matches_alnLen_allQ.first << ", " << avgbaseLineTP_notool2Matches_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where " + tool2Name + " didn't find match: " << avgbaseLineTP_notool2Matches_ed_allQ.first << ", " << avgbaseLineTP_notool2Matches_ed_allQ.second << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries contains N: " << numberOfQueryContainN << endl;

        cmp << "------------------------------------------------------------------------------------------" << endl;
        pair<uint32_t, uint32_t> baseLineReadPerQuery_allQ  = util.calculateStatistics2(baseLineReadPerQuerySet);
        pair<uint32_t, uint32_t> tool2ReadPerQuery_allQ  = util.calculateStatistics2(tool2ReadPerQuerySet);
        cmp << left << setw(80) << "(avg, median) baseLineS's Matches per Query: " << baseLineReadPerQuery_allQ.first << ", " << baseLineReadPerQuery_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s Matches per Query: " << tool2ReadPerQuery_allQ.first << ", " << tool2ReadPerQuery_allQ.second << endl;
        cmp << "------------------------------------------------------------------------------------------" << endl;

        for(uint32_t i = 1; i < cfg.contigSize; i++) {
            tpAlnCmp << i << "\t" << baselineTPOutperform_aln_sz_difference[i] << "\t" << tool2TPOutperform_aln_sz_difference[i] << endl;
            bestAlnCmp << i << "\t" << baselineBestOutperform_aln_sz_difference[i] << "\t" << tool2BestOutperform_aln_sz_difference[i] << endl;
            tpEdCmp << i << "\t" << baselineTPOutperform_edit_difference[i] << "\t" << tool2TPOutperform_edit_difference[i] << endl;
            bestEdCmp << i << "\t" << baselineBestOutperform_edit_difference[i] << "\t" << tool2BestOutperform_edit_difference[i] << endl;
        }

        for(uint32_t i = 1; i < cfg.contigSize + 20; i++) {
            tpAlnSz << i << "\t" << baselineTP_aln_sz[i] << "\t" << tool2TP_aln_sz[i] << endl;
            bestAlnSz << i << "\t" << baselineBest_aln_sz[i] << "\t" << tool2Best_aln_sz[i] << endl;
            fnAlnSz << i << "\t" << baselineFN_aln_sz[i] << "\t" << tool2FN_aln_sz[i] << endl;
        }
        cmp.close();
        if(REPORT_ALN_PER_Q) alnPerQ.close();
        bestAlnCmp.close();
        tpAlnCmp.close();
        tpEdCmp.close();
        bestEdCmp.close();
        tpAlnSz.close();
        bestAlnSz.close();
        fnAlnSz.close();
        if (REPORT_baseLine_FN) baseLineFn.close();
    }
};

#endif