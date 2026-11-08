#ifndef COMAPREWITHBASELINE_H
#define COMAPREWITHBASELINE_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include <exception>
#include <sstream>
#include <cmath>
#include <omp.h>
#include "BlastReader.h"
#include "Aligner.h"
#include "Alignment.h"
#include "SamReader.h"
#include "Minimap2Reader.h"
#include "Container.h"

#define REPORT_BEST_ALN 0
#define REPORT_baseLine_FN 0
#define REPORT_ALN_PER_Q 0
#define CHECK_EXACT_MATCH_CRITERION 1

using namespace std;


class CompareWithBaseLine {
private:
    double percentageIdentity;
    struct Statistics {
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
        map<uint32_t, int64_t> thresholdFNChanges;
        uint64_t numberOfQueryContainN = 0;
        uint64_t queriestool2FoundMatch = 0;
        uint64_t queriesblfoundMatch = 0;
        uint64_t baseLineTotalNumberOfReadIDs = 0;
        uint64_t tool2TotalNumberOfReadIDs = 0;
        uint64_t numnberOfBaseLine_FN = 0;
        uint64_t tool2FN = 0;
        uint64_t tool2FN_noCriteriaFittedMatches = 0;
        uint64_t totalTool2TN = 0;
        uint64_t totalbaseLineTN = 0;
        uint64_t totaltool2TP = 0;
        uint64_t totaltool2FP = 0;
        uint64_t tool2TP_nobaseLineMatches_allQ = 0;
        uint64_t tool2TP_tool2Outperfomed_allQ = 0;
        uint64_t tool2TP_baseLineOutperfomed_allQ = 0;
        uint64_t tool2TP_tool2EqualbaseLine_allQ = 0;
        uint64_t baseLineTP_notool2Matches_allQ = 0;
        uint64_t tool2NoExactRun_allQ = 0;
        uint64_t tool2FP_lowPercentageIdentity_allQ = 0;
        uint64_t tool2Best_tool2Outperfomed_allQ = 0;
        uint64_t tool2Best_baseLineOutperfomed_allQ = 0;
        uint64_t tool2Best_tool2EqualbaseLine_allQ = 0;
        multiset<uint32_t> baseLineReadPerQuerySet;
        multiset<uint32_t> tool2ReadPerQuerySet;
        multiset<uint32_t> tool2TPPerQuery;
        multiset<uint32_t> tool2FPPerQuery;
        multiset<uint32_t> tool2TP_nobaseLineMatches_alnLen_allQ;
        multiset<uint32_t> tool2TP_tool2Outperfomed_alnLen_allQ;
        multiset<uint32_t> tool2TP_baseLineOutperfomed_alnLen_allQ;
        multiset<uint32_t> tool2TP_tool2EqualbaseLine_alnLen_allQ;
        multiset<uint32_t> tool2TP_nobaseLineMatches_ed_allQ;
        multiset<uint32_t> tool2TP_tool2Outperfomed_ed_allQ;
        multiset<uint32_t> tool2TP_baseLineOutperfomed_ed_allQ;
        multiset<uint32_t> tool2TP_tool2EqualbaseLine_ed_allQ;
        multiset<uint32_t> baseLineTP_notool2Matches_alnLen_allQ;
        multiset<uint32_t> baseLineTP_notool2Matches_ed_allQ;
        multiset<uint32_t> baseLineTP_tool2Outperfomed_alnLen_allQ;
        multiset<uint32_t> baseLineTP_baseLineOutperfomed_alnLen_allQ;
        multiset<uint32_t> baseLineTP_tool2EqualbaseLine_alnLen_allQ;
        multiset<uint32_t> baseLineTP_tool2Outperfomed_ed_allQ;
        multiset<uint32_t> baseLineTP_baseLineOutperfomed_ed_allQ;
        multiset<uint32_t> baseLineTP_tool2EqualbaseLine_ed_allQ;
        multiset<uint32_t> tool2Best_tool2Outperfomed_alnLen_allQ;
        multiset<uint32_t> tool2Best_baseLineOutperfomed_alnLen_allQ;
        multiset<uint32_t> tool2Best_tool2EqualbaseLine_alnLen_allQ;
        multiset<uint32_t> tool2Best_tool2Outperfomed_ed_allQ;
        multiset<uint32_t> tool2Best_baseLineOutperfomed_ed_allQ;
        multiset<uint32_t> tool2Best_tool2EqualbaseLine_ed_allQ;
        multiset<uint32_t> baseLineBest_tool2Outperfomed_alnLen_allQ;
        multiset<uint32_t> baseLineBest_baseLineOutperfomed_alnLen_allQ;
        multiset<uint32_t> baseLineBest_tool2EqualbaseLine_alnLen_allQ;
        multiset<uint32_t> baseLineBest_tool2Outperfomed_ed_allQ;
        multiset<uint32_t> baseLineBest_baseLineOutperfomed_ed_allQ;
        multiset<uint32_t> baseLineBest_tool2EqualbaseLine_ed_allQ;
        multiset<uint32_t> tool2NoExactRun_alnLen_allQ;
        multiset<uint32_t> tool2NoExactRun_matchsize_allQ;
        multiset<uint32_t> tool2FP_lowerExactMatchKmers_alnLen_allQ;
        multiset<uint32_t> tool2FP_lowPercentageIdentity_alnLen_allQ;
        multiset<uint32_t> tool2NoExactRun_ed_allQ;
        multiset<uint32_t> tool2FP_lowerExactMatchKmers_ed_allQ;
        multiset<uint32_t> tool2FP_lowPercentageIdentity_ed_allQ;
        multiset<uint32_t> tool2FN_baseLine_alnLen_allQ;
        multiset<uint32_t> tool2FN_baseLine_ed_allQ;
        multiset<uint32_t> tool2FN_noCriteria_baseLine_alnLen_allQ;
        multiset<uint32_t> tool2FN_noCriteria_tool1_ed_allQ;
        multiset<uint32_t> tool2FN_noCriteria_baseLine_ed_allQ;
        multiset<uint32_t> baseLineFN_tool2_alnLen_allQ;
        multiset<uint32_t> baseLineFN_tool2_ed_allQ;

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

        void merge(Statistics& other) {
            for (const auto& [key, value] : other.tool2FN_aln_sz) tool2FN_aln_sz[key] += value;
            for (const auto& [key, value] : other.baselineFN_aln_sz) baselineFN_aln_sz[key] += value;
            for (const auto& [key, value] : other.tool2TP_aln_sz) tool2TP_aln_sz[key] += value;
            for (const auto& [key, value] : other.baselineTP_aln_sz) baselineTP_aln_sz[key] += value;
            for (const auto& [key, value] : other.tool2Best_aln_sz) tool2Best_aln_sz[key] += value;
            for (const auto& [key, value] : other.baselineBest_aln_sz) baselineBest_aln_sz[key] += value;
            for (const auto& [key, value] : other.tool2TPOutperform_aln_sz_difference) tool2TPOutperform_aln_sz_difference[key] += value;
            for (const auto& [key, value] : other.baselineTPOutperform_aln_sz_difference) baselineTPOutperform_aln_sz_difference[key] += value;
            for (const auto& [key, value] : other.tool2TPOutperform_edit_difference) tool2TPOutperform_edit_difference[key] += value;
            for (const auto& [key, value] : other.baselineTPOutperform_edit_difference) baselineTPOutperform_edit_difference[key] += value;
            for (const auto& [key, value] : other.tool2BestOutperform_aln_sz_difference) tool2BestOutperform_aln_sz_difference[key] += value;
            for (const auto& [key, value] : other.baselineBestOutperform_aln_sz_difference) baselineBestOutperform_aln_sz_difference[key] += value;
            for (const auto& [key, value] : other.tool2BestOutperform_edit_difference) tool2BestOutperform_edit_difference[key] += value;
            for (const auto& [key, value] : other.baselineBestOutperform_edit_difference) baselineBestOutperform_edit_difference[key] += value;
            for (const auto& [key, value] : other.thresholdFNChanges) thresholdFNChanges[key] += value;
            numberOfQueryContainN += other.numberOfQueryContainN;
            queriestool2FoundMatch += other.queriestool2FoundMatch;
            queriesblfoundMatch += other.queriesblfoundMatch;
            baseLineTotalNumberOfReadIDs += other.baseLineTotalNumberOfReadIDs;
            tool2TotalNumberOfReadIDs += other.tool2TotalNumberOfReadIDs;
            numnberOfBaseLine_FN += other.numnberOfBaseLine_FN;
            tool2FN += other.tool2FN;
            tool2FN_noCriteriaFittedMatches += other.tool2FN_noCriteriaFittedMatches;
            totalTool2TN += other.totalTool2TN;
            totalbaseLineTN += other.totalbaseLineTN;
            totaltool2TP += other.totaltool2TP;
            totaltool2FP += other.totaltool2FP;
            tool2TP_nobaseLineMatches_allQ += other.tool2TP_nobaseLineMatches_allQ;
            tool2TP_tool2Outperfomed_allQ += other.tool2TP_tool2Outperfomed_allQ;
            tool2TP_baseLineOutperfomed_allQ += other.tool2TP_baseLineOutperfomed_allQ;
            tool2TP_tool2EqualbaseLine_allQ += other.tool2TP_tool2EqualbaseLine_allQ;
            baseLineTP_notool2Matches_allQ += other.baseLineTP_notool2Matches_allQ;
            tool2NoExactRun_allQ += other.tool2NoExactRun_allQ;
            tool2FP_lowPercentageIdentity_allQ += other.tool2FP_lowPercentageIdentity_allQ;
            tool2Best_tool2Outperfomed_allQ += other.tool2Best_tool2Outperfomed_allQ;
            tool2Best_baseLineOutperfomed_allQ += other.tool2Best_baseLineOutperfomed_allQ;
            tool2Best_tool2EqualbaseLine_allQ += other.tool2Best_tool2EqualbaseLine_allQ;
            // Transfer observation nodes without copying or allocating them.
            baseLineReadPerQuerySet.merge(other.baseLineReadPerQuerySet);
            tool2ReadPerQuerySet.merge(other.tool2ReadPerQuerySet);
            tool2TPPerQuery.merge(other.tool2TPPerQuery);
            tool2FPPerQuery.merge(other.tool2FPPerQuery);
            tool2TP_nobaseLineMatches_alnLen_allQ.merge(other.tool2TP_nobaseLineMatches_alnLen_allQ);
            tool2TP_tool2Outperfomed_alnLen_allQ.merge(other.tool2TP_tool2Outperfomed_alnLen_allQ);
            tool2TP_baseLineOutperfomed_alnLen_allQ.merge(other.tool2TP_baseLineOutperfomed_alnLen_allQ);
            tool2TP_tool2EqualbaseLine_alnLen_allQ.merge(other.tool2TP_tool2EqualbaseLine_alnLen_allQ);
            tool2TP_nobaseLineMatches_ed_allQ.merge(other.tool2TP_nobaseLineMatches_ed_allQ);
            tool2TP_tool2Outperfomed_ed_allQ.merge(other.tool2TP_tool2Outperfomed_ed_allQ);
            tool2TP_baseLineOutperfomed_ed_allQ.merge(other.tool2TP_baseLineOutperfomed_ed_allQ);
            tool2TP_tool2EqualbaseLine_ed_allQ.merge(other.tool2TP_tool2EqualbaseLine_ed_allQ);
            baseLineTP_notool2Matches_alnLen_allQ.merge(other.baseLineTP_notool2Matches_alnLen_allQ);
            baseLineTP_notool2Matches_ed_allQ.merge(other.baseLineTP_notool2Matches_ed_allQ);
            baseLineTP_tool2Outperfomed_alnLen_allQ.merge(other.baseLineTP_tool2Outperfomed_alnLen_allQ);
            baseLineTP_baseLineOutperfomed_alnLen_allQ.merge(other.baseLineTP_baseLineOutperfomed_alnLen_allQ);
            baseLineTP_tool2EqualbaseLine_alnLen_allQ.merge(other.baseLineTP_tool2EqualbaseLine_alnLen_allQ);
            baseLineTP_tool2Outperfomed_ed_allQ.merge(other.baseLineTP_tool2Outperfomed_ed_allQ);
            baseLineTP_baseLineOutperfomed_ed_allQ.merge(other.baseLineTP_baseLineOutperfomed_ed_allQ);
            baseLineTP_tool2EqualbaseLine_ed_allQ.merge(other.baseLineTP_tool2EqualbaseLine_ed_allQ);
            tool2Best_tool2Outperfomed_alnLen_allQ.merge(other.tool2Best_tool2Outperfomed_alnLen_allQ);
            tool2Best_baseLineOutperfomed_alnLen_allQ.merge(other.tool2Best_baseLineOutperfomed_alnLen_allQ);
            tool2Best_tool2EqualbaseLine_alnLen_allQ.merge(other.tool2Best_tool2EqualbaseLine_alnLen_allQ);
            tool2Best_tool2Outperfomed_ed_allQ.merge(other.tool2Best_tool2Outperfomed_ed_allQ);
            tool2Best_baseLineOutperfomed_ed_allQ.merge(other.tool2Best_baseLineOutperfomed_ed_allQ);
            tool2Best_tool2EqualbaseLine_ed_allQ.merge(other.tool2Best_tool2EqualbaseLine_ed_allQ);
            baseLineBest_tool2Outperfomed_alnLen_allQ.merge(other.baseLineBest_tool2Outperfomed_alnLen_allQ);
            baseLineBest_baseLineOutperfomed_alnLen_allQ.merge(other.baseLineBest_baseLineOutperfomed_alnLen_allQ);
            baseLineBest_tool2EqualbaseLine_alnLen_allQ.merge(other.baseLineBest_tool2EqualbaseLine_alnLen_allQ);
            baseLineBest_tool2Outperfomed_ed_allQ.merge(other.baseLineBest_tool2Outperfomed_ed_allQ);
            baseLineBest_baseLineOutperfomed_ed_allQ.merge(other.baseLineBest_baseLineOutperfomed_ed_allQ);
            baseLineBest_tool2EqualbaseLine_ed_allQ.merge(other.baseLineBest_tool2EqualbaseLine_ed_allQ);
            tool2NoExactRun_alnLen_allQ.merge(other.tool2NoExactRun_alnLen_allQ);
            tool2NoExactRun_matchsize_allQ.merge(other.tool2NoExactRun_matchsize_allQ);
            tool2FP_lowerExactMatchKmers_alnLen_allQ.merge(other.tool2FP_lowerExactMatchKmers_alnLen_allQ);
            tool2FP_lowPercentageIdentity_alnLen_allQ.merge(other.tool2FP_lowPercentageIdentity_alnLen_allQ);
            tool2NoExactRun_ed_allQ.merge(other.tool2NoExactRun_ed_allQ);
            tool2FP_lowerExactMatchKmers_ed_allQ.merge(other.tool2FP_lowerExactMatchKmers_ed_allQ);
            tool2FP_lowPercentageIdentity_ed_allQ.merge(other.tool2FP_lowPercentageIdentity_ed_allQ);
            tool2FN_baseLine_alnLen_allQ.merge(other.tool2FN_baseLine_alnLen_allQ);
            tool2FN_baseLine_ed_allQ.merge(other.tool2FN_baseLine_ed_allQ);
            tool2FN_noCriteria_baseLine_alnLen_allQ.merge(other.tool2FN_noCriteria_baseLine_alnLen_allQ);
            tool2FN_noCriteria_tool1_ed_allQ.merge(other.tool2FN_noCriteria_tool1_ed_allQ);
            tool2FN_noCriteria_baseLine_ed_allQ.merge(other.tool2FN_noCriteria_baseLine_ed_allQ);
            baseLineFN_tool2_alnLen_allQ.merge(other.baseLineFN_tool2_alnLen_allQ);
            baseLineFN_tool2_ed_allQ.merge(other.baseLineFN_tool2_ed_allQ);
        }
    };
public:
    CompareWithBaseLine(double pi) : percentageIdentity(pi) {
        if (!std::isfinite(pi) || pi < 0 || pi > 1)
            throw invalid_argument("Percentage identity must be between 0 and 1");
    }

    void determineEditsLocationsAndType(string read, string query, vector<int>& editPos) {
        // Check if the sequences have the same length
        if (read.length() != query.length()) {
            // cout << "q : " << query << endl;
            // cout << "r : " << read << endl;
            cerr << "Error: Sequences have different lengths." << '\n';
            return;
        }
        // Iterate through the sequences and compare corresponding characters
        for (size_t i = 0; i < read.length(); ++i) {
            if (read[i] != query[i]) {
                editPos.push_back(i);
            }
        }

    }

    string blastGetCigarStr(uint32_t queryS, string queryAligned, string readAligned, uint32_t contigSize, bool skipClips = false) const
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
        cout << "Sam: Read " << baseLineSamAlignments.size() << " alignments." << '\n';
    }

    string convertCigarToStr(const string& cigar, bool skipClips = false) const {
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
                case 'x':
                    // Match or mismatch (substitution)
                    simplifiedCigar << string(num, 'X');
                    break;
                case 'I':
                case 'i':
                    // Insertion
                    simplifiedCigar << string(num, 'I');
                    break;
                case 'D':
                case 'd':
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

    string getCigarStr(Alignment aln, Config cfg, string toolName) const{
        if(toolName == "BLAST" || toolName == "blast") {
            return blastGetCigarStr(aln.queryRegionStartPos, aln.alignedQuery, aln.alignedRead, cfg.contigSize, true);
        } else {
            return convertCigarToStr(aln.cigar, true);
        }
    }

    bool hasConsecutiveMatches(const std::string& s, int k) const {
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

    uint32_t getMatchesCount(string cigarStr) const {
        uint32_t cnt = 0;
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == '=') {
                cnt++;
            }
        }
        return cnt;
    }

    string trimClips(string cigarStr) const{
        string trimmedCigarStr = "";
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'S' || cigarStr[i] == 'H') {
                continue;
            }
            trimmedCigarStr += cigarStr[i];
        }
        return trimmedCigarStr;
    }

    bool checkIdentityPercentange(string cigarStr, bool isTrimClips = false) const {
        if (isTrimClips) cigarStr = trimClips(cigarStr);
        if (cigarStr.empty()) return false;
        uint32_t matches = getMatchesCount(cigarStr);
        double identity =  (double) matches / (double) cigarStr.size();
        if(DEBUG_MODE) cout << "matches: " << matches << ", cigarStr.size : " << cigarStr.size() << ", identity: " << identity << '\n';
        if(identity < percentageIdentity)
            return false;
        return true;
    }

    bool checkIdentityPercentangeWithAln(Alignment aln) const {
        if (aln.partialMatchSize == 0) return false;
        double identity =  (double) aln.matches / (double) aln.partialMatchSize;
        if(identity < percentageIdentity)
            return false;
        return true;
    }

    int longestRunOfMatches(const std::string& cigarStr) const {
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
private:
    // Inputs are immutable; each worker owns its accumulated statistics and output.
    void compareQuery(const Config& cfg, uint32_t queryInd,
                      const tsl::robin_map<uint32_t, string>& queries,
                      const IndexContainer<uint32_t, Alignment>& queryTool2Alignment,
                      const IndexContainer<uint32_t, Alignmen_Baseline>& queryBaselineAlignment,
                      Statistics& stats, ostream& cmp, ostream& alnPerQ,
                      ostream& baseLineFn) const {
        const string& tool2Name = cfg.otherTool;
        const bool isBwotie = tool2Name == "bowtie" || tool2Name == "BOWTIE";
        const bool isMinimap2 = tool2Name == "minimap2" || tool2Name == "MINIMAP2";
        uint32_t baseLineBestAlnSize = 0;
        const auto queryIt = queries.find(queryInd);
        if (queryIt == queries.end())
            throw runtime_error("Missing query ID: " + to_string(queryInd));
        const string& query = queryIt->second;
        if(query.find('N') != string::npos || query.find('n') != string::npos)
        {
            if (!cfg.comparisonSummaryOnly) cmp << "query contains N!" << '\n';
            stats.numberOfQueryContainN++;
            if(REPORT_ALN_PER_Q) alnPerQ << queryInd << " 0 0"<< '\n';
            return;
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
                    stats.addBaselineTPAlnSz(aln.matches + aln.inDels + aln.substitutions);
                } else {
                    auto readAln = readIt->second;
                    if ((readAln.matches + readAln.inDels + readAln.substitutions < aln.matches + aln.inDels + aln.substitutions) ||
                    ((readAln.matches + readAln.inDels + readAln.substitutions == aln.matches + aln.inDels + aln.substitutions) &&
                    (readAln.inDels + readAln.substitutions > aln.inDels + aln.substitutions))) {
                        baselineReadID_Aln[aln.readID] = aln;
                        stats.removeBaselineTPAlnSz(readAln.matches + readAln.inDels + readAln.substitutions);
                        stats.addBaselineTPAlnSz(aln.matches + aln.inDels + aln.substitutions);
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
        stats.baseLineTotalNumberOfReadIDs += baseLineReadPerQuery;
        stats.baseLineReadPerQuerySet.insert(baseLineReadPerQuery);
        if(baseLineReadPerQuery > 0){
            stats.queriesblfoundMatch++;
            baseLineBestAlnSize = baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions;
            stats.addBaselineBestAlnSz(baseLineBestAlnSize);
        }
        //query level parameters
        uint32_t tool2TP = 0, tool2FP = 0;
        uint32_t tool2TP_nobaseLineMatches = 0, tool2TP_tool2Outperfomed = 0, tool2TP_baseLineOutperfomed = 0, tool2TP_tool2EqualbaseLine = 0;
        uint32_t baseLineTP_notool2Matches = 0;
        uint32_t tool2Best_tool2Outperfomed = 0, tool2Best_baseLineOutperfomed = 0, tool2Best_tool2EqualbaseLine = 0;
        uint32_t tool2NoExactRun = 0, tool2FP_lowPercentageIdentity = 0;
        if (!cfg.comparisonSummaryOnly) cmp << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        if (!cfg.comparisonSummaryOnly) cmp << "Q : " << query << ", queryInd: " << queryInd << '\n';
        Alignment bestAlntool2;

        struct Candidate {
            Alignment alignment;
            bool passesIdentity;
        };
        tsl::robin_map<uint32_t, Candidate> query_t2Alignments;
        auto tool2QueryRrange = queryTool2Alignment.getRange(queryInd);
        for (auto it = tool2QueryRrange.first; it != tool2QueryRrange.second; it++){
            const Alignment& aln = it->second;
            const bool passesIdentity = (isMinimap2 || isBwotie)
                ? checkIdentityPercentangeWithAln(aln)
                : checkIdentityPercentange(getCigarStr(aln, cfg, tool2Name));
            auto readIt = query_t2Alignments.find(aln.readID);
            if(readIt == query_t2Alignments.end()) {
                query_t2Alignments[aln.readID] = {aln, passesIdentity};
            } else {
                const auto& previous = readIt->second;
                const auto& readAln = previous.alignment;
                // Choose the longest accepted alignment for each pair. If all
                // candidates fail identity, retain the best rejected one for FP reporting.
                if ((passesIdentity && !previous.passesIdentity) ||
                    (passesIdentity == previous.passesIdentity &&
                     (aln.partialMatchSize > readAln.partialMatchSize ||
                      (aln.partialMatchSize == readAln.partialMatchSize &&
                       aln.inDels + aln.substitutions < readAln.inDels + readAln.substitutions)))) {
                    readIt.value() = {aln, passesIdentity};
                }
            }
        }
        size_t tool2ReadPerQuery = query_t2Alignments.size();
        stats.tool2TotalNumberOfReadIDs += tool2ReadPerQuery;
        stats.tool2ReadPerQuerySet.insert(tool2ReadPerQuery);
        if (!cfg.comparisonSummaryOnly) cmp << "# of reads found by " + tool2Name + " : " << tool2ReadPerQuery << '\n';
        if (!cfg.comparisonSummaryOnly) cmp << "# of reads found by Baseline : " << baseLineReadPerQuery << '\n';

        unordered_map<uint32_t, uint32_t> acceptedLengths;
        if(tool2ReadPerQuery == 0 && baseLineReadPerQuery == 0){
            //TN
            stats.totalTool2TN++;
            stats.totalbaseLineTN++;
            if (!cfg.comparisonSummaryOnly) cmp << "TN for " + tool2Name << '\n';
        } else if(tool2ReadPerQuery == 0 && baseLineReadPerQuery > 0) {
            //tool2 FN
            stats.tool2FN++;
            stats.tool2FN_baseLine_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
            stats.tool2FN_baseLine_ed_allQ.insert(baseLineBestAln.inDels + baseLineBestAln.substitutions);
            if (!cfg.comparisonSummaryOnly) cmp << "FN for " + tool2Name + ", Baseline alnlen : " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", ed : " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", readID: " <<  baseLineBestAln.readID << '\n';
        } else {    // if both tool2 and baseLine have alignments or only tool2 has alignments
            for (auto it = query_t2Alignments.begin(); it != query_t2Alignments.end(); it++)
            {
                Alignment aln = it->second.alignment;
                // string tool2R = reads[aln.readID];
                // if(tool2R.find('N') != string::npos || tool2R.find('n') != string::npos)
                // {
                //     numberOfTool2ReadsContainingN++;
                //     continue;
                // }
                bool isFP = false;
                string cigarStr = getCigarStr(aln, cfg, tool2Name);
                // M + NM gives identity but not mismatch positions. Exact-run
                // diagnostics must not interpret ambiguous M columns as exact matches.
                const bool exactPositionsKnown = !isMinimap2 &&
                    !(isBwotie && aln.cigar.find('M') != string::npos);
                if(exactPositionsKnown && CHECK_EXACT_MATCH_CRITERION && !hasConsecutiveMatches(cigarStr, cfg.kmerLength)) {
                    // A diagnostic only; acceptance is determined by identity.
                    aln.criteriaCode = 0x10;
                    tool2NoExactRun++;
                    uint16_t ltm = longestRunOfMatches(cigarStr);
                    stats.tool2NoExactRun_alnLen_allQ.insert(aln.partialMatchSize);
                    stats.tool2NoExactRun_ed_allQ.insert(aln.substitutions + aln.inDels);
                    stats.tool2NoExactRun_matchsize_allQ.insert(ltm);
                    if (!cfg.comparisonSummaryOnly) cmp << "Diagnostic for " + tool2Name + ", no exact match run of length " + to_string(cfg.kmerLength) + ", alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", match size: " << ltm << '\n';
                }
                // The exact-run diagnostic does not exempt a pair from identity filtering.
                if (!it->second.passesIdentity) {
                    isFP = true;
                    aln.criteriaCode |= 0x80;
                    tool2FP_lowPercentageIdentity++;
                    stats.tool2FP_lowPercentageIdentity_alnLen_allQ.insert(aln.partialMatchSize);
                    stats.tool2FP_lowPercentageIdentity_ed_allQ.insert(aln.substitutions + aln.inDels);
                    if (!cfg.comparisonSummaryOnly) cmp << "FP for " + tool2Name + ", due to low PercentageIdentity, alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << '\n';
                }

                if (isFP){
                    tool2FP++;
                }
                else {
                    tool2TP++;
                    acceptedLengths[aln.readID] = aln.partialMatchSize;
                    stats.addTool2TPAlnSz(aln.partialMatchSize);
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
                        stats.tool2TP_nobaseLineMatches_alnLen_allQ.insert(aln.partialMatchSize);
                        stats.tool2TP_nobaseLineMatches_ed_allQ.insert(aln.substitutions + aln.inDels);
                        if (!cfg.comparisonSummaryOnly) cmp << "TP for " + tool2Name + ", no baseLine match, alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << '\n';
                        stats.addBaselineFNAlnSz(aln.partialMatchSize);
                        if(REPORT_baseLine_FN) baseLineFn << queryInd << "\t" << aln.readID << "\t" << aln.flag << '\n';
                    } else {
                        if (baselineAln.matches + baselineAln.inDels + baselineAln.substitutions > aln.partialMatchSize){ //baseLine outperformed
                            tool2TP_baseLineOutperfomed++;
                            stats.tool2TP_baseLineOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                            stats.tool2TP_baseLineOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                            stats.baseLineTP_baseLineOutperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                            stats.baseLineTP_baseLineOutperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                            stats.addBaselineTPOutperform(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions - aln.partialMatchSize);
                        } else if (baselineAln.matches + baselineAln.inDels + baselineAln.substitutions < aln.partialMatchSize){//tool2 outperformed
                            tool2TP_tool2Outperfomed++;
                            stats.tool2TP_tool2Outperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                            stats.tool2TP_tool2Outperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                            stats.baseLineTP_tool2Outperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                            stats.baseLineTP_tool2Outperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                            if (!cfg.comparisonSummaryOnly) cmp << "TP for " + tool2Name + ", " + tool2Name + " outperformed in terms of alnlen, " + tool2Name + " alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", baseLine alnlen: " << baselineAln.matches + baselineAln.inDels + baselineAln.substitutions << ", ed: " << baselineAln.substitutions + baselineAln.inDels << ", readID: " << baselineAln.readID << '\n';
                            stats.addTool2TPOutperform(aln.partialMatchSize - baselineAln.matches - baselineAln.inDels - baselineAln.substitutions);
                        } else {
                            if (baselineAln.substitutions + baselineAln.inDels < aln.substitutions + aln.inDels)//baseLine outperformed
                            {
                                tool2TP_baseLineOutperfomed++;
                                stats.tool2TP_baseLineOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                stats.tool2TP_baseLineOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                stats.baseLineTP_baseLineOutperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                stats.baseLineTP_baseLineOutperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                                if (!cfg.comparisonSummaryOnly) cmp << "TP for " + tool2Name + ", baseLine outperformed in terms of ed, " + tool2Name + " alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", baseLine alnlen: " << baselineAln.matches + baselineAln.inDels << ", ed: " << baselineAln.substitutions + baselineAln.inDels << ", readID: " << baselineAln.readID << '\n';
                                stats.addBaselineTPOutperform_edit(aln.substitutions + aln.inDels - baselineAln.substitutions - baselineAln.inDels);
                            } else if (baselineAln.substitutions + baselineAln.inDels > aln.substitutions + aln.inDels)//tool2 outperformed
                            {
                                tool2TP_tool2Outperfomed++;
                                stats.tool2TP_tool2Outperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                stats.tool2TP_tool2Outperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                stats.baseLineTP_tool2Outperfomed_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                stats.baseLineTP_tool2Outperfomed_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
                                if (!cfg.comparisonSummaryOnly) cmp << "TP for " + tool2Name + ", " + tool2Name + " outperformed in terms of ed, " + tool2Name + " alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", baseLine alnlen: " << baselineAln.matches + baselineAln.inDels << ", ed: " << baselineAln.substitutions + baselineAln.inDels << ", readID: " << baselineAln.readID << '\n';
                                stats.addTool2TPOutperform_edit(baselineAln.substitutions + baselineAln.inDels - aln.substitutions - aln.inDels);
                            } // tool2 and baseLine performed equally
                            else{
                                tool2TP_tool2EqualbaseLine++;
                                stats.tool2TP_tool2EqualbaseLine_alnLen_allQ.insert(aln.partialMatchSize);
                                stats.tool2TP_tool2EqualbaseLine_ed_allQ.insert(aln.substitutions + aln.inDels);
                                stats.baseLineTP_tool2EqualbaseLine_alnLen_allQ.insert(baselineAln.matches + baselineAln.inDels + baselineAln.substitutions);
                                stats.baseLineTP_tool2EqualbaseLine_ed_allQ.insert(baselineAln.substitutions + baselineAln.inDels);
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
            if (tool2TP > 0) stats.addTool2BestAlnSz(bestAlntool2.partialMatchSize);
            if(baseLineReadPerQuery > 0 && tool2ReadPerQuery > 0 && tool2TP == 0){//baseLine TP
                stats.tool2FN_noCriteriaFittedMatches++;
                stats.tool2FN_noCriteria_baseLine_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                stats.tool2FN_noCriteria_baseLine_ed_allQ.insert(baseLineBestAln.inDels + baseLineBestAln.substitutions);
                if (!cfg.comparisonSummaryOnly) cmp << tool2Name + "FN_noCriteriaFitted and baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine readID: " << baseLineBestAln.readID << '\n';
            }
            if(baseLineReadPerQuery == 0 && tool2TP > 0){//baseLine FN
                stats.baseLineFN_tool2_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                stats.baseLineFN_tool2_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                if (!cfg.comparisonSummaryOnly) cmp << "baseLine FN, and " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels << ", " + tool2Name + " readID: " << bestAlntool2.readID << '\n';
                stats.numnberOfBaseLine_FN++;
            } else if(baseLineReadPerQuery > 0 && tool2TP > 0){//compare the best alignmnets for this queyry
                if (baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions > bestAlntool2.partialMatchSize){ //baseLine outperformed
                    tool2Best_baseLineOutperfomed++;
                    stats.tool2Best_baseLineOutperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                    stats.tool2Best_baseLineOutperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                    stats.baseLineBest_baseLineOutperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                    stats.baseLineBest_baseLineOutperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                    if (!cfg.comparisonSummaryOnly) cmp << "Best: baseLine outperformed " + tool2Name + ", baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions <<  ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels <<  ", " + tool2Name + " readID: " << bestAlntool2.readID << '\n';
                    stats.addBaselineBestOutperform((baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions) - (bestAlntool2.partialMatchSize));
                } else if (baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions < bestAlntool2.partialMatchSize){//tool2 outperformed
                    tool2Best_tool2Outperfomed++;
                    stats.tool2Best_tool2Outperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                    stats.tool2Best_tool2Outperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                    stats.baseLineBest_tool2Outperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                    stats.baseLineBest_tool2Outperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                    if (!cfg.comparisonSummaryOnly) cmp << "Best: " + tool2Name + " outperformed baseLine in terms of alnlen, baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions <<  ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels << ", " + tool2Name + " readID: " << bestAlntool2.readID << '\n';
                    stats.addTool2BestOutperform((bestAlntool2.partialMatchSize) - (baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions));
                } else {
                    if (baseLineBestAln.substitutions + baseLineBestAln.inDels < bestAlntool2.substitutions + bestAlntool2.inDels)//baseLine outperformed
                    {
                        tool2Best_baseLineOutperfomed++;
                        stats.tool2Best_baseLineOutperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                        stats.tool2Best_baseLineOutperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                        stats.baseLineBest_baseLineOutperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                        stats.baseLineBest_baseLineOutperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                        if (!cfg.comparisonSummaryOnly) cmp << "Best: baseLine outperformed " + tool2Name + ", baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels <<  ", " + tool2Name + " readID: " << bestAlntool2.readID << '\n';
                        stats.addBaselineBestOutperform_edit((bestAlntool2.substitutions + bestAlntool2.inDels) - (baseLineBestAln.substitutions + baseLineBestAln.inDels));
                    } else if (baseLineBestAln.substitutions + baseLineBestAln.inDels > bestAlntool2.substitutions + bestAlntool2.inDels)//tool2 outperformed
                    {
                        tool2Best_tool2Outperfomed++;
                        stats.tool2Best_tool2Outperfomed_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                        stats.tool2Best_tool2Outperfomed_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                        stats.baseLineBest_tool2Outperfomed_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                        stats.baseLineBest_tool2Outperfomed_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                        if (!cfg.comparisonSummaryOnly) cmp << "Best: " + tool2Name + " outperformed baseLine in terms of ed, baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions <<  ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels <<  ", " + tool2Name + " readID: " << bestAlntool2.readID << '\n';
                        stats.addTool2BestOutperform_edit((baseLineBestAln.substitutions + baseLineBestAln.inDels) - (bestAlntool2.substitutions + bestAlntool2.inDels));
                    } // tool2 and baseLine performed equally
                    else{
                        tool2Best_tool2EqualbaseLine++;
                        stats.tool2Best_tool2EqualbaseLine_alnLen_allQ.insert(bestAlntool2.partialMatchSize);
                        stats.tool2Best_tool2EqualbaseLine_ed_allQ.insert(bestAlntool2.substitutions + bestAlntool2.inDels);
                        stats.baseLineBest_tool2EqualbaseLine_alnLen_allQ.insert(baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions);
                        stats.baseLineBest_tool2EqualbaseLine_ed_allQ.insert(baseLineBestAln.substitutions + baseLineBestAln.inDels);
                        if (!cfg.comparisonSummaryOnly) cmp << "Best: " + tool2Name + " & baseLine perfromed equally, baseLine alnsize: " << baseLineBestAln.matches + baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine ed: " << baseLineBestAln.inDels + baseLineBestAln.substitutions << ", baseLine readID: " << baseLineBestAln.readID << ", " + tool2Name + " alnsize: " << bestAlntool2.partialMatchSize << ", " + tool2Name + " ed: " << bestAlntool2.substitutions + bestAlntool2.inDels << ", " + tool2Name + " readID: " << bestAlntool2.readID << '\n';
                    }
                }
            }
        }
        // Run for every query, including queries with no tool candidates.
        for (const auto& entry : baselineReadID_Aln) {
            const auto& baseline = entry.second;
            const uint32_t baselineLength = baseline.matches + baseline.inDels + baseline.substitutions;
            const auto found = acceptedLengths.find(entry.first);
            const uint32_t toolLength = found == acceptedLengths.end() ? 0 : found->second;
            if (found == acceptedLengths.end()) {
                stats.addTool2FNAlnSz(baselineLength);
                ++baseLineTP_notool2Matches;
                stats.baseLineTP_notool2Matches_alnLen_allQ.insert(baselineLength);
                stats.baseLineTP_notool2Matches_ed_allQ.insert(baseline.substitutions + baseline.inDels);
            }
            if (toolLength < baselineLength) {
                ++stats.thresholdFNChanges[toolLength + 1];
                --stats.thresholdFNChanges[baselineLength + 1];
            }
        }
        // Include all evaluated queries, even when the tool has no candidates.
        stats.baseLineTP_notool2Matches_allQ += baseLineTP_notool2Matches;
        stats.tool2TPPerQuery.insert(tool2TP);
        stats.tool2FPPerQuery.insert(tool2FP);
        stats.totaltool2TP += tool2TP;
        stats.totaltool2FP += tool2FP;
        if (tool2TP > 0) ++stats.queriestool2FoundMatch;
        if (!cfg.comparisonSummaryOnly) {
            cmp << "baseLine TP that " + tool2Name + " did not find: " << baseLineTP_notool2Matches << '\n';
            cmp << tool2Name + " TP: " << tool2TP << "\t" << tool2Name + " FP: " << tool2FP << '\n';
        }
        stats.tool2TP_nobaseLineMatches_allQ += tool2TP_nobaseLineMatches;
        stats.tool2TP_tool2Outperfomed_allQ += tool2TP_tool2Outperfomed;
        stats.tool2TP_baseLineOutperfomed_allQ += tool2TP_baseLineOutperfomed;
        stats.tool2TP_tool2EqualbaseLine_allQ += tool2TP_tool2EqualbaseLine;
        if (!cfg.comparisonSummaryOnly) cmp << tool2Name + "TP_nobaseLineMatches: " << tool2TP_nobaseLineMatches << ", " + tool2Name + "TP_tool2Outperfomed: " << tool2TP_tool2Outperfomed << ", " + tool2Name + "TP_baseLineOutperfomed: " << tool2TP_baseLineOutperfomed << ", " + tool2Name + "TP_tool2EqualbaseLine: " << tool2TP_tool2EqualbaseLine << '\n';
        stats.tool2NoExactRun_allQ += tool2NoExactRun;
        stats.tool2FP_lowPercentageIdentity_allQ += tool2FP_lowPercentageIdentity;
        if (!cfg.comparisonSummaryOnly) cmp << tool2Name + "NoExactRun (diagnostic): " << tool2NoExactRun << ", " + tool2Name << "FP_lowPercentageIdentity: " << tool2FP_lowPercentageIdentity << '\n';
        stats.tool2Best_baseLineOutperfomed_allQ += tool2Best_baseLineOutperfomed;
        stats.tool2Best_tool2Outperfomed_allQ += tool2Best_tool2Outperfomed;
        stats.tool2Best_tool2EqualbaseLine_allQ += tool2Best_tool2EqualbaseLine;
        if (!cfg.comparisonSummaryOnly) cmp << tool2Name + "Best_baseLineOutperfomed: " << tool2Best_baseLineOutperfomed << ", " + tool2Name + "Best_tool2Outperfomed: " << tool2Best_tool2Outperfomed << ", " + tool2Name + "Best_tool2EqualbaseLine: " << tool2Best_tool2EqualbaseLine << '\n';
        // best_aln_sz.push_back(make_pair(baseLineBestAlnSize, tool2BestAlnSize));
        //baseLine and tool2 alignments per query
        if(REPORT_ALN_PER_Q) alnPerQ << queryInd << " " << baselineReadID_Aln.size() << " " << tool2TP << '\n';
    }

public:
    void compareWithBaseLine(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress,
    const uint32_t queryCount, const string alnReportAddressBase, string baseLineFilePrefixAddress)
    {
        if (cfg.comparisonThreads == 0)
            throw invalid_argument("comparisonThreads must be positive");
        // Chunk naming and reader bounds require consecutive zero-based IDs.
        // Reject incompatible inputs before replacing any comparison reports.
        if (queries.size() != queryCount)
            throw invalid_argument("Query count does not match the number of unique query IDs");
        for (uint32_t q = 0; q < queryCount; ++q) {
            const auto query = queries.find(q);
            if (query == queries.end() || query->second.empty())
                throw invalid_argument("Comparison requires nonempty queries with IDs 0..queryCount-1; missing or empty ID " + to_string(q));
        }
        const string tool2Name = cfg.otherTool;
        cout << "Comparing Baseline and " + tool2Name << endl;
        // Validate tool input before creating comparison reports.
        {
            ifstream input(cfg.otherToolOutputFileAddress);
            if (!input)
                throw runtime_error("Cannot open comparison alignment input: " +
                                    cfg.otherToolOutputFileAddress);
        }
        vector<Alignment> tool2Alignments;
        bool isBwotie = (tool2Name == "bowtie" || tool2Name == "BOWTIE");
        bool isParmik = (tool2Name == "parmik" || tool2Name == "PARMIK");
        bool isMinimap2 = (tool2Name == "minimap2" || tool2Name == "MINIMAP2");
        if (tool2Name == "BLAST" || tool2Name == "blast") {
            if (queryCount > 0)
                BlastReader(cfg.otherToolOutputFileAddress).parseFile(queryCount - 1, tool2Alignments);
        } else if (isMinimap2) {
            Minimap2Reader(cfg.otherToolOutputFileAddress).parseFile(queryCount, tool2Alignments);
        } else if (isParmik || isBwotie) {
            if (queryCount > 0)
                SamReader(cfg.otherToolOutputFileAddress).parseFile(queryCount - 1, tool2Alignments, false, isBwotie);
        } else {
            throw runtime_error("Unsupported baseline comparison tool: " + tool2Name);
        }
        ofstream cmp(comparisonResultsFileAddress);
        cout << "cmp file address: " << comparisonResultsFileAddress << '\n';
        if(!cmp.is_open())
            throw runtime_error("Cannot create comparison report: " + comparisonResultsFileAddress);
        ofstream alnPerQ;
        ofstream baseLineFn;
        ofstream fnAlnSz;
        // A missed pair contributes FN for toolLength < R <= baselineLength.
        // Signed range deltas retain that relationship without storing all pairs.
        ofstream thresholdRecall(alnReportAddressBase + "RecallByThreshold.txt");
        if (!thresholdRecall)
            throw runtime_error("Cannot create threshold recall report");
        ofstream baselineRecall(alnReportAddressBase + "BaselineRecallByThreshold.txt");
        ofstream additionalHits(alnReportAddressBase + "AdditionalHitsByThreshold.txt");
        if (!baselineRecall || !additionalHits)
            throw runtime_error("Cannot create baseline recall reports: " + alnReportAddressBase);
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
        if (!bestAlnCmp || !tpAlnCmp || !bestEdCmp || !tpEdCmp ||
            !tpAlnSz || !bestAlnSz || !fnAlnSz)
            throw runtime_error("Cannot create comparison count reports: " + alnReportAddressBase);
        IndexContainer<uint32_t, Alignment> queryTool2Alignment;
        for (const Alignment& aln : tool2Alignments)
        {
            queryTool2Alignment.put(aln.queryID, aln);
        }
        tool2Alignments.clear();


        const unsigned workerCount = std::min<unsigned>(cfg.comparisonThreads,
                                                        std::max<uint32_t>(1, std::min<uint32_t>(1000, queryCount)));
        vector<Statistics> workerStats(workerCount);
        vector<Alignmen_Baseline> baseLineSamAlignments;
        IndexContainer<uint32_t, Alignmen_Baseline> queryBaselineAlignment;
        // One baseline chunk stays resident. Batches keep detailed log buffers
        // bounded to a few queries per worker, instead of the entire run.
        const uint32_t batchSize = workerCount * 4;
        for (uint32_t chunkBase = 0; chunkBase < queryCount;) {
            const uint32_t chunkEnd = chunkBase + std::min<uint32_t>(1000, queryCount - chunkBase);
            queryBaselineAlignment.clear();
            cout << "queries processed: " << chunkBase << " / " << queryCount << endl;
            readChunkOfBaseLine(baseLineSamAlignments, baseLineFilePrefixAddress, chunkBase, 1000);
            for (const auto& aln : baseLineSamAlignments)
                queryBaselineAlignment.put(aln.queryID, aln);
            baseLineSamAlignments.clear();

            for (uint32_t batchBase = chunkBase; batchBase < chunkEnd;) {
                const uint32_t batchCount = std::min(batchSize, chunkEnd - batchBase);
                struct QueryOutput {
                    ostringstream comparison, perQuery, baselineFN;
                    exception_ptr error;
                };
                vector<QueryOutput> output(batchCount);
                #pragma omp parallel for schedule(dynamic, 1) num_threads(workerCount)
                for (uint32_t offset = 0; offset < batchCount; ++offset) {
                    auto& result = output[offset];
                    try {
                        compareQuery(cfg, batchBase + offset, queries, queryTool2Alignment,
                                     queryBaselineAlignment, workerStats[omp_get_thread_num()],
                                     result.comparison, result.perQuery, result.baselineFN);
                    } catch (...) {
                        // Rethrow outside OpenMP on the caller's thread.
                        result.error = current_exception();
                    }
                }
                for (auto& result : output) {
                    if (result.error) rethrow_exception(result.error);
                    cmp << result.comparison.str();
                    if (REPORT_ALN_PER_Q) alnPerQ << result.perQuery.str();
                    if (REPORT_baseLine_FN) baseLineFn << result.baselineFN.str();
                }
                batchBase += batchCount;
            }
            chunkBase = chunkEnd;
        }
        Statistics stats;
        for (auto& local : workerStats) stats.merge(local);
        Utilities<uint32_t> util;
        cmp << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Overall Comparison Results>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << '\n';
        cmp << "BaselineRecallByThreshold reports conventional baseline-pair recall.\n";
        cmp << "TP/FP summaries below count accepted/rejected tool pairs across all lengths.\n";
        cmp << left << setw(80) << "# Of queries : " + to_string(queryCount) << '\n';
        cmp << left << setw(80) << "# Of evaluated queries (excluding N): " << queryCount - stats.numberOfQueryContainN << '\n';
        cmp << left << setw(80) << "# Of queries that baseLine found match : " << stats.queriesblfoundMatch << '\n';
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " found match (TP) : " << stats.queriestool2FoundMatch << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of queries that Baseline didn't find TN : " << stats.totalbaseLineTN << '\n';
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find TN : " << stats.totalTool2TN << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find FN (Total) : " << stats.tool2FN + stats.tool2FN_noCriteriaFittedMatches << '\n';
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find any match FN : " << stats.tool2FN << '\n';
        pair<uint32_t, uint32_t> avgTool2FN_baseLine_alnLen_allQ = util.calculateStatistics2(stats.tool2FN_baseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgTool2FN_baseLine_ed_allQ = util.calculateStatistics2(stats.tool2FN_baseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) Baseline's alignment length where " + tool2Name + " FN: " << avgTool2FN_baseLine_alnLen_allQ.first << ", " << avgTool2FN_baseLine_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) Baseline's edit distance where " + tool2Name + " FN: " << avgTool2FN_baseLine_ed_allQ.first << ", " << avgTool2FN_baseLine_ed_allQ.second << '\n';
        cmp << left << setw(80) << "# Of queries that " + tool2Name + " didn't find based on criteria FN : " << stats.tool2FN_noCriteriaFittedMatches << '\n';
        pair<uint32_t, uint32_t> avgtool2FN_noCriteria_baseLine_alnLen_allQ = util.calculateStatistics2(stats.tool2FN_noCriteria_baseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2FN_noCriteria_baseLine_ed_allQ = util.calculateStatistics2(stats.tool2FN_noCriteria_baseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where " + tool2Name + " FN based on criteria: " << avgtool2FN_noCriteria_baseLine_alnLen_allQ.first << ", " << avgtool2FN_noCriteria_baseLine_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where " + tool2Name + " FN based on criteria: " << avgtool2FN_noCriteria_baseLine_ed_allQ.first << ", " << avgtool2FN_noCriteria_baseLine_ed_allQ.second << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of queries that baseLine didn't find FN : " << stats.numnberOfBaseLine_FN << '\n';
        pair<uint32_t, uint32_t> avgBaselineFN_tool2_alnLen_allQ = util.calculateStatistics2(stats.baseLineFN_tool2_alnLen_allQ);
        pair<uint32_t, uint32_t> avgBaselineFN_tool2_ed_allQ = util.calculateStatistics2(stats.baseLineFN_tool2_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where baseLine FN: " << avgBaselineFN_tool2_alnLen_allQ.first << ", " << avgBaselineFN_tool2_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where baseLine FN: " << avgBaselineFN_tool2_ed_allQ.first << ", " << avgBaselineFN_tool2_ed_allQ.second << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of readIDs found by Baseline (total): " << stats.baseLineTotalNumberOfReadIDs << '\n';
        cmp << left << setw(80) << "# Of readIDs found by " + tool2Name + " (total): " << stats.tool2TotalNumberOfReadIDs << '\n';
        // cmp << left << setw(80) << "# Of readIDs found by " + tool2Name + " containing N (total): " << numberOfTool2ReadsContainingN << endl;

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of total " + tool2Name + " (TP): " << stats.totaltool2TP << '\n';
        pair<uint32_t, uint32_t> avgtool2TPPerQuery = util.calculateStatistics2(stats.tool2TPPerQuery);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s TP per query: " << avgtool2TPPerQuery.first << ", " << avgtool2TPPerQuery.second << '\n';

        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where baseLine didn't find match: " << stats.tool2TP_nobaseLineMatches_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2TP_nobaseLineMatches_alnLen_allQ = util.calculateStatistics2(stats.tool2TP_nobaseLineMatches_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2TP_nobaseLineMatches_ed_allQ = util.calculateStatistics2(stats.tool2TP_nobaseLineMatches_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where baseLine didn't find match: " << avgtool2TP_nobaseLineMatches_alnLen_allQ.first << ", " << avgtool2TP_nobaseLineMatches_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where baseLine didn't find match: " << avgtool2TP_nobaseLineMatches_ed_allQ.first << ", " << avgtool2TP_nobaseLineMatches_ed_allQ.second << '\n';

        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where baseLine outperformed: " << stats.tool2TP_baseLineOutperfomed_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2TP_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(stats.tool2TP_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2TP_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(stats.tool2TP_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where baseLine outperformed: " << avgtool2TP_baseLineOutperfomed_alnLen_allQ.first << ", " << avgtool2TP_baseLineOutperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where baseLine outperformed: " << avgtool2TP_baseLineOutperfomed_ed_allQ.first << ", " << avgtool2TP_baseLineOutperfomed_ed_allQ.second << '\n';
        pair<uint32_t, uint32_t> avgbaseLineTP_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(stats.baseLineTP_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineTP_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(stats.baseLineTP_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where baseLine outperformed: " << avgbaseLineTP_baseLineOutperfomed_alnLen_allQ.first << ", " << avgbaseLineTP_baseLineOutperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where baseLine outperformed: " << avgbaseLineTP_baseLineOutperfomed_ed_allQ.first << ", " << avgbaseLineTP_baseLineOutperfomed_ed_allQ.second << '\n';

        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where " + tool2Name + " outperformed: " << stats.tool2TP_tool2Outperfomed_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2TP_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(stats.tool2TP_tool2Outperfomed_alnLen_allQ);

        pair<uint32_t, uint32_t> avgtool2TP_tool2Outperfomed_ed_allQ = util.calculateStatistics2(stats.tool2TP_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where " + tool2Name + " outperformed: " << avgtool2TP_tool2Outperfomed_alnLen_allQ.first << ", " << avgtool2TP_tool2Outperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where " + tool2Name + " outperformed: " << avgtool2TP_tool2Outperfomed_ed_allQ.first << ", " << avgtool2TP_tool2Outperfomed_ed_allQ.second << '\n';
        pair<uint32_t, uint32_t> avgbaseLineTP_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(stats.baseLineTP_tool2Outperfomed_alnLen_allQ);

        pair<uint32_t, uint32_t> avgbaseLineTP_tool2Outperfomed_ed_allQ = util.calculateStatistics2(stats.baseLineTP_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where " + tool2Name + " outperformed: " << avgbaseLineTP_tool2Outperfomed_alnLen_allQ.first << ", " << avgbaseLineTP_tool2Outperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where " + tool2Name + " outperformed: " << avgbaseLineTP_tool2Outperfomed_ed_allQ.first << ", " << avgbaseLineTP_tool2Outperfomed_ed_allQ.second << '\n';

        cmp << left << setw(80) << "# Of " + tool2Name + " (TP) where they perfromed equaly: " << stats.tool2TP_tool2EqualbaseLine_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2TP_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(stats.tool2TP_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2TP_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(stats.tool2TP_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where they perfromed equaly: " << avgtool2TP_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgtool2TP_tool2EqualbaseLine_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where they perfromed equaly: " << avgtool2TP_tool2EqualbaseLine_ed_allQ.first << ", " << avgtool2TP_tool2EqualbaseLine_ed_allQ.second << '\n';
        pair<uint32_t, uint32_t> avgbaseLineTP_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(stats.baseLineTP_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineTP_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(stats.baseLineTP_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where they perfromed equaly: " << avgbaseLineTP_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgbaseLineTP_tool2EqualbaseLine_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where they perfromed equaly: " << avgbaseLineTP_tool2EqualbaseLine_ed_allQ.first << ", " << avgbaseLineTP_tool2EqualbaseLine_ed_allQ.second << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of " + tool2Name + " best where baseLine outperformed: " << stats.tool2Best_baseLineOutperfomed_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2Best_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(stats.tool2Best_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2Best_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(stats.tool2Best_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best alignment length where baseLine outperformed: " << avgtool2Best_baseLineOutperfomed_alnLen_allQ.first << ", " << avgtool2Best_baseLineOutperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best edit distance where baseLine outperformed: " << avgtool2Best_baseLineOutperfomed_ed_allQ.first << ", " << avgtool2Best_baseLineOutperfomed_ed_allQ.second << '\n';
        pair<uint32_t, uint32_t> avgbaseLineBest_baseLineOutperfomed_alnLen_allQ = util.calculateStatistics2(stats.baseLineBest_baseLineOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineBest_baseLineOutperfomed_ed_allQ = util.calculateStatistics2(stats.baseLineBest_baseLineOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's best alignment length where baseLine outperformed: " << avgbaseLineBest_baseLineOutperfomed_alnLen_allQ.first << ", " << avgbaseLineBest_baseLineOutperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's best edit distance where baseLine outperformed: " << avgbaseLineBest_baseLineOutperfomed_ed_allQ.first << ", " << avgbaseLineBest_baseLineOutperfomed_ed_allQ.second << '\n';

        cmp << left << setw(80) << "# Of " + tool2Name + " best where " + tool2Name + " outperformed: " << stats.tool2Best_tool2Outperfomed_allQ  << '\n';
        pair<uint32_t, uint32_t> avgtool2Best_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(stats.tool2Best_tool2Outperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2Best_tool2Outperfomed_ed_allQ = util.calculateStatistics2(stats.tool2Best_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best alignment length where " + tool2Name + " outperformed: " << avgtool2Best_tool2Outperfomed_alnLen_allQ.first << ", " << avgtool2Best_tool2Outperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best edit distance where " + tool2Name + " outperformed: " << avgtool2Best_tool2Outperfomed_ed_allQ.first << ", " << avgtool2Best_tool2Outperfomed_ed_allQ.second << '\n';
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2Outperfomed_alnLen_allQ = util.calculateStatistics2(stats.baseLineBest_tool2Outperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2Outperfomed_ed_allQ = util.calculateStatistics2(stats.baseLineBest_tool2Outperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's best alignment length where " + tool2Name + " outperformed: " << avgbaseLineBest_tool2Outperfomed_alnLen_allQ.first << ", " << avgbaseLineBest_tool2Outperfomed_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's best edit distance where " + tool2Name + " outperformed: " << avgbaseLineBest_tool2Outperfomed_ed_allQ.first << ", " << avgbaseLineBest_tool2Outperfomed_ed_allQ.second << '\n';

        cmp << left << setw(80) << "# Of " + tool2Name + " best where they perfromed equaly: " << stats.tool2Best_tool2EqualbaseLine_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2Best_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(stats.tool2Best_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2Best_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(stats.tool2Best_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best alignment length where they perfromed equaly: " << avgtool2Best_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgtool2Best_tool2EqualbaseLine_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s best edit distance where they perfromed equaly: " << avgtool2Best_tool2EqualbaseLine_ed_allQ.first << ", " << avgtool2Best_tool2EqualbaseLine_ed_allQ.second << '\n';
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2EqualbaseLine_alnLen_allQ = util.calculateStatistics2(stats.baseLineBest_tool2EqualbaseLine_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineBest_tool2EqualbaseLine_ed_allQ = util.calculateStatistics2(stats.baseLineBest_tool2EqualbaseLine_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's best alignment length where they perfromed equaly: " << avgbaseLineBest_tool2EqualbaseLine_alnLen_allQ.first << ", " << avgbaseLineBest_tool2EqualbaseLine_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's best edit distance where they perfromed equaly: " << avgbaseLineBest_tool2EqualbaseLine_ed_allQ.first << ", " << avgbaseLineBest_tool2EqualbaseLine_ed_allQ.second << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of total " + tool2Name + " (FP): " << stats.totaltool2FP << '\n';
        pair<uint32_t, uint32_t> avgtool2FPPerQuery = util.calculateStatistics2(stats.tool2FPPerQuery);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s FP per query: " << avgtool2FPPerQuery.first << ", " << avgtool2FPPerQuery.second << '\n';

        cmp << left << setw(80) << "# Of " + tool2Name + " candidates with no exact match run of length " + to_string(cfg.kmerLength) + " (diagnostic): " << stats.tool2NoExactRun_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2NoExactRun_alnLen_allQ = util.calculateStatistics2(stats.tool2NoExactRun_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2NoExactRun_ed_allQ = util.calculateStatistics2(stats.tool2NoExactRun_ed_allQ);
        pair<uint32_t, uint32_t> avgtool2NoExactRun_matchsize_allQ = util.calculateStatistics2(stats.tool2NoExactRun_matchsize_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where no exact match run meets k (diagnostic): " << avgtool2NoExactRun_alnLen_allQ.first << ", " << avgtool2NoExactRun_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where no exact match run meets k (diagnostic): " << avgtool2NoExactRun_ed_allQ.first << ", " << avgtool2NoExactRun_ed_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s longest exact run where no exact match run meets k (diagnostic): " << avgtool2NoExactRun_matchsize_allQ.first << ", " << avgtool2NoExactRun_matchsize_allQ.second << '\n';

        // cmp << left << setw(80) << "# Of " + tool2Name + " (FP) that # of kmers is lower than regionSize" << tool2FP_lowerExactMatchKmers_allQ << endl;
        // pair<uint32_t, uint32_t> avgtool2FP_lowerExactMatchKmers_alnLen_allQ = util.calculateStatistics2(tool2FP_lowerExactMatchKmers_alnLen_allQ);
        // pair<uint32_t, uint32_t> avgtool2FP_lowerExactMatchKmers_ed_allQ = util.calculateStatistics2(tool2FP_lowerExactMatchKmers_ed_allQ);
        // cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where # of kmers < MinExactNumberOfKmers" <<  avgtool2FP_lowerExactMatchKmers_alnLen_allQ.first << ", " << avgtool2FP_lowerExactMatchKmers_alnLen_allQ.second << endl;
        // cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where # of kmers < MinExactNumberOfKmers" << avgtool2FP_lowerExactMatchKmers_ed_allQ.first << ", " << avgtool2FP_lowerExactMatchKmers_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of " + tool2Name + " (FP) that percentage identity is lower than PI" << stats.tool2FP_lowPercentageIdentity_allQ << '\n';
        pair<uint32_t, uint32_t> avgtool2FP_lowPercentageIdentity_alnLen_allQ = util.calculateStatistics2(stats.tool2FP_lowPercentageIdentity_alnLen_allQ);
        pair<uint32_t, uint32_t> avgtool2FP_lowPercentageIdentity_ed_allQ = util.calculateStatistics2(stats.tool2FP_lowPercentageIdentity_ed_allQ);
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s alignment length where percentage identity is lower than PI" <<  avgtool2FP_lowPercentageIdentity_alnLen_allQ.first << ", " << avgtool2FP_lowPercentageIdentity_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s edit distance where percentage identity is lower than PI" << avgtool2FP_lowPercentageIdentity_ed_allQ.first << ", " << avgtool2FP_lowPercentageIdentity_ed_allQ.second << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of baseLine (TP) where " + tool2Name + " didn't find match: " << stats.baseLineTP_notool2Matches_allQ << '\n';
        pair<uint32_t, uint32_t> avgbaseLineTP_notool2Matches_alnLen_allQ = util.calculateStatistics2(stats.baseLineTP_notool2Matches_alnLen_allQ);
        pair<uint32_t, uint32_t> avgbaseLineTP_notool2Matches_ed_allQ = util.calculateStatistics2(stats.baseLineTP_notool2Matches_ed_allQ);
        cmp << left << setw(80) << "(avg, median) baseLine's alignment length where " + tool2Name + " didn't find match: " << avgbaseLineTP_notool2Matches_alnLen_allQ.first << ", " << avgbaseLineTP_notool2Matches_alnLen_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) baseLine's edit distance where " + tool2Name + " didn't find match: " << avgbaseLineTP_notool2Matches_ed_allQ.first << ", " << avgbaseLineTP_notool2Matches_ed_allQ.second << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        cmp << left << setw(80) << "# Of queries contains N: " << stats.numberOfQueryContainN << '\n';

        cmp << "------------------------------------------------------------------------------------------" << '\n';
        pair<uint32_t, uint32_t> baseLineReadPerQuery_allQ  = util.calculateStatistics2(stats.baseLineReadPerQuerySet);
        pair<uint32_t, uint32_t> tool2ReadPerQuery_allQ  = util.calculateStatistics2(stats.tool2ReadPerQuerySet);
        cmp << left << setw(80) << "(avg, median) baseLineS's Matches per Query: " << baseLineReadPerQuery_allQ.first << ", " << baseLineReadPerQuery_allQ.second << '\n';
        cmp << left << setw(80) << "(avg, median) " + tool2Name + "'s Matches per Query: " << tool2ReadPerQuery_allQ.first << ", " << tool2ReadPerQuery_allQ.second << '\n';
        cmp << "------------------------------------------------------------------------------------------" << '\n';

        for(uint32_t i = 1; i < cfg.contigSize; i++) {
            tpAlnCmp << i << "\t" << stats.baselineTPOutperform_aln_sz_difference[i] << "\t" << stats.tool2TPOutperform_aln_sz_difference[i] << '\n';
            bestAlnCmp << i << "\t" << stats.baselineBestOutperform_aln_sz_difference[i] << "\t" << stats.tool2BestOutperform_aln_sz_difference[i] << '\n';
            tpEdCmp << i << "\t" << stats.baselineTPOutperform_edit_difference[i] << "\t" << stats.tool2TPOutperform_edit_difference[i] << '\n';
            bestEdCmp << i << "\t" << stats.baselineBestOutperform_edit_difference[i] << "\t" << stats.tool2BestOutperform_edit_difference[i] << '\n';
        }

        for(uint32_t i = 1; i < cfg.contigSize + 20; i++) {
            tpAlnSz << i << "\t" << stats.baselineTP_aln_sz[i] << "\t" << stats.tool2TP_aln_sz[i] << '\n';
            bestAlnSz << i << "\t" << stats.baselineBest_aln_sz[i] << "\t" << stats.tool2Best_aln_sz[i] << '\n';
            fnAlnSz << i << "\t" << stats.baselineFN_aln_sz[i] << "\t" << stats.tool2FN_aln_sz[i] << '\n';
        }
        // Preserve all accepted tool hits in the legacy report; restrict
        // conventional recall to qualifying baseline pairs at each threshold.
        uint64_t thresholdTP = 0, thresholdBaseline = 0;
        for (const auto& bin : stats.tool2TP_aln_sz) thresholdTP += bin.second;
        for (const auto& bin : stats.baselineTP_aln_sz) thresholdBaseline += bin.second;
        int64_t thresholdFN = 0;
        const uint32_t maxThreshold = std::max<uint32_t>(255, cfg.contigSize + 19);
        for (uint32_t r = 1; r <= maxThreshold; ++r) {
            thresholdFN += stats.thresholdFNChanges[r];
            if (thresholdFN < 0 || uint64_t(thresholdFN) > thresholdBaseline)
                throw logic_error("Inconsistent baseline false negatives at threshold " + to_string(r));
            const uint64_t recoveredBaseline = thresholdBaseline - uint64_t(thresholdFN);
            if (recoveredBaseline > thresholdTP)
                throw logic_error("Recovered baseline pairs exceed accepted tool pairs at threshold " + to_string(r));
            thresholdRecall << r << "\t" << thresholdTP << "\t" << thresholdFN << "\n";
            baselineRecall << r << "\t" << recoveredBaseline << "\t" << thresholdFN << "\n";
            additionalHits << r << "\t" << thresholdTP - recoveredBaseline << "\n";
            thresholdTP -= stats.tool2TP_aln_sz[r];
            thresholdBaseline -= stats.baselineTP_aln_sz[r];
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
        thresholdRecall.close();
        baselineRecall.close();
        additionalHits.close();
        if (!cmp || !bestAlnCmp || !tpAlnCmp || !bestEdCmp || !tpEdCmp ||
            !tpAlnSz || !bestAlnSz || !fnAlnSz || !thresholdRecall || !baselineRecall || !additionalHits)
            throw runtime_error("Failed writing comparison reports: " + alnReportAddressBase);
        if (REPORT_baseLine_FN) baseLineFn.close();
        cout << "queries processed: " << queryCount << " / " << queryCount << endl;
    }
};

#endif
