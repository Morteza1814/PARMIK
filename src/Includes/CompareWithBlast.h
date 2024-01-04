#ifndef COMAPREWITHBLAST_H
#define COMAPREWITHBLAST_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include "BlastReader.h"

using namespace std;

class CompareWithBlast {
public: 

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

    bool hasMinConsecutiveMatches(uint32_t queryS, const string& str1, const string& str2, Config &cfg) {
        uint32_t consecutiveMatchCount = 0;

        // Ensure that both strings have the same length for exact matching
        if (str1.length() != str2.length()) {
            cerr << "Error: Sequences have different lengths." << endl;
            return false;
        }

        // Iterate through each character in the strings and count consecutive exact matches
        uint32_t startPos = 0;
        for (size_t i = 0; i < str1.length(); i++) {
            if (str1[i] == str2[i]) {
                consecutiveMatchCount++;
                if (consecutiveMatchCount >= cfg.minExactMatchLen) {
                    if ((startPos + queryS + cfg.minExactMatchLen <= cfg.regionSize) || ((startPos + queryS >= cfg.contigSize - cfg.regionSize) && (startPos + queryS + cfg.minExactMatchLen <= cfg.contigSize))) {
                        return true;  // Found enough consecutive matches inth front or back region
                    }
                }
            } else {
                consecutiveMatchCount = 0;  // Reset count if consecutive match is broken
                startPos = i+1;
            }
        }

        // Return false if the number of consecutive exact matches is less than minConsecutiveMatch
        return false;
    }

    bool checkBlastEditPositions(BlastReader::Blast& blastAlignment, Config cfg)
    {
        // vector<int> editPos;
        // determineEditsLocationsAndType(blastAlignment.readAligned, blastAlignment.queryAligned, editPos);
        // if (editPos.size() != blastAlignment.Mismatches) 
        // {
        //     cout << "Error: Edit distance not equals mismatch." << endl;
        //     return false;
        // }
        // bool firstKmerInFrontRegionDismissed = false, lastKmerInFrontRegionDismissed = false,
        //     firstKmerInBackRegionDismissed = false, lastKmerInBacktRegionDismissed = false;
        bool frontRegionDismissed = false, backRegionDismissed = false;
        uint32_t queryS = blastAlignment.queryS - 1;
        // for(auto v : editPos)
        // {
        //     if(((v + queryS) >= 0 && (v + queryS) <= cfg.minExactMatchLen))
        //     {
        //         firstKmerInFrontRegionDismissed = true;
        //     }
        //     if(((v + queryS) >= (cfg.regionSize - cfg.minExactMatchLen) && (v + queryS) <= cfg.regionSize))
        //     {
        //         lastKmerInFrontRegionDismissed = true;
        //     }
        //     if(((v + queryS) >= (cfg.contigSize - cfg.regionSize)) && ((v + queryS) <= (cfg.contigSize - cfg.regionSize + cfg.minExactMatchLen)))
        //     {
        //         firstKmerInBackRegionDismissed = true;
        //     }
        //     if(((v + queryS) >= (cfg.contigSize - cfg.minExactMatchLen)) && ((v + queryS) <= cfg.contigSize))
        //     {
        //         lastKmerInBacktRegionDismissed = true;
        //     }
        // }
        //check whether the region starts at most from 100 + 2
        if(queryS > (cfg.contigSize - cfg.regionSize + cfg.editDistance)) //first bp of back region starts from 100
            return false;
        //check whether the regions has at least minExactMatchLen consecutive matches
        if(!hasMinConsecutiveMatches(queryS, blastAlignment.queryAligned, blastAlignment.readAligned, cfg))
            return false;
        // if(firstKmerInFrontRegionDismissed && lastKmerInFrontRegionDismissed)
        //     frontRegionDismissed = true;
        // if(firstKmerInBackRegionDismissed && lastKmerInBacktRegionDismissed)
        //     backRegionDismissed = true;
        if(cfg.editDistance >= queryS){
            auto editsAllwedInFrontRegion = cfg.editDistance - queryS;
            if(blastAlignment.Mismatches + blastAlignment.InDels > editsAllwedInFrontRegion)
                frontRegionDismissed = true;
        } else {
            frontRegionDismissed = true;
        }
        if(queryS + blastAlignment.AlignmentLength >= cfg.contigSize - cfg.editDistance){ //query start position + alignment len should be larger than conig size - allowed edit
            auto editsAllwedInBackRegion = queryS + blastAlignment.AlignmentLength - (cfg.contigSize - cfg.editDistance);
            if(blastAlignment.Mismatches + blastAlignment.InDels > editsAllwedInBackRegion)
                backRegionDismissed = true;
        }else {
            backRegionDismissed = true;
        }
        if(frontRegionDismissed && backRegionDismissed)
            return false;
        return true;
    }

    void comparePmWithBlast(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress, IndexContainer<uint32_t, LevAlign>& pmAlignments, vector<pair<uint32_t, uint32_t>>& alnPmBLASTHisto, const uint32_t queryCount, string alnPerQueryFileAddress)
    {
        ofstream cmp(comparisonResultsFileAddress);
        ofstream alnPerQ(alnPerQueryFileAddress);
        SamReader sam(cfg.otherToolOutputFileAddress);
        BlastReader blastReader(cfg.otherToolOutputFileAddress);
        IndexContainer<uint32_t, BlastReader::Blast> blastAlignments;
        blastReader.parseFile(queryCount, blastAlignments);
        uint32_t numberOfQueryContainN = 0, numberOfBlastReadsContainingN = 0, 
        queriesBLASTFoundMatch = 0, queriesPMFoundMatch = 0,
        blastTotalNumberOfReadIDs = 0, pmTotalNumberOfReadIDs = 0,
        numnberOfPM_FN = 0, blastFN = 0, blastFN_noCriteriaFittedMatches = 0, // places where BLAST found no matches or none of the matches conform with the criteria
        totalBlastTN = 0,
        totalBlastTP = 0, 
        totalBlastFP = 0,
        blastTP_noParmikMatches_allQ = 0, blastTP_blastOutperfomed_allQ = 0, blastTP_parmikOutperfomed_allQ = 0, blastTP_blastEqualParmik_allQ = 0,
        parmikTP_noBlastMatches_allQ = 0,
        blastFP_editsExceed_allQ = 0, blastFP_lowAlnLen_allQ = 0, blastFP_editPos_allQ = 0,
        blastBest_blastOutperfomed_allQ = 0, blastBest_parmikOutperfomed_allQ = 0, blastBest_blastEqualParmik_allQ = 0;
        //all query level parameters sets (if total is needed, just add all the elements in the set)
        set<uint32_t> pmReadPerQuerySet, blastReadPerQuerySet;
        set<uint32_t> blastTPPerQuery, blastFPPerQuery;
        set<uint32_t> blastTP_noParmikMatches_alnLen_allQ, blastTP_blastOutperfomed_alnLen_allQ, blastTP_parmikOutperfomed_alnLen_allQ, blastTP_blastEqualParmik_alnLen_allQ; 
        set<uint32_t> blastTP_noParmikMatches_ed_allQ, blastTP_blastOutperfomed_ed_allQ, blastTP_parmikOutperfomed_ed_allQ, blastTP_blastEqualParmik_ed_allQ; 
        set<uint32_t> parmikTP_noBlastMatches_alnLen_allQ, parmikTP_noBlastMatches_ed_allQ;
        set<uint32_t> parmikTP_blastOutperfomed_alnLen_allQ, parmikTP_parmikOutperfomed_alnLen_allQ, parmikTP_blastEqualParmik_alnLen_allQ; 
        set<uint32_t> parmikTP_blastOutperfomed_ed_allQ, parmikTP_parmikOutperfomed_ed_allQ, parmikTP_blastEqualParmik_ed_allQ; 
        // best alignment comparisons
        set<uint32_t> blastBest_blastOutperfomed_alnLen_allQ, blastBest_parmikOutperfomed_alnLen_allQ, blastBest_blastEqualParmik_alnLen_allQ; 
        set<uint32_t> blastBest_blastOutperfomed_ed_allQ, blastBest_parmikOutperfomed_ed_allQ, blastBest_blastEqualParmik_ed_allQ; 
        set<uint32_t> parmikBest_blastOutperfomed_alnLen_allQ, parmikBest_parmikOutperfomed_alnLen_allQ, parmikBest_blastEqualParmik_alnLen_allQ; 
        set<uint32_t> parmikBest_blastOutperfomed_ed_allQ, parmikBest_parmikOutperfomed_ed_allQ, parmikBest_blastEqualParmik_ed_allQ; 
        set<uint32_t> blastFP_editsExceed_alnLen_allQ, blastFP_lowAlnLen_alnLen_allQ, blastFP_editPos_alnLen_allQ;
        set<uint32_t> blastFP_editsExceed_ed_allQ, blastFP_lowAlnLen_ed_allQ, blastFP_editPos_ed_allQ;
        set<uint32_t> blastFN_parmik_alnLen_allQ, blastFN_parmik_ed_allQ;// PARMIK alignment characteristics for all queries when BLAST did not find a match
        set<uint32_t> blastFN_noCriteria_parmik_alnLen_allQ, blastFN_noCriteria_parmik_ed_allQ;
        set<uint32_t> pmFN_blast_alnLen_allQ, pmFN_blast_ed_allQ;// BLAST alignment characteristics for all queries when PARMIK did not find a match
        for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            uint32_t blastMatchSize = 0, pmMatchSize = 0;
            string query = queries[queryInd];
            if(query.find('N') != string::npos || query.find('n') != string::npos)
            {
                cmp << "query contains N!" << endl;
                numberOfQueryContainN++; 
                alnPerQ << queryInd << " 0 0"<< endl;
                continue;
            }
            //query level parameters
            uint32_t blastTP = 0, blastFP = 0;
            uint32_t blastTP_noParmikMatches = 0, blastTP_blastOutperfomed = 0, blastTP_parmikOutperfomed = 0, blastTP_blastEqualParmik = 0;
            uint32_t parmikTP_noBlastMatches = 0;
            uint32_t blastBest_blastOutperfomed = 0, blastBest_parmikOutperfomed = 0, blastBest_blastEqualParmik = 0;
            uint32_t blastFP_editsExceed = 0, blastFP_lowAlnLen = 0, blastFP_editPos = 0;
            cmp << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            cmp << "Q : " << query << ", queryInd: " << queryInd << endl;
            BlastReader::Blast bestAlnBlast;
            auto pmRange = pmAlignments.getRange(queryInd);
            size_t pmReadPerQuery = distance(pmRange.first, pmRange.second);
            pmTotalNumberOfReadIDs += pmReadPerQuery;
            pmReadPerQuerySet.insert(pmReadPerQuery);
            auto blRrange = blastAlignments.getRange(queryInd);
            size_t blastReadPerQuery = distance(blRrange.first, blRrange.second);
            blastTotalNumberOfReadIDs += blastReadPerQuery;
            blastReadPerQuerySet.insert(blastReadPerQuery);
            cmp << "# of reads found by BLAST : " << blastReadPerQuery << endl;
            cmp << "# of reads found by PARMIK : " << pmReadPerQuery << endl;
            LevAlign bestAlnPm;
            vector<uint32_t> blastTPReadIds;
            //get the best PARMIK alignment
            if(pmReadPerQuery > 0){
                queriesPMFoundMatch++;
                for (auto it = pmRange.first; it != pmRange.second; it++) {
                    LevAlign pmAlnn = it->second;
                    if (pmAlnn.numberOfMatches + pmAlnn.numberOfInDel > bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel) // only exact matches
                    {
                        bestAlnPm = pmAlnn;
                    }
                    else if (pmAlnn.numberOfMatches + pmAlnn.numberOfInDel == bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel)
                    {
                        if (pmAlnn.numberOfInDel + pmAlnn.numberOfSub < bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub) // InDel has the same wight as substitution
                        {
                            bestAlnPm = pmAlnn;
                        }
                    }
                }
                pmMatchSize = bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel;
            }
            if(blastReadPerQuery == 0 && pmReadPerQuery == 0){
                //TN
                totalBlastTN++;
                cmp << "TN for BLAST" << endl;
            } else if(blastReadPerQuery == 0 && pmReadPerQuery > 0) {
                //FN
                blastFN++;
                blastFN_parmik_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
                blastFN_parmik_ed_allQ.insert(bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub);
                cmp << "FN for BLAST, PARMIK alnlen : " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", ed : " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", readID: " <<  bestAlnPm.readID << endl;
            } else {
                for (auto it = blRrange.first; it != blRrange.second; it++) 
                {
                    BlastReader::Blast aln = it->second;
                    string blastR = reads[aln.readId]; 
                    if(blastR.find('N') != string::npos || blastR.find('n') != string::npos)
                    {
                        numberOfBlastReadsContainingN++;
                        continue;
                    }
                    bool isFP = false;
                    if (aln.Mismatches + aln.InDels > cfg.editDistance)
                    {
                        // cmp << "with edits (Indel+Subs) [" << aln.Mismatches + aln.InDels << "] > " << cfg.editDistance << endl;
                        isFP = true;
                        blastFP_editsExceed++;
                        blastFP_editsExceed_alnLen_allQ.insert(aln.AlignmentLength);
                        blastFP_editsExceed_ed_allQ.insert(aln.Mismatches + aln.InDels);
                    } else if(aln.AlignmentLength < cfg.regionSize) {    
                        // cmp << "with low match size : " << blastMatchSize << endl;
                        isFP = true;
                        blastFP_lowAlnLen++;
                        blastFP_lowAlnLen_alnLen_allQ.insert(aln.AlignmentLength);
                        blastFP_lowAlnLen_ed_allQ.insert(aln.Mismatches + aln.InDels);
                    } else if (!checkBlastEditPositions(aln, cfg)){
                        isFP = true;
                        blastFP_editPos++;
                        blastFP_editPos_alnLen_allQ.insert(aln.AlignmentLength);
                        blastFP_editPos_ed_allQ.insert(aln.Mismatches + aln.InDels);
                        // cmp << "the read was supposed to be discarded based on our criteria (edit pos)" << endl;
                    } 
                    if (isFP){
                        blastFP++;
                    }
                    else {
                        blastTP++;
                        blastTPReadIds.push_back(aln.readId);
                        bool pmfound = false;
                        LevAlign pmaln;
                        for (auto itt = pmRange.first; itt != pmRange.second; itt++) 
                        {
                            pmaln = itt->second;
                            if((uint32_t) pmaln.readID == aln.readId)
                            {
                                pmfound = true;
                                break;
                            }
                        }
                        if (!pmfound)
                        {
                            blastTP_noParmikMatches++;
                            blastTP_noParmikMatches_alnLen_allQ.insert(aln.AlignmentLength);
                            blastTP_noParmikMatches_ed_allQ.insert(aln.Mismatches + aln.InDels);
                            cmp << "TP for BLAST, no PARMIK match, alnlen : " << aln.AlignmentLength << ", ed : " << aln.Mismatches + aln.InDels << ", readID: " <<  aln.readId << endl;
                        } else {
                            if (pmaln.numberOfMatches + pmaln.numberOfInDel > aln.AlignmentLength){ //PARMIK outperformed
                                blastTP_parmikOutperfomed++;
                                blastTP_parmikOutperfomed_alnLen_allQ.insert(aln.AlignmentLength);
                                blastTP_parmikOutperfomed_ed_allQ.insert(aln.Mismatches + aln.InDels);
                                parmikTP_parmikOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
                                parmikTP_parmikOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
                            } else if (pmaln.numberOfMatches + pmaln.numberOfInDel < aln.AlignmentLength){//BLAST outperformed
                                blastTP_blastOutperfomed++;
                                blastTP_blastOutperfomed_alnLen_allQ.insert(aln.AlignmentLength);
                                blastTP_blastOutperfomed_ed_allQ.insert(aln.Mismatches + aln.InDels);
                                parmikTP_blastOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
                                parmikTP_blastOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
                                cmp << "TP for BLAST, BLAST outperformed in terms of alnlen, Blast alnlen : " << aln.AlignmentLength << ", ed : " << aln.Mismatches + aln.InDels << ", readID: " <<  aln.readId << ", PARMIK alnlen: " << pmaln.numberOfMatches + pmaln.numberOfInDel << ", ed: " << pmaln.numberOfSub + pmaln.numberOfInDel << ", readID: " << pmaln.readID << endl;
                            } else {
                                if (pmaln.numberOfSub + pmaln.numberOfInDel < aln.Mismatches + aln.InDels)//PARMIK outperformed
                                {
                                    blastTP_parmikOutperfomed++;
                                    blastTP_parmikOutperfomed_alnLen_allQ.insert(aln.AlignmentLength);
                                    blastTP_parmikOutperfomed_ed_allQ.insert(aln.Mismatches + aln.InDels);
                                    parmikTP_parmikOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
                                    parmikTP_parmikOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
                                } else if (pmaln.numberOfSub + pmaln.numberOfInDel > aln.Mismatches + aln.InDels)//BLAST outperformed
                                {
                                    blastTP_blastOutperfomed++;
                                    blastTP_blastOutperfomed_alnLen_allQ.insert(aln.AlignmentLength);
                                    blastTP_blastOutperfomed_ed_allQ.insert(aln.Mismatches + aln.InDels);
                                    parmikTP_blastOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
                                    parmikTP_blastOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
                                    cmp << "TP for BLAST, BLAST outperformed in terms of ed, Blast alnlen : " << aln.AlignmentLength << ", ed : " << aln.Mismatches + aln.InDels << ", readID: " <<  aln.readId << ", PARMIK alnlen: " << pmaln.numberOfMatches + pmaln.numberOfInDel << ", ed: " << pmaln.numberOfSub + pmaln.numberOfInDel << ", readID: " << pmaln.readID << endl;
                                } // BLAST and PARMIK performed equally
                                else{
                                    blastTP_blastEqualParmik++;
                                    blastTP_blastEqualParmik_alnLen_allQ.insert(aln.AlignmentLength);
                                    blastTP_blastEqualParmik_ed_allQ.insert(aln.Mismatches + aln.InDels);
                                    parmikTP_blastEqualParmik_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
                                    parmikTP_blastEqualParmik_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
                                }
                            }
                        }
                        //find the best alignment for BLAST
                        if (aln.AlignmentLength  > bestAlnBlast.AlignmentLength)
                        {
                            bestAlnBlast = aln;
                        } else if (aln.AlignmentLength == bestAlnBlast.AlignmentLength)
                        {
                            if (aln.Mismatches + aln.InDels < bestAlnBlast.Mismatches + bestAlnBlast.InDels)
                            {
                                bestAlnBlast = aln;
                            }
                        }
                    }
                }
                for (auto itt = pmRange.first; itt != pmRange.second; itt++) 
                {
                    LevAlign aaln = itt->second;
                    bool blastfound = false;
                    for (const auto& element : blastTPReadIds)
                    {
                        if(element == (uint32_t)aaln.readID)
                        {
                            blastfound = true;
                            break;
                        }
                    }
                    if(!blastfound)
                    {
                        parmikTP_noBlastMatches++;
                        parmikTP_noBlastMatches_alnLen_allQ.insert(aaln.numberOfMatches + aaln.numberOfInDel);
                        parmikTP_noBlastMatches_ed_allQ.insert(aaln.numberOfSub + aaln.numberOfInDel);
                    }
                }
                parmikTP_noBlastMatches_allQ += parmikTP_noBlastMatches;
                cmp << "parmik TP that blast did not find: " << parmikTP_noBlastMatches << endl;
                cmp << "BlastTP: " << blastTP << "\t" << "BlastFP: " << blastFP << endl;
                blastMatchSize = bestAlnBlast.AlignmentLength;
                blastTPPerQuery.insert(blastTP);
                totalBlastTP += blastTP;
                blastFPPerQuery.insert(blastFP);
                totalBlastFP += blastFP;
                if(blastTP > 0){ 
                    queriesBLASTFoundMatch++;
                }
                if(pmReadPerQuery > 0 && blastReadPerQuery > 0 && blastTP == 0){//PARMIK TP
                    blastFN_noCriteriaFittedMatches++;
                    blastFN_noCriteria_parmik_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
                    blastFN_noCriteria_parmik_ed_allQ.insert(bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub);
                    cmp << "BlastFN_noCriteriaFitted and PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", PARMIK readID: " << bestAlnPm.readID << endl;
                }
                if(pmReadPerQuery == 0 && blastTP > 0){//PARMIK FN
                    pmFN_blast_alnLen_allQ.insert(bestAlnBlast.AlignmentLength);
                    pmFN_blast_ed_allQ.insert(bestAlnBlast.Mismatches + bestAlnBlast.InDels);
                    cmp << "PARMIK FN, and BLAST alnsize: " << bestAlnBlast.AlignmentLength << ", BLAST ed: " << bestAlnBlast.Mismatches + bestAlnBlast.InDels << ", BLAST readID: " << bestAlnBlast.readId << endl;
                    numnberOfPM_FN++;
                } else if(pmReadPerQuery > 0 && blastTP > 0){//compare the best alignmnets for this queyry
                    if (bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel > bestAlnBlast.AlignmentLength){ //PARMIK outperformed
                        blastBest_parmikOutperfomed++;
                        blastBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnBlast.AlignmentLength);
                        blastBest_parmikOutperfomed_ed_allQ.insert(bestAlnBlast.Mismatches + bestAlnBlast.InDels);
                        parmikBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
                        parmikBest_parmikOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
                        cmp << "Best: PARMIK outperformed PARMIK, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub <<  ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.AlignmentLength << ", BLAST ed: " << bestAlnBlast.Mismatches + bestAlnBlast.InDels <<  ", BLAST readID: " << bestAlnBlast.readId << endl;
                    } else if (bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel < bestAlnBlast.AlignmentLength){//BLAST outperformed
                        blastBest_blastOutperfomed++;
                        blastBest_blastOutperfomed_alnLen_allQ.insert(bestAlnBlast.AlignmentLength);
                        blastBest_blastOutperfomed_ed_allQ.insert(bestAlnBlast.Mismatches + bestAlnBlast.InDels);
                        parmikBest_blastOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
                        parmikBest_blastOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
                        cmp << "Best: BLAST outperformed PARMIK in terms of alnlen, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub <<  ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.AlignmentLength << ", BLAST ed: " << bestAlnBlast.Mismatches + bestAlnBlast.InDels << ", BLAST readID: " << bestAlnBlast.readId << endl;
                    } else {
                        if (bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel < bestAlnBlast.Mismatches + bestAlnBlast.InDels)//PARMIK outperformed
                        {
                            blastBest_parmikOutperfomed++;
                            blastBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnBlast.AlignmentLength);
                            blastBest_parmikOutperfomed_ed_allQ.insert(bestAlnBlast.Mismatches + bestAlnBlast.InDels);
                            parmikBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
                            parmikBest_parmikOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
                            cmp << "Best: PARMIK outperformed PARMIK, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.AlignmentLength << ", BLAST ed: " << bestAlnBlast.Mismatches + bestAlnBlast.InDels <<  ", BLAST readID: " << bestAlnBlast.readId << endl;
                        } else if (bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel > bestAlnBlast.Mismatches + bestAlnBlast.InDels)//BLAST outperformed
                        {
                            blastBest_blastOutperfomed++;
                            blastBest_blastOutperfomed_alnLen_allQ.insert(bestAlnBlast.AlignmentLength);
                            blastBest_blastOutperfomed_ed_allQ.insert(bestAlnBlast.Mismatches + bestAlnBlast.InDels);
                            parmikBest_blastOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
                            parmikBest_blastOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
                            cmp << "Best: BLAST outperformed PARMIK in terms of ed, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub <<  ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.AlignmentLength << ", BLAST ed: " << bestAlnBlast.Mismatches + bestAlnBlast.InDels <<  ", BLAST readID: " << bestAlnBlast.readId << endl;
                        } // BLAST and PARMIK performed equally
                        else{
                            blastBest_blastEqualParmik++;
                            blastBest_blastEqualParmik_alnLen_allQ.insert(bestAlnBlast.AlignmentLength);
                            blastBest_blastEqualParmik_ed_allQ.insert(bestAlnBlast.Mismatches + bestAlnBlast.InDels);
                            parmikBest_blastEqualParmik_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
                            parmikBest_blastEqualParmik_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
                            cmp << "Best: BLAST & PARMIK perfromed equally, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.AlignmentLength << ", BLAST ed: " << bestAlnBlast.Mismatches + bestAlnBlast.InDels << ", BLAST readID: " << bestAlnBlast.readId << endl;
                        }
                    }
                }
            }
            blastTP_noParmikMatches_allQ += blastTP_noParmikMatches; 
            blastTP_blastOutperfomed_allQ += blastTP_blastOutperfomed; 
            blastTP_parmikOutperfomed_allQ += blastTP_parmikOutperfomed;
            blastTP_blastEqualParmik_allQ += blastTP_blastEqualParmik;
            cmp << "blastTP_noParmikMatches: " << blastTP_noParmikMatches << ", blastTP_blastOutperfomed: " << blastTP_blastOutperfomed << ", blastTP_parmikOutperfomed: " << blastTP_parmikOutperfomed << ", blastTP_blastEqualParmik: " << blastTP_blastEqualParmik << endl;
            blastFP_editsExceed_allQ += blastFP_editsExceed;
            blastFP_lowAlnLen_allQ += blastFP_lowAlnLen; 
            blastFP_editPos_allQ += blastFP_editPos;
            cmp << "blastFP_editsExceed: " << blastFP_editsExceed << ", blastFP_lowAlnLen: " << blastFP_lowAlnLen << ", blastFP_editPos: " << blastFP_editPos << endl;
            blastBest_parmikOutperfomed_allQ += blastBest_parmikOutperfomed;
            blastBest_blastOutperfomed_allQ += blastBest_blastOutperfomed;
            blastBest_blastEqualParmik_allQ += blastBest_blastEqualParmik;
            cmp << "blastBest_parmikOutperfomed: " << blastBest_parmikOutperfomed << ", blastBest_blastOutperfomed: " << blastBest_blastOutperfomed << ", blastBest_blastEqualParmik: " << blastBest_blastEqualParmik << endl;
            alnPmBLASTHisto.push_back(make_pair(pmMatchSize, blastMatchSize)); 
            //PM and BLAST alignments per query
            alnPerQ << queryInd << " " << pmAlignments.container_.count(queryInd) << " " << blastTP << endl;
        }
        Utilities<uint32_t> util;   
        cmp << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Overall Comparison Results>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cmp << left << setw(80) << "# Of queries : " << queryCount << endl;
        cmp << left << setw(80) << "# Of queries that PARMIK found match : " << queriesPMFoundMatch << endl;
        cmp << left << setw(80) << "# Of queries that BLAST found match (TP) : " << queriesBLASTFoundMatch << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries that BLAST didn't find TN : " << totalBlastTN << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries that BLAST didn't find FN (Total) : " << blastFN + blastFN_noCriteriaFittedMatches << endl;
        cmp << left << setw(80) << "# Of queries that BLAST didn't find any match FN : " << blastFN << endl;
        pair<uint32_t, uint32_t> avgBlastFN_parmik_alnLen_allQ = util.calculateStatistics2(blastFN_parmik_alnLen_allQ);
        pair<uint32_t, uint32_t> avgBlastFN_parmik_ed_allQ = util.calculateStatistics2(blastFN_parmik_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BLAST FN: " << avgBlastFN_parmik_alnLen_allQ.first << ", " << avgBlastFN_parmik_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BLAST FN: " << avgBlastFN_parmik_ed_allQ.first << ", " << avgBlastFN_parmik_ed_allQ.second << endl;
        cmp << left << setw(80) << "# Of queries that BLAST didn't find based on criteria FN : " << blastFN_noCriteriaFittedMatches << endl;
        pair<uint32_t, uint32_t> avgblastFN_noCriteria_parmik_alnLen_allQ = util.calculateStatistics2(blastFN_noCriteria_parmik_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastFN_noCriteria_parmik_ed_allQ = util.calculateStatistics2(blastFN_noCriteria_parmik_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BLAST FN based on criteria: " << avgblastFN_noCriteria_parmik_alnLen_allQ.first << ", " << avgblastFN_noCriteria_parmik_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BLAST FN based on criteria: " << avgblastFN_noCriteria_parmik_ed_allQ.first << ", " << avgblastFN_noCriteria_parmik_ed_allQ.second << endl;
       
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries that PARMIK didn't find FN : " << numnberOfPM_FN << endl;
        pair<uint32_t, uint32_t> avgPmFN_blast_alnLen_allQ = util.calculateStatistics2(pmFN_blast_alnLen_allQ);
        pair<uint32_t, uint32_t> avgPmFN_blast_ed_allQ = util.calculateStatistics2(pmFN_blast_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where PARMIK FN: " << avgPmFN_blast_alnLen_allQ.first << ", " << avgPmFN_blast_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where PARMIK FN: " << avgPmFN_blast_ed_allQ.first << ", " << avgPmFN_blast_ed_allQ.second << endl;
       
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of readIDs found by PARMIK (total): " << pmTotalNumberOfReadIDs << endl;
        cmp << left << setw(80) << "# Of readIDs found by BLAST (total): " << blastTotalNumberOfReadIDs << endl;
        cmp << left << setw(80) << "# Of readIDs found by BLAST containing N (total): " << numberOfBlastReadsContainingN << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of total BLAST (TP): " << totalBlastTP << endl;
        pair<uint32_t, uint32_t> avgblastTPPerQuery = util.calculateStatistics2(blastTPPerQuery);
        cmp << left << setw(80) << "(avg, median) BLAST's TP per query: " << avgblastTPPerQuery.first << ", " << avgblastTPPerQuery.second << endl;

        cmp << left << setw(80) << "# Of BLAST (TP) where PARMIK didn't find match: " << blastTP_noParmikMatches_allQ << endl;
        pair<uint32_t, uint32_t> avgBlastTP_noParmikMatches_alnLen_allQ = util.calculateStatistics2(blastTP_noParmikMatches_alnLen_allQ);
        pair<uint32_t, uint32_t> avgBlastTP_noParmikMatches_ed_allQ = util.calculateStatistics2(blastTP_noParmikMatches_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where PARMIK didn't find match: " << avgBlastTP_noParmikMatches_alnLen_allQ.first << ", " << avgBlastTP_noParmikMatches_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where PARMIK didn't find match: " << avgBlastTP_noParmikMatches_ed_allQ.first << ", " << avgBlastTP_noParmikMatches_ed_allQ.second << endl;
        
        cmp << left << setw(80) << "# Of BLAST (TP) where PARMIK outperformed: " << blastTP_parmikOutperfomed_allQ << endl;
        pair<uint32_t, uint32_t> avgblastTP_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(blastTP_parmikOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastTP_parmikOutperfomed_ed_allQ = util.calculateStatistics2(blastTP_parmikOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where PARMIK outperformed: " << avgblastTP_parmikOutperfomed_alnLen_allQ.first << ", " << avgblastTP_parmikOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where PARMIK outperformed: " << avgblastTP_parmikOutperfomed_ed_allQ.first << ", " << avgblastTP_parmikOutperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgparmikTP_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikTP_parmikOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgparmikTP_parmikOutperfomed_ed_allQ = util.calculateStatistics2(parmikTP_parmikOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where PARMIK outperformed: " << avgparmikTP_parmikOutperfomed_alnLen_allQ.first << ", " << avgparmikTP_parmikOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where PARMIK outperformed: " << avgparmikTP_parmikOutperfomed_ed_allQ.first << ", " << avgparmikTP_parmikOutperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of BLAST (TP) where BLAST outperformed: " << blastTP_blastOutperfomed_allQ << endl;
        pair<uint32_t, uint32_t> avgblastTP_blastOutperfomed_alnLen_allQ = util.calculateStatistics2(blastTP_blastOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastTP_blastOutperfomed_ed_allQ = util.calculateStatistics2(blastTP_blastOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where BLAST outperformed: " << avgblastTP_blastOutperfomed_alnLen_allQ.first << ", " << avgblastTP_blastOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where BLAST outperformed: " << avgblastTP_blastOutperfomed_ed_allQ.first << ", " << avgblastTP_blastOutperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgparmikTP_blastOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikTP_blastOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgparmikTP_blastOutperfomed_ed_allQ = util.calculateStatistics2(parmikTP_blastOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BLAST outperformed: " << avgparmikTP_blastOutperfomed_alnLen_allQ.first << ", " << avgparmikTP_blastOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BLAST outperformed: " << avgparmikTP_blastOutperfomed_ed_allQ.first << ", " << avgparmikTP_blastOutperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of BLAST (TP) where they perfromed equaly: " << blastTP_blastEqualParmik_allQ << endl;
        pair<uint32_t, uint32_t> avgblastTP_blastEqualParmik_alnLen_allQ = util.calculateStatistics2(blastTP_blastEqualParmik_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastTP_blastEqualParmik_ed_allQ = util.calculateStatistics2(blastTP_blastEqualParmik_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where they perfromed equaly: " << avgblastTP_blastEqualParmik_alnLen_allQ.first << ", " << avgblastTP_blastEqualParmik_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where they perfromed equaly: " << avgblastTP_blastEqualParmik_ed_allQ.first << ", " << avgblastTP_blastEqualParmik_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgparmikTP_blastEqualParmik_alnLen_allQ = util.calculateStatistics2(parmikTP_blastEqualParmik_alnLen_allQ);
        pair<uint32_t, uint32_t> avgparmikTP_blastEqualParmik_ed_allQ = util.calculateStatistics2(parmikTP_blastEqualParmik_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where they perfromed equaly: " << avgparmikTP_blastEqualParmik_alnLen_allQ.first << ", " << avgparmikTP_blastEqualParmik_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where they perfromed equaly: " << avgparmikTP_blastEqualParmik_ed_allQ.first << ", " << avgparmikTP_blastEqualParmik_ed_allQ.second << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of BLAST best where PARMIK outperformed: " << blastBest_parmikOutperfomed_allQ << endl;
        pair<uint32_t, uint32_t> avgblastBest_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(blastBest_parmikOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastBest_parmikOutperfomed_ed_allQ = util.calculateStatistics2(blastBest_parmikOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's best alignment length where PARMIK outperformed: " << avgblastBest_parmikOutperfomed_alnLen_allQ.first << ", " << avgblastBest_parmikOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's best edit distance where PARMIK outperformed: " << avgblastBest_parmikOutperfomed_ed_allQ.first << ", " << avgblastBest_parmikOutperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgparmikBest_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikBest_parmikOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgparmikBest_parmikOutperfomed_ed_allQ = util.calculateStatistics2(parmikBest_parmikOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's best alignment length where PARMIK outperformed: " << avgparmikBest_parmikOutperfomed_alnLen_allQ.first << ", " << avgparmikBest_parmikOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's best edit distance where PARMIK outperformed: " << avgparmikBest_parmikOutperfomed_ed_allQ.first << ", " << avgparmikBest_parmikOutperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of BLAST best where BLAST outperformed: " << blastBest_blastOutperfomed_allQ  << endl;
        pair<uint32_t, uint32_t> avgblastBest_blastOutperfomed_alnLen_allQ = util.calculateStatistics2(blastBest_blastOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastBest_blastOutperfomed_ed_allQ = util.calculateStatistics2(blastBest_blastOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's best alignment length where BLAST outperformed: " << avgblastBest_blastOutperfomed_alnLen_allQ.first << ", " << avgblastBest_blastOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's best edit distance where BLAST outperformed: " << avgblastBest_blastOutperfomed_ed_allQ.first << ", " << avgblastBest_blastOutperfomed_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgparmikBest_blastOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikBest_blastOutperfomed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgparmikBest_blastOutperfomed_ed_allQ = util.calculateStatistics2(parmikBest_blastOutperfomed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's best alignment length where BLAST outperformed: " << avgparmikBest_blastOutperfomed_alnLen_allQ.first << ", " << avgparmikBest_blastOutperfomed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's best edit distance where BLAST outperformed: " << avgparmikBest_blastOutperfomed_ed_allQ.first << ", " << avgparmikBest_blastOutperfomed_ed_allQ.second << endl;

        cmp << left << setw(80) << "# Of BLAST best where they perfromed equaly: " << blastBest_blastEqualParmik_allQ << endl;
        pair<uint32_t, uint32_t> avgblastBest_blastEqualParmik_alnLen_allQ = util.calculateStatistics2(blastBest_blastEqualParmik_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastBest_blastEqualParmik_ed_allQ = util.calculateStatistics2(blastBest_blastEqualParmik_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's best alignment length where they perfromed equaly: " << avgblastBest_blastEqualParmik_alnLen_allQ.first << ", " << avgblastBest_blastEqualParmik_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's best edit distance where they perfromed equaly: " << avgblastBest_blastEqualParmik_ed_allQ.first << ", " << avgblastBest_blastEqualParmik_ed_allQ.second << endl;
        pair<uint32_t, uint32_t> avgparmikBest_blastEqualParmik_alnLen_allQ = util.calculateStatistics2(parmikBest_blastEqualParmik_alnLen_allQ);
        pair<uint32_t, uint32_t> avgparmikBest_blastEqualParmik_ed_allQ = util.calculateStatistics2(parmikBest_blastEqualParmik_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's best alignment length where they perfromed equaly: " << avgparmikBest_blastEqualParmik_alnLen_allQ.first << ", " << avgparmikBest_blastEqualParmik_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's best edit distance where they perfromed equaly: " << avgparmikBest_blastEqualParmik_ed_allQ.first << ", " << avgparmikBest_blastEqualParmik_ed_allQ.second << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of total BLAST (FP): " << totalBlastFP << endl;
        pair<uint32_t, uint32_t> avgblastFPPerQuery = util.calculateStatistics2(blastFPPerQuery);
        cmp << left << setw(80) << "(avg, median) BLAST's FP per query: " << avgblastFPPerQuery.first << ", " << avgblastFPPerQuery.second << endl;

        cmp << left << setw(80) << "# Of BLAST (FP) that exceed max edit distance: " << blastFP_editsExceed_allQ << endl;
        pair<uint32_t, uint32_t> avgblastFP_editsExceed_alnLen_allQ = util.calculateStatistics2(blastFP_editsExceed_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastFP_editsExceed_ed_allQ = util.calculateStatistics2(blastFP_editsExceed_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where exceed max edit distance: " << avgblastFP_editsExceed_alnLen_allQ.first << ", " << avgblastFP_editsExceed_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where exceed max edit distance: " << avgblastFP_editsExceed_ed_allQ.first << ", " << avgblastFP_editsExceed_ed_allQ.second << endl;
        
        cmp << left << setw(80) << "# Of BLAST (FP) that alignment length < R: " << blastFP_lowAlnLen_allQ << endl;
        pair<uint32_t, uint32_t> avgblastFP_lowAlnLen_alnLen_allQ = util.calculateStatistics2(blastFP_lowAlnLen_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastFP_lowAlnLen_ed_allQ = util.calculateStatistics2(blastFP_lowAlnLen_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where alignment length < R: " << avgblastFP_lowAlnLen_alnLen_allQ.first << ", " << avgblastFP_lowAlnLen_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where alignment length < R: " << avgblastFP_lowAlnLen_ed_allQ.first << ", " << avgblastFP_lowAlnLen_ed_allQ.second << endl;
        
        cmp << left << setw(80) << "# Of BLAST (FP) that edit pos breaks criteria: " << blastFP_editPos_allQ << endl;
        pair<uint32_t, uint32_t> avgblastFP_editPos_alnLen_allQ = util.calculateStatistics2(blastFP_editPos_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastFP_editPos_ed_allQ = util.calculateStatistics2(blastFP_editPos_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where edit pos breaks criteria: " << avgblastFP_editPos_alnLen_allQ.first << ", " << avgblastFP_editPos_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where edit pos breaks criteria: " << avgblastFP_editPos_ed_allQ.first << ", " << avgblastFP_editPos_ed_allQ.second << endl;

        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of PARMIK (TP) where BLAST didn't find match: " << parmikTP_noBlastMatches_allQ << endl;
        pair<uint32_t, uint32_t> avgparmikTP_noBlastMatches_alnLen_allQ = util.calculateStatistics2(parmikTP_noBlastMatches_alnLen_allQ);
        pair<uint32_t, uint32_t> avgparmikTP_noBlastMatches_ed_allQ = util.calculateStatistics2(parmikTP_noBlastMatches_ed_allQ);
        cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BLAST didn't find match: " << avgparmikTP_noBlastMatches_alnLen_allQ.first << ", " << avgparmikTP_noBlastMatches_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BLAST didn't find match: " << avgparmikTP_noBlastMatches_ed_allQ.first << ", " << avgparmikTP_noBlastMatches_ed_allQ.second << endl;
        
        cmp << "------------------------------------------------------------------------------------------" << endl;
        cmp << left << setw(80) << "# Of queries contains N: " << numberOfQueryContainN << endl;

        cmp << "------------------------------------------------------------------------------------------" << endl;
        pair<uint32_t, uint32_t> pmReadPerQuery_allQ  = util.calculateStatistics2(pmReadPerQuerySet);
        pair<uint32_t, uint32_t> blastReadPerQuery_allQ  = util.calculateStatistics2(blastReadPerQuerySet);
        cmp << left << setw(80) << "(avg, median) PARMIKS's Matches per Query: " << pmReadPerQuery_allQ.first << ", " << pmReadPerQuery_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's Matches per Query: " << blastReadPerQuery_allQ.first << ", " << blastReadPerQuery_allQ.second << endl;

        cmp.close();
        alnPerQ.close();
    }
};

#endif