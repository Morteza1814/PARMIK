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

    void comparePmWithBlast(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress, IndexContainer<uint32_t, Alignment>& pmAlignments, vector<pair<uint32_t, uint32_t>>& alnPmBLASTHisto, const uint32_t queryCount, const string alnPerQueryFileAddress, const string parmikFnReadsFileAddress, const string bestAlnCmpFileAddress)
    {
        ofstream cmp(comparisonResultsFileAddress);
        cout << "cmp file address: " << comparisonResultsFileAddress << endl;
        if(!cmp.is_open())
            cout << "cmp file not opened" << endl;
        ofstream alnPerQ(alnPerQueryFileAddress);
        ofstream pmFn(parmikFnReadsFileAddress);
        ofstream bestAlnCmp(bestAlnCmpFileAddress);
        BlastReader blastReader(cfg.otherToolOutputFileAddress);
        IndexContainer<uint32_t, Alignment> blastAlignments;
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
            Alignment bestAlnBlast;
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
            Alignment bestAlnPm;
            vector<uint32_t> blastTPReadIds;
            //get the best PARMIK alignment
            if(pmReadPerQuery > 0){
                queriesPMFoundMatch++;
                for (auto it = pmRange.first; it != pmRange.second; it++) {
                    Alignment pmAlnn = it->second;
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
                pmMatchSize = bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions;
            }
            if(blastReadPerQuery == 0 || pmReadPerQuery == 0){
                bestAlnCmp << "- -" << endl;
            }
            if(blastReadPerQuery == 0 && pmReadPerQuery == 0){
                //TN
                totalBlastTN++;
                cmp << "TN for BLAST" << endl;
            } else if(blastReadPerQuery == 0 && pmReadPerQuery > 0) {
                //FN
                blastFN++;
                blastFN_parmik_alnLen_allQ.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions);
                blastFN_parmik_ed_allQ.insert(bestAlnPm.inDels + bestAlnPm.substitutions);
                cmp << "FN for BLAST, PARMIK alnlen : " << bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << ", ed : " << bestAlnPm.inDels + bestAlnPm.substitutions << ", readID: " <<  bestAlnPm.readID << endl;
            } else {
                for (auto it = blRrange.first; it != blRrange.second; it++) 
                {
                    Alignment aln = it->second;
                    string blastR = reads[aln.readID]; 
                    if(blastR.find('N') != string::npos || blastR.find('n') != string::npos)
                    {
                        numberOfBlastReadsContainingN++;
                        continue;
                    }
                    bool isFP = false;
                    PostFilter pf(cfg.regionSize, cfg.editDistance, cfg.contigSize, cfg.identityPercentage);
                    string cigarStr = getCigarStr(aln.queryRegionStartPos - 1, aln.alignedQuery, aln.alignedRead, cfg.contigSize);
                    //TODO: check for the new implementation
                    bool criteriaCheck = pf.checkAndUpdateBasedOnAlingmentCriteria(aln);
                    if(!criteriaCheck){
                        isFP = true;
                        if(aln.criteriaCode == 0x10) {
                            blastFP_lowAlnLen++;
                            blastFP_lowAlnLen_alnLen_allQ.insert(aln.partialMatchSize);
                            blastFP_lowAlnLen_ed_allQ.insert(aln.substitutions + aln.inDels);
                        }
                    }
                    // if (aln.criteriaCode == 0x08)
                    // {
                    //     // cmp << "with edits (Indel+Subs) [" << aln.substitutions + aln.inDels << "] > " << cfg.editDistance << endl;
                    //     isFP = true;
                    //     blastFP_editsExceed++;
                    //     blastFP_editsExceed_alnLen_allQ.insert(aln.partialMatchSize);
                    //     blastFP_editsExceed_ed_allQ.insert(aln.substitutions + aln.inDels);
                    // } else if(aln.criteriaCode == 0x10) {    
                    //     // cmp << "with low match size : " << blastMatchSize << endl;
                    //     isFP = true;
                    //     blastFP_lowAlnLen++;
                    //     blastFP_lowAlnLen_alnLen_allQ.insert(aln.partialMatchSize);
                    //     blastFP_lowAlnLen_ed_allQ.insert(aln.substitutions + aln.inDels);
                    // } else if (aln.criteriaCode >= 0x04) {
                    //     isFP = true;
                    //     blastFP_editPos++;
                    //     blastFP_editPos_alnLen_allQ.insert(aln.partialMatchSize);
                    //     blastFP_editPos_ed_allQ.insert(aln.substitutions + aln.inDels);
                    //     // cmp << "the read was supposed to be discarded based on our criteria (edit pos)" << endl;
                    // } 
        
                    if (isFP){
                        blastFP++;
                    }
                    else {
                        blastTP++;
                        blastTPReadIds.push_back(aln.readID);
                        bool pmfound = false;
                        Alignment pmaln;
                        for (auto itt = pmRange.first; itt != pmRange.second; itt++) 
                        {
                            Alignment tmpAln = itt->second;
                            if(tmpAln.readID == aln.readID)
                            {
                                if (!pmfound)
                                    pmaln = itt->second;
                                pmfound = true;
                                if ((tmpAln.matches + tmpAln.inDels + tmpAln.substitutions > pmaln.matches + pmaln.inDels + pmaln.substitutions) || 
                                ((tmpAln.matches + tmpAln.inDels + tmpAln.substitutions == pmaln.matches + pmaln.inDels + pmaln.substitutions) && 
                                (tmpAln.inDels + tmpAln.substitutions <  pmaln.inDels + pmaln.substitutions))) 
                                {
                                    pmaln = tmpAln;
                                }
                            }
                        }
                        if (!pmfound)
                        {
                            blastTP_noParmikMatches++;
                            blastTP_noParmikMatches_alnLen_allQ.insert(aln.partialMatchSize);
                            blastTP_noParmikMatches_ed_allQ.insert(aln.substitutions + aln.inDels);
                            cmp << "TP for BLAST, no PARMIK match, alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << endl;
                            pmFn << queryInd << "\t" << aln.readID << "\t" << aln.flag << endl;
                        } else {
                            if (pmaln.matches + pmaln.inDels + pmaln.substitutions > aln.partialMatchSize){ //PARMIK outperformed
                                blastTP_parmikOutperfomed++;
                                blastTP_parmikOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                blastTP_parmikOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                parmikTP_parmikOutperfomed_alnLen_allQ.insert(pmaln.matches + pmaln.inDels + pmaln.substitutions);
                                parmikTP_parmikOutperfomed_ed_allQ.insert(pmaln.substitutions + pmaln.inDels);
                            } else if (pmaln.matches + pmaln.inDels + pmaln.substitutions < aln.partialMatchSize){//BLAST outperformed
                                blastTP_blastOutperfomed++;
                                blastTP_blastOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                blastTP_blastOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                parmikTP_blastOutperfomed_alnLen_allQ.insert(pmaln.matches + pmaln.inDels + pmaln.substitutions);
                                parmikTP_blastOutperfomed_ed_allQ.insert(pmaln.substitutions + pmaln.inDels);
                                cmp << "TP for BLAST, BLAST outperformed in terms of alnlen, Blast alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", PARMIK alnlen: " << pmaln.matches + pmaln.inDels + pmaln.substitutions << ", ed: " << pmaln.substitutions + pmaln.inDels << ", readID: " << pmaln.readID << endl;
                            } else {
                                if (pmaln.substitutions + pmaln.inDels < aln.substitutions + aln.inDels)//PARMIK outperformed
                                {
                                    blastTP_parmikOutperfomed++;
                                    blastTP_parmikOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                    blastTP_parmikOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                    parmikTP_parmikOutperfomed_alnLen_allQ.insert(pmaln.matches + pmaln.inDels + pmaln.substitutions);
                                    parmikTP_parmikOutperfomed_ed_allQ.insert(pmaln.substitutions + pmaln.inDels);
                                    cmp << "TP for BLAST, PARMIK outperformed in terms of ed, Blast alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", PARMIK alnlen: " << pmaln.matches + pmaln.inDels << ", ed: " << pmaln.substitutions + pmaln.inDels << ", readID: " << pmaln.readID << endl;
                                } else if (pmaln.substitutions + pmaln.inDels > aln.substitutions + aln.inDels)//BLAST outperformed
                                {
                                    blastTP_blastOutperfomed++;
                                    blastTP_blastOutperfomed_alnLen_allQ.insert(aln.partialMatchSize);
                                    blastTP_blastOutperfomed_ed_allQ.insert(aln.substitutions + aln.inDels);
                                    parmikTP_blastOutperfomed_alnLen_allQ.insert(pmaln.matches + pmaln.inDels + pmaln.substitutions);
                                    parmikTP_blastOutperfomed_ed_allQ.insert(pmaln.substitutions + pmaln.inDels);
                                    cmp << "TP for BLAST, BLAST outperformed in terms of ed, Blast alnlen : " << aln.partialMatchSize << ", ed : " << aln.substitutions + aln.inDels << ", readID: " <<  aln.readID << ", PARMIK alnlen: " << pmaln.matches + pmaln.inDels << ", ed: " << pmaln.substitutions + pmaln.inDels << ", readID: " << pmaln.readID << endl;
                                } // BLAST and PARMIK performed equally
                                else{
                                    blastTP_blastEqualParmik++;
                                    blastTP_blastEqualParmik_alnLen_allQ.insert(aln.partialMatchSize);
                                    blastTP_blastEqualParmik_ed_allQ.insert(aln.substitutions + aln.inDels);
                                    parmikTP_blastEqualParmik_alnLen_allQ.insert(pmaln.matches + pmaln.inDels + pmaln.substitutions);
                                    parmikTP_blastEqualParmik_ed_allQ.insert(pmaln.substitutions + pmaln.inDels);
                                }
                            }
                        }
                        //find the best alignment for BLAST
                        if (aln.partialMatchSize  > bestAlnBlast.partialMatchSize)
                        {
                            bestAlnBlast = aln;
                        } else if (aln.partialMatchSize == bestAlnBlast.partialMatchSize)
                        {
                            if (aln.substitutions + aln.inDels < bestAlnBlast.substitutions + bestAlnBlast.inDels)
                            {
                                bestAlnBlast = aln;
                            }
                        }
                    }
                }
                for (auto itt = pmRange.first; itt != pmRange.second; itt++) 
                {
                    Alignment aaln = itt->second;
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
                        parmikTP_noBlastMatches_alnLen_allQ.insert(aaln.matches + aaln.inDels + aaln.substitutions);
                        parmikTP_noBlastMatches_ed_allQ.insert(aaln.substitutions + aaln.inDels);
                    }
                }
                parmikTP_noBlastMatches_allQ += parmikTP_noBlastMatches;
                cmp << "parmik TP that blast did not find: " << parmikTP_noBlastMatches << endl;
                cmp << "BlastTP: " << blastTP << "\t" << "BlastFP: " << blastFP << endl;
                blastMatchSize = bestAlnBlast.partialMatchSize;
                blastTPPerQuery.insert(blastTP);
                totalBlastTP += blastTP;
                blastFPPerQuery.insert(blastFP);
                totalBlastFP += blastFP;
                if(blastTP > 0){ 
                    queriesBLASTFoundMatch++;
                }
                if(pmReadPerQuery > 0 && blastReadPerQuery > 0 && blastTP == 0){//PARMIK TP
                    blastFN_noCriteriaFittedMatches++;
                    blastFN_noCriteria_parmik_alnLen_allQ.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions);
                    blastFN_noCriteria_parmik_ed_allQ.insert(bestAlnPm.inDels + bestAlnPm.substitutions);
                    cmp << "BlastFN_noCriteriaFitted and PARMIK alnsize: " << bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK ed: " << bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK readID: " << bestAlnPm.readID << endl;
                }
                if(pmReadPerQuery == 0 && blastTP > 0){//PARMIK FN
                    pmFN_blast_alnLen_allQ.insert(bestAlnBlast.partialMatchSize);
                    pmFN_blast_ed_allQ.insert(bestAlnBlast.substitutions + bestAlnBlast.inDels);
                    cmp << "PARMIK FN, and BLAST alnsize: " << bestAlnBlast.partialMatchSize << ", BLAST ed: " << bestAlnBlast.substitutions + bestAlnBlast.inDels << ", BLAST readID: " << bestAlnBlast.readID << endl;
                    numnberOfPM_FN++;
                } else if(pmReadPerQuery > 0 && blastTP > 0){//compare the best alignmnets for this queyry
                    if (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions > bestAlnBlast.partialMatchSize){ //PARMIK outperformed
                        blastBest_parmikOutperfomed++;
                        blastBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnBlast.partialMatchSize);
                        blastBest_parmikOutperfomed_ed_allQ.insert(bestAlnBlast.substitutions + bestAlnBlast.inDels);
                        parmikBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions);
                        parmikBest_parmikOutperfomed_ed_allQ.insert(bestAlnPm.substitutions + bestAlnPm.inDels);
                        cmp << "Best: PARMIK outperformed BLAST, PARMIK alnsize: " << bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK ed: " << bestAlnPm.inDels + bestAlnPm.substitutions <<  ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.partialMatchSize << ", BLAST ed: " << bestAlnBlast.substitutions + bestAlnBlast.inDels <<  ", BLAST readID: " << bestAlnBlast.readID << endl;
                        bestAlnCmp << (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) - (bestAlnBlast.partialMatchSize) << " 0" << endl;
                    } else if (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions < bestAlnBlast.partialMatchSize){//BLAST outperformed
                        blastBest_blastOutperfomed++;
                        blastBest_blastOutperfomed_alnLen_allQ.insert(bestAlnBlast.partialMatchSize);
                        blastBest_blastOutperfomed_ed_allQ.insert(bestAlnBlast.substitutions + bestAlnBlast.inDels);
                        parmikBest_blastOutperfomed_alnLen_allQ.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions);
                        parmikBest_blastOutperfomed_ed_allQ.insert(bestAlnPm.substitutions + bestAlnPm.inDels);
                        cmp << "Best: BLAST outperformed PARMIK in terms of alnlen, PARMIK alnsize: " << bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK ed: " << bestAlnPm.inDels + bestAlnPm.substitutions <<  ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.partialMatchSize << ", BLAST ed: " << bestAlnBlast.substitutions + bestAlnBlast.inDels << ", BLAST readID: " << bestAlnBlast.readID << endl;
                        bestAlnCmp << "0 " << (bestAlnBlast.partialMatchSize) - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) << endl;
                    } else {
                        if (bestAlnPm.substitutions + bestAlnPm.inDels < bestAlnBlast.substitutions + bestAlnBlast.inDels)//PARMIK outperformed
                        {
                            blastBest_parmikOutperfomed++;
                            blastBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnBlast.partialMatchSize);
                            blastBest_parmikOutperfomed_ed_allQ.insert(bestAlnBlast.substitutions + bestAlnBlast.inDels);
                            parmikBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions);
                            parmikBest_parmikOutperfomed_ed_allQ.insert(bestAlnPm.substitutions + bestAlnPm.inDels);
                            cmp << "Best: PARMIK outperformed BLAST, PARMIK alnsize: " << bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK ed: " << bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.partialMatchSize << ", BLAST ed: " << bestAlnBlast.substitutions + bestAlnBlast.inDels <<  ", BLAST readID: " << bestAlnBlast.readID << endl;
                            bestAlnCmp << (bestAlnBlast.substitutions + bestAlnBlast.inDels) - (bestAlnPm.substitutions + bestAlnPm.inDels) << " 0" << endl;
                        } else if (bestAlnPm.substitutions + bestAlnPm.inDels > bestAlnBlast.substitutions + bestAlnBlast.inDels)//BLAST outperformed
                        {
                            blastBest_blastOutperfomed++;
                            blastBest_blastOutperfomed_alnLen_allQ.insert(bestAlnBlast.partialMatchSize);
                            blastBest_blastOutperfomed_ed_allQ.insert(bestAlnBlast.substitutions + bestAlnBlast.inDels);
                            parmikBest_blastOutperfomed_alnLen_allQ.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions);
                            parmikBest_blastOutperfomed_ed_allQ.insert(bestAlnPm.substitutions + bestAlnPm.inDels);
                            cmp << "Best: BLAST outperformed PARMIK in terms of ed, PARMIK alnsize: " << bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK ed: " << bestAlnPm.inDels + bestAlnPm.substitutions <<  ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.partialMatchSize << ", BLAST ed: " << bestAlnBlast.substitutions + bestAlnBlast.inDels <<  ", BLAST readID: " << bestAlnBlast.readID << endl;
                            bestAlnCmp << "0 " << (bestAlnPm.substitutions + bestAlnPm.inDels) - (bestAlnBlast.substitutions + bestAlnBlast.inDels) << endl;
                        } // BLAST and PARMIK performed equally
                        else{
                            blastBest_blastEqualParmik++;
                            blastBest_blastEqualParmik_alnLen_allQ.insert(bestAlnBlast.partialMatchSize);
                            blastBest_blastEqualParmik_ed_allQ.insert(bestAlnBlast.substitutions + bestAlnBlast.inDels);
                            parmikBest_blastEqualParmik_alnLen_allQ.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions);
                            parmikBest_blastEqualParmik_ed_allQ.insert(bestAlnPm.substitutions + bestAlnPm.inDels);
                            cmp << "Best: BLAST & PARMIK perfromed equally, PARMIK alnsize: " << bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK ed: " << bestAlnPm.inDels + bestAlnPm.substitutions << ", PARMIK readID: " << bestAlnPm.readID << ", BLAST alnsize: " << bestAlnBlast.partialMatchSize << ", BLAST ed: " << bestAlnBlast.substitutions + bestAlnBlast.inDels << ", BLAST readID: " << bestAlnBlast.readID << endl;
                            bestAlnCmp << "0 0" << endl;
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

        // cmp << left << setw(80) << "# Of BLAST (FP) that exceed max edit distance: " << blastFP_editsExceed_allQ << endl;
        // pair<uint32_t, uint32_t> avgblastFP_editsExceed_alnLen_allQ = util.calculateStatistics2(blastFP_editsExceed_alnLen_allQ);
        // pair<uint32_t, uint32_t> avgblastFP_editsExceed_ed_allQ = util.calculateStatistics2(blastFP_editsExceed_ed_allQ);
        // cmp << left << setw(80) << "(avg, median) BLAST's alignment length where exceed max edit distance: " << avgblastFP_editsExceed_alnLen_allQ.first << ", " << avgblastFP_editsExceed_alnLen_allQ.second << endl;
        // cmp << left << setw(80) << "(avg, median) BLAST's edit distance where exceed max edit distance: " << avgblastFP_editsExceed_ed_allQ.first << ", " << avgblastFP_editsExceed_ed_allQ.second << endl;
        
        cmp << left << setw(80) << "# Of BLAST (FP) that alignment length < R: " << blastFP_lowAlnLen_allQ << endl;
        pair<uint32_t, uint32_t> avgblastFP_lowAlnLen_alnLen_allQ = util.calculateStatistics2(blastFP_lowAlnLen_alnLen_allQ);
        pair<uint32_t, uint32_t> avgblastFP_lowAlnLen_ed_allQ = util.calculateStatistics2(blastFP_lowAlnLen_ed_allQ);
        cmp << left << setw(80) << "(avg, median) BLAST's alignment length where alignment length < R: " << avgblastFP_lowAlnLen_alnLen_allQ.first << ", " << avgblastFP_lowAlnLen_alnLen_allQ.second << endl;
        cmp << left << setw(80) << "(avg, median) BLAST's edit distance where alignment length < R: " << avgblastFP_lowAlnLen_ed_allQ.first << ", " << avgblastFP_lowAlnLen_ed_allQ.second << endl;
        
        // cmp << left << setw(80) << "# Of BLAST (FP) that edit pos breaks criteria: " << blastFP_editPos_allQ << endl;
        // pair<uint32_t, uint32_t> avgblastFP_editPos_alnLen_allQ = util.calculateStatistics2(blastFP_editPos_alnLen_allQ);
        // pair<uint32_t, uint32_t> avgblastFP_editPos_ed_allQ = util.calculateStatistics2(blastFP_editPos_ed_allQ);
        // cmp << left << setw(80) << "(avg, median) BLAST's alignment length where edit pos breaks criteria: " << avgblastFP_editPos_alnLen_allQ.first << ", " << avgblastFP_editPos_alnLen_allQ.second << endl;
        // cmp << left << setw(80) << "(avg, median) BLAST's edit distance where edit pos breaks criteria: " << avgblastFP_editPos_ed_allQ.first << ", " << avgblastFP_editPos_ed_allQ.second << endl;

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
        bestAlnCmp.close();
    }
};

#endif