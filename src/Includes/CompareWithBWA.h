#ifndef COMAPREWITHBWA_H
#define COMAPREWITHBWA_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>

using namespace std;

class ComparatorWithBWA {
public: 
    pair<uint32_t, uint32_t> countBwaClippedBases(const string& cigarString) {
        uint32_t softClippedStart = 0;
        uint32_t softClippedEnd = 0;
        uint32_t currentCount = 0;
        bool frontClip = true;
        for (char c : cigarString) {
            if (isdigit(c)) {
                currentCount = currentCount * 10 + (c - '0');
            } else if (c == 'S' || c == 'H') {
                if (frontClip) {
                    softClippedStart = currentCount;
                } else {
                    softClippedEnd = currentCount;
                }
                currentCount = 0;
            } else { // c == M, I, D etc
                frontClip = false;
                currentCount = 0;
            }
        }
        
        return make_pair(softClippedStart, softClippedEnd);
    }

    vector<int> editPositionsInBwaAlignment(string md, int frontClip)
    {
        vector<int> edits;
        int currentCount = 0;
        int cumCount = 0;
        for (char c : md) {
            if (isdigit(c)) {
                currentCount = currentCount * 10 + (c - '0');
            } else {
                currentCount++;
                cumCount += currentCount;
                edits.push_back(cumCount + frontClip);
                currentCount = 0;
            }
        }
        return edits;
    }

   vector<int> extractInDelLoc(string cigar, int frontClip)
    {
        vector<int> indels;
        int currentCount = 0, cumCount = frontClip;
        for (char c : cigar) {
            if (isdigit(c)) {
                currentCount = currentCount * 10 + (c - '0');
            } else {
                // cout << "cumCount : " << cumCount << endl;
                if (c == 'I' || c=='D')
                {    
                    uint32_t pos = cumCount;
                    for (int i = 0; i < currentCount; i++)
                    {
                        indels.push_back(pos);
                        pos++;
                    }
                    // cout << "indel in : " << cumCount + frontClip << endl;
                }
                cumCount += currentCount;
                currentCount = 0;
            }
        }
        return indels;
    }

    // bool hasMinConsecutiveMatches(SamReader::Sam bwaAln, Config cfg, set<int> editset, uint32_t frontClip)
    // {
    //     int currentCount = 0, cumCount = frontClip;
    //     for (char c : bwaAln.cigar) {
    //         if (isdigit(c)) {
    //             currentCount = currentCount * 10 + (c - '0');
    //         } else {
    //             // cout << "cumCount : " << cumCount << endl;
    //             if (c == 'M' && currentCount >= cfg.minExactMatchLen)
    //             {    
    //                 uint32_t startPos = cumCount;
    //                 uint32_t endPos = cumCount + currentCount - 1;
    //                 uint32_t exactMatchSize = 0;
    //                 // cout << "start : " << startPos << " end: " << endPos << endl;
    //                 for (int i = startPos; i <= endPos; i++)
    //                 {
    //                     if (i < cfg.regionSize || i >= cfg.contigSize - cfg.regionSize)
    //                     {
    //                         bool editfound = false;
    //                         for (auto j : edits)
    //                         {
    //                             // cout << "i : " << i << " j : " << j << endl;
    //                             if (i == j)
    //                             {
    //                                 // cout << "found" << endl;
    //                                 editfound = true;
    //                                 break;
    //                             }
    //                         }
    //                         if (!editfound)
    //                         {
    //                             exactMatchSize++;
    //                             if (exactMatchSize >= cfg.minExactMatchLen)
    //                             {
    //                                 // cout << "exactMatchSize : " << exactMatchSize << endl;
    //                                 return true;
    //                             }
    //                         } else
    //                             exactMatchSize = 0;
    //                     } else
    //                         exactMatchSize = 0;
    //                 }
    //             }
    //             cumCount += currentCount;
    //             currentCount = 0;
    //         }
    //     }
    //     return false;
    // }

    // bool checkBwaEditPositions(SamReader::Sam bwaAln, Config cfg, uint32_t bwaAlignmentLen)
    // {
    //     bool frontRegionDismissed = false, backRegionDismissed = false;
    //     pair<uint32_t, uint32_t> clips = countBwaClippedBases(bwaAln.cigar);
    //     vector<int> edits = editPositionsInBwaAlignment(bwaAln.mismatchPositions, clips.first);
    //     vector<int> Indels = extractInDelLoc(bwaAln.cigar, clips.first);
    //     uint32_t totalEdits = edits.size() + Indels.size();
    //     set<int> editset(edits.begin(), edits.end());
    //     editset.insert(Indels.begin(), Indels.end());
    //     uint32_t queryS = clips.first;

    //     //check whether the region starts at most from 100 + 2
    //     if(queryS > (cfg.contigSize - cfg.regionSize + cfg.editDistance)){ //first bp of back region starts from 100
    //         return false;
    //     }
    //     //check whether the regions has at least minExactMatchLen consecutive matches
    //     if(!hasMinConsecutiveMatches(bwaAln, cfg, editset, clips.first)){
    //         // cout << "!hasMinConsecutiveMatches" << endl;
    //         return false;
    //     }
    //     if(cfg.editDistance >= queryS){
    //         auto editsAllwedInFrontRegion = cfg.editDistance - queryS;
    //         if(bwaAln.editDistance + bwaAln.InDels > editsAllwedInFrontRegion)
    //             frontRegionDismissed = true;
    //     } else {
    //         frontRegionDismissed = true;
    //     }
    //     if(queryS + bwaAlignmentLen >= cfg.contigSize - cfg.editDistance){ //query start position + alignment len should be larger than conig size - allowed edit
    //         auto editsAllwedInBackRegion = queryS + bwaAlignmentLen - (cfg.contigSize - cfg.editDistance);
    //         if(totalEdits > editsAllwedInBackRegion)
    //             backRegionDismissed = true;
    //     }else {
    //         backRegionDismissed = true;
    //     }
    //     if(frontRegionDismissed && backRegionDismissed)
    //         return false;
    //     return true;
    // }

    // // bool checkBwaAlignmentClips(const string& cigarString, uint32_t maxEditAllowed)
    // // {
    // //     bool frontDismissed = false;//the alignment does not have min exact size front of query in it
    // //     bool backDismissed = false;//the alignment does not have min exact size back of query in it
    // //     pair<uint32_t, uint32_t> clips = countBwaClippedBases(cigarString);
    // //     // cout << "front clips: " << clips.first <<  ", back clips: " << clips.second << endl;
    // //     if (clips.first > maxEditAllowed)
    // //     {
    // //         frontDismissed = true;
    // //         // cout << "frontDismissed\n";
    // //     }
    // //     if (clips.second > maxEditAllowed)
    // //     {
    // //         backDismissed = true;
    // //         // cout << "backDismissed\n";
    // //     }
    // //     if (frontDismissed & backDismissed)
    // //         return false;
    // //     return true;
    // // }

    // // bool checkBwaAlignmentBasedOnOurCriteria(uint32_t maxEditAllowed, int regionSize, int minExactMatchSize, int contigSize, string cigar, string md, ofstream &cmp)
    // // {
    // //     //check if the clips covered the front and back region
    // //     if(checkBwaAlignmentClips(cigar, maxEditAllowed))
    // //     {
    // //         //check if clips + edits covered the front and back region
    // //         pair<int, int> clips = countBwaClippedBases(cigar);
    // //         vector<int> edits = editPositionsInBwaAlignment(md, clips.first);
    // //         vector<int>Indels = extractInDelLoc(cigar, clips.first);
    // //         bool firstKmerInFrontRegionDismissed = false, lastKmerInFrontRegionDismissed = false,
    // //         firstKmerInBackRegionDismissed = false, lastKmerInBacktRegionDismissed = false;
    // //         for(auto v : edits)
    // //         {
    // //             // cout << "loc : " << v << endl;
    // //             if((v >= 0 && v <= minExactMatchSize) || clips.first > 0)
    // //             {
    // //                 firstKmerInFrontRegionDismissed = true;
    // //                 // cout << "firstKmerInFrontRegionDismissed" << endl;
    // //             }
    // //             if((v >= regionSize - minExactMatchSize && v <= regionSize) || (clips.first > regionSize - minExactMatchSize))
    // //             {
    // //                 lastKmerInFrontRegionDismissed = true;
    // //                 // cout << "lastKmerInFrontRegionDismissed" << endl;
    // //             }
    // //             if((v >= contigSize - regionSize && v <= contigSize - regionSize + minExactMatchSize)  || (clips.second > regionSize - minExactMatchSize))
    // //             {
    // //                 firstKmerInBackRegionDismissed = true;
    // //                 // cout << "firstKmerInBackRegionDismissed" << endl;
    // //             }
    // //             if((v >= contigSize - minExactMatchSize && v <= contigSize) || clips.second > 0)
    // //             {
    // //                 lastKmerInBacktRegionDismissed = true;
    // //                 // cout << "lastKmerInBacktRegionDismissed" << endl;
    // //             }
    // //         }
    // //         for(auto v : Indels)
    // //         {
    // //             // cout << "indel loc : " << v << endl;
    // //             if(v >= 0 && v <= minExactMatchSize)
    // //             {
    // //                 firstKmerInFrontRegionDismissed = true;
    // //                 // cout << "firstKmerInFrontRegionDismissed" << endl;
    // //             }
    // //             if(v >= regionSize - minExactMatchSize && v <= regionSize)
    // //             {
    // //                 lastKmerInFrontRegionDismissed = true;
    // //                 // cout << "lastKmerInFrontRegionDismissed" << endl;
    // //             }
    // //             if((v >= contigSize - regionSize && v <= contigSize - regionSize + minExactMatchSize)  || (clips.second > regionSize - minExactMatchSize))
    // //             {
    // //                 firstKmerInBackRegionDismissed = true;
    // //                 // cout << "firstKmerInBackRegionDismissed" << endl;
    // //             }
    // //             if((v >= contigSize - minExactMatchSize && v <= contigSize) || clips.second > 0)
    // //             {
    // //                 lastKmerInBacktRegionDismissed = true;
    // //                 // cout << "lastKmerInBacktRegionDismissed" << endl;
    // //             }
    // //         }
    // //         if(firstKmerInFrontRegionDismissed && lastKmerInFrontRegionDismissed && firstKmerInBackRegionDismissed && lastKmerInBacktRegionDismissed)
    // //         {
    // //             cmp << "the read was supposed to be discarded based on our criteria (edit pos)";
    // //             return false;
    // //         }
    // //     } else 
    // //     {
    // //         cmp << "the read was supposed to be discarded based on our criteria (clips)";
    // //         return false;
    // //     }
    // //     return true;
    // // }

    // void comparePmWithBWA(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress, IndexContainer<uint32_t, LevAlign>& pmAlignments, vector<std::pair<uint32_t, uint32_t>>& alnPmBwaHisto, const uint32_t queryCount, string alnPerQueryFileAddress)
    // {
    //     ofstream cmp(comparisonResultsFileAddress);
    //     ofstream alnPerQ(alnPerQueryFileAddress);
    //     SamReader bwaSam(cfg.otherToolOutputFileAddress);
    //     vector<SamReader::Sam> bwaAlignments = bwaSam.parseFile(queryCount);
    //     uint64_t numberOfEqualalignment = 0, numberOfBWAbetter = 0, numberOfPMbetter = 0, 
    //     numberOfBWAbetterExceedMaxEdits = 0, numberOfBwaClippedAlignments = 0, 
    //     numberOfBWAbetterWithLowMatchSize = 0, numberOf BWAbetterNotObserveOurCriteria = 0, 
    //     onlyBwaFoundMatchForQuery = 0, bwaReadIdFoundInPM = 0, bwaReadIdNotFoundInPM = 0, 
    //     onlyPmFoundMatchForQuery = 0, nonePmBWAFoundAlignmentForQuery = 0, bwaClipCount = 0, 
    //     numberOfQueryContainN = 0, pmReadIdNotFoundInBWA = 0, totalNumberOfBWATruePositive = 0, totalNumberOfBWAFalsePositive = 0;
    //     LevAlign pmAlignment;
    //     uint32_t pmQueriesFound = 0, bwaQueriesFound = 0;
    //     SamReader::Sam bwaAlignment;
    //     for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
    //     {
    //         bool bwaFound = false, pmFound = false;
    //         uint32_t bwaMatchSize = 0, bwaInDels = 0, pmMatchSize = 0;
    //         string query = queries[queryInd];
    //         if(query.find('N') != string::npos || query.find('n') != string::npos)
    //         {
    //             cmp << "query contains N!" << endl;
    //             numberOfQueryContainN++; 
    //             alnPerQ << queryInd << " 0 0"<< endl;
    //             continue;
    //         }
    //         cmp << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    //         cmp << "Q : " << query << ", queryInd: " << queryInd << endl;

    //         for (const SamReader::Sam& aln : bwaAlignments) 
    //         {
    //             if(aln.queryId == queryInd)
    //             {
    //                 string bwaR = reads[bwaAlignment.readId];
    //                 if(bwaR.find('N') != string::npos || bwaR.find('n') != string::npos)
    //                     break;
    //                 bwaQueriesFound++;
    //                 bwaFound = true;
    //                 bwaAlignment = aln;
    //                 bwaMatchSize = bwaSam.countMatches(bwaAlignment.cigar);
    //                 bwaClipCount = bwaSam.countClips(bwaAlignment.cigar);
    //                 bwaInDels = bwaSam.countInsertions(bwaAlignment.cigar) + bwaSam.countDeletions(bwaAlignment.cigar);
    //                 if (bwaClipCount > 0)
    //                 {
    //                     numberOfBwaClippedAlignments++;
    //                 }
    //                 break;
    //             }
    //         }
    //         string bwaRead = "";
    //         if(bwaFound)
    //             bwaRead = reads[bwaAlignment.readId];
    //         auto range = pmAlignments.getRange(queryInd);
    //         size_t pmReadPerQuery = distance(range.first, range.second);
    //         size_t bwaReadPerQuery = bwaFound ? 1 : 0 ;
    //         cmp << "# of reads found by BWA : " << bwaReadPerQuery << endl;
    //         cmp << "# of reads found by PARMIK : " << pmReadPerQuery << endl;
    //         LevAlign bestAlnPm;
    //         for (auto it = range.first; it != range.second; it++) 
    //         {
    //             pmFound = true;
    //             LevAlign aln = it->second;
    //             if (aln.numberOfMatches + aln.numberOfInDel > bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel)
    //             {
    //                 bestAlnPm = aln;
    //             } else if (aln.numberOfMatches + aln.numberOfInDel == bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel)
    //             {
    //                 if (aln.numberOfInDel + aln.numberOfSub < bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub) // InDel has the same wight as substitution
    //                 {
    //                     bestAlnPm = aln;
    //                 }
    //             }
    //         }
    //         string pmRead = "";
    //         if (pmFound)
    //         {
    //             pmRead = reads[bestAlnPm.readID];
    //             pmQueriesFound++;
    //             pmMatchSize = bwaSam.countMatches(bestAlnPm.cigar);
    //         }
    //         alnPmBwaHisto.push_back(std::make_pair(pmMatchSize, bwaMatchSize));
    //         if(bwaReadPerQuery == 0 && pmReadPerQuery == 0){
    //             //TN
    //             totalBWATN++;
    //             cmp << "TN for BWA" << endl;
    //         } else if(bwaReadPerQuery == 0 && pmReadPerQuery > 0) {
    //             //FN
    //             bwaFN++;
    //             bwaFN_parmik_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
    //             bwaFN_parmik_ed_allQ.insert(bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub);
    //             cmp << "FN for BWA, PARMIK alnlen : " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", ed : " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", readID: " <<  bestAlnPm.readID << endl;
    //         } else {
    //             bool isFP = false;
    //             if (bwaAlignment.editDistance + bwaAlignment.InDels > cfg.editDistance)
    //             {
    //                 // cmp << "with edits (Indel+Subs) [" << bwaAlignment.editDistance + bwaInDels << "] > " << cfg.editDistance << endl;
    //                 isFP = true;
    //                 bwaFP_editsExceed++;
    //                 bwaFP_editsExceed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                 bwaFP_editsExceed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //             } else if(bwaMatchSize + bwaInDels < cfg.regionSize) {    
    //                 // cmp << "with low match size : " << bwaMatchSize << endl;
    //                 isFP = true;
    //                 bwaFP_lowAlnLen++;
    //                 bwaFP_lowAlnLen_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                 bwaFP_lowAlnLen_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //             } else if (!hasMinConsecutiveMatches(bwaAlignment, cfg, set<int> edits)checkBwaEditPositions(bwaAlignment, cfg, bwaMatchSize + bwaInDelsn)){
    //                 isFP = true;
    //                 bwaFP_editPos++;
    //                 bwaFP_editPos_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                 bwaFP_editPos_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                 // cmp << "the read was supposed to be discarded based on our criteria (edit pos)" << endl;
    //             } 
    //             if (isFP){
    //                 bwaFP++;
    //             } else {
    //                 bwaTP++;
    //                 bool pmfound = false;
    //                 LevAlign pmaln;
    //                 for (auto itt = pmRange.first; itt != pmRange.second; itt++) 
    //                 {
    //                     pmaln = itt->second;
    //                     if((uint32_t) pmaln.readID == A.readId)
    //                     {
    //                         pmfound = true;
    //                         break;
    //                     }
    //                 }
    //                 if (!pmfound)
    //                 {
    //                     bwaTP_noParmikMatches++;
    //                     bwaTP_noParmikMatches_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                     bwaTP_noParmikMatches_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                     cmp << "TP for BWA, no PARMIK match, alnlen : " << bwaMatchSize + bwaInDels << ", ed : " << bwaAlignment.editDistance + bwaInDels << ", readID: " <<  bwaAlignment.readId << endl;
    //                     pmFn << queryInd << "\t" <<  bwaAlignment.readId << "\t" <<  bwaAlignment.flag << endl;
    //                 } else {
    //                     if (pmaln.numberOfMatches + pmaln.numberOfInDel > bwaMatchSize + bwaInDels){ //PARMIK outperformed
    //                         bwaTP_parmikOutperfomed++;
    //                         bwaTP_parmikOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                         bwaTP_parmikOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                         parmikTP_parmikOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
    //                         parmikTP_parmikOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
    //                     } else if (pmaln.numberOfMatches + pmaln.numberOfInDel < bwaMatchSize + bwaInDels){//BWA outperformed
    //                         bwaTP_bwaOutperfomed++;
    //                         bwaTP_bwaOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                         bwaTP_bwaOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                         parmikTP_bwaOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
    //                         parmikTP_bwaOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
    //                         cmp << "TP for BWA, BWA outperformed in terms of alnlen, Bwa alnlen : " << bwaMatchSize + bwaInDels << ", ed : " << bwaAlignment.editDistance + bwaInDels << ", readID: " <<  aln.readId << ", PARMIK alnlen: " << pmaln.numberOfMatches + pmaln.numberOfInDel << ", ed: " << pmaln.numberOfSub + pmaln.numberOfInDel << ", readID: " << pmaln.readID << endl;
    //                     } else {
    //                         if (pmaln.numberOfSub + pmaln.numberOfInDel < bwaAlignment.editDistance + bwaInDels)//PARMIK outperformed
    //                         {
    //                             bwaTP_parmikOutperfomed++;
    //                             bwaTP_parmikOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                             bwaTP_parmikOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                             parmikTP_parmikOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
    //                             parmikTP_parmikOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
    //                         } else if (pmaln.numberOfSub + pmaln.numberOfInDel > bwaAlignment.editDistance + bwaInDels)//BWA outperformed
    //                         {
    //                             bwaTP_bwaOutperfomed++;
    //                             bwaTP_bwaOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                             bwaTP_bwaOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                             parmikTP_bwaOutperfomed_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
    //                             parmikTP_bwaOutperfomed_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
    //                             cmp << "TP for BWA, BWA outperformed in terms of ed, Bwa alnlen : " << bwaMatchSize + bwaInDels << ", ed : " << bwaAlignment.editDistance + bwaInDels << ", readID: " <<  aln.readId << ", PARMIK alnlen: " << pmaln.numberOfMatches + pmaln.numberOfInDel << ", ed: " << pmaln.numberOfSub + pmaln.numberOfInDel << ", readID: " << pmaln.readID << endl;
    //                         } // BWA and PARMIK performed equally
    //                         else{
    //                             bwaTP_bwaEqualParmik++;
    //                             bwaTP_bwaEqualParmik_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                             bwaTP_bwaEqualParmik_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                             parmikTP_bwaEqualParmik_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
    //                             parmikTP_bwaEqualParmik_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
    //                         }
    //                     }
    //                 }
    //             }
    //             for (auto itt = pmRange.first; itt != pmRange.second; itt++) 
    //             {
    //                 LevAlign aaln = itt->second;
    //                 bool bwafound = false;
    //                 if(bwaAlignment.readId == (uint32_t)aaln.readID)
    //                 {
    //                     bwafound = true;
    //                     break;
    //                 }
    //                 if(!bwafound)
    //                 {
    //                     parmikTP_noBwaMatches++;
    //                     parmikTP_noBwaMatches_alnLen_allQ.insert(pmaln.numberOfMatches + pmaln.numberOfInDel);
    //                     parmikTP_noBwaMatches_ed_allQ.insert(pmaln.numberOfSub + pmaln.numberOfInDel);
    //                 }
    //             }
    //             parmikTP_noBwaMatches_allQ += parmikTP_noBwaMatches;
    //             cmp << "parmik TP that bwa did not find: " << parmikTP_noBwatMatches << endl;
    //             cmp << "BwaTP: " << bwaTP << "\t" << "BwaFP: " << bwaFP << endl;
    //             bwaTPPerQuery.insert(bwaTP);
    //             totalBwaTP += bwaTP;
    //             bwaFPPerQuery.insert(bwaFP);
    //             totalBwaFP += bwaFP;
    //             if(bwaTP > 0){ 
    //                 queriesBWAFoundMatch++;
    //             }
    //             if(pmReadPerQuery > 0 && bwaReadPerQuery > 0 && bwaTP == 0){//PARMIK TP
    //                 bwaFN_noCriteriaFittedMatches++;
    //                 bwaFN_noCriteria_parmik_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
    //                 bwaFN_noCriteria_parmik_ed_allQ.insert(bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub);
    //                 cmp << "BwaFN_noCriteriaFitted and PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", PARMIK readID: " << bestAlnPm.readID << endl;
    //             }
    //             if(pmReadPerQuery == 0 && bwaTP > 0){//PARMIK FN
    //                 pmFN_bwa_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                 pmFN_bwa_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                 cmp << "PARMIK FN, and BWA alnsize: " << bwaMatchSize + bwaInDels << ", BWA ed: " << bwaAlignment.editDistance + bwaInDels << ", BWA readID: " << bestAlnBwa.readId << endl;
    //                 numnberOfPM_FN++;
    //             } else if(pmReadPerQuery > 0 && bwaTP > 0){//compare the best alignmnets for this queyry
    //                 if (bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel > bwaMatchSize + bwaInDels){ //PARMIK outperformed
    //                     bwaBest_parmikOutperfomed++;
    //                     bwaBest_parmikOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                     bwaBest_parmikOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                     parmikBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
    //                     parmikBest_parmikOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
    //                     cmp << "Best: PARMIK outperformed BWA, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub <<  ", PARMIK readID: " << bestAlnPm.readID << ", BWA alnsize: " << bwaMatchSize + bwaInDels << ", BWA ed: " << bwaAlignment.editDistance + bwaInDels <<  ", BWA readID: " << bestAlnBwa.readId << endl;
    //                 } else if (bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel < bwaMatchSize + bwaInDels){//BWA outperformed
    //                     bwaBest_bwaOutperfomed++;
    //                     bwaBest_bwaOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                     bwaBest_bwaOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                     parmikBest_bwaOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
    //                     parmikBest_bwaOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
    //                     cmp << "Best: BWA outperformed PARMIK in terms of alnlen, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub <<  ", PARMIK readID: " << bestAlnPm.readID << ", BWA alnsize: " << bwaMatchSize + bwaInDels << ", BWA ed: " << bwaAlignment.editDistance + bwaInDels << ", BWA readID: " << bestAlnBwa.readId << endl;
    //                 } else {
    //                     if (bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel < bwaAlignment.editDistance + bwaInDels)//PARMIK outperformed
    //                     {
    //                         bwaBest_parmikOutperfomed++;
    //                         bwaBest_parmikOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                         bwaBest_parmikOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                         parmikBest_parmikOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
    //                         parmikBest_parmikOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
    //                         cmp << "Best: PARMIK outperformed BWA, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", PARMIK readID: " << bestAlnPm.readID << ", BWA alnsize: " << bwaMatchSize + bwaInDels << ", BWA ed: " << bwaAlignment.editDistance + bwaInDels <<  ", BWA readID: " << bestAlnBwa.readId << endl;
    //                     } else if (bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel > bwaAlignment.editDistance + bwaInDels)//BWA outperformed
    //                     {
    //                         bwaBest_bwaOutperfomed++;
    //                         bwaBest_bwaOutperfomed_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                         bwaBest_bwaOutperfomed_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                         parmikBest_bwaOutperfomed_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
    //                         parmikBest_bwaOutperfomed_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
    //                         cmp << "Best: BWA outperformed PARMIK in terms of ed, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub <<  ", PARMIK readID: " << bestAlnPm.readID << ", BWA alnsize: " << bwaMatchSize + bwaInDels << ", BWA ed: " << bwaAlignment.editDistance + bwaInDels <<  ", BWA readID: " << bestAlnBwa.readId << endl;
    //                     } // BWA and PARMIK performed equally
    //                     else{
    //                         bwaBest_bwaEqualParmik++;
    //                         bwaBest_bwaEqualParmik_alnLen_allQ.insert(bwaMatchSize + bwaInDels);
    //                         bwaBest_bwaEqualParmik_ed_allQ.insert(bwaAlignment.editDistance + bwaInDels);
    //                         parmikBest_bwaEqualParmik_alnLen_allQ.insert(bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel);
    //                         parmikBest_bwaEqualParmik_ed_allQ.insert(bestAlnPm.numberOfSub + bestAlnPm.numberOfInDel);
    //                         cmp << "Best: BWA & PARMIK perfromed equally, PARMIK alnsize: " << bestAlnPm.numberOfMatches + bestAlnPm.numberOfInDel << ", PARMIK ed: " << bestAlnPm.numberOfInDel + bestAlnPm.numberOfSub << ", PARMIK readID: " << bestAlnPm.readID << ", BWA alnsize: " << bwaMatchSize + bwaInDels << ", BWA ed: " << bwaAlignment.editDistance + bwaInDels << ", BWA readID: " << bestAlnBwa.readId << endl;
    //                     }
    //                 }
    //             }
    //         }
    //         bwaTP_noParmikMatches_allQ += bwaTP_noParmikMatches; 
    //         bwaTP_bwaOutperfomed_allQ += bwaTP_bwaOutperfomed; 
    //         bwaTP_parmikOutperfomed_allQ += bwaTP_parmikOutperfomed;
    //         bwaTP_bwaEqualParmik_allQ += bwaTP_bwaEqualParmik;
    //         cmp << "bwaTP_noParmikMatches: " << bwaTP_noParmikMatches << ", bwaTP_bwaOutperfomed: " << bwaTP_bwaOutperfomed << ", bwaTP_parmikOutperfomed: " << bwaTP_parmikOutperfomed << ", bwaTP_bwaEqualParmik: " << bwaTP_bwaEqualParmik << endl;
    //         bwaFP_editsExceed_allQ += bwaFP_editsExceed;
    //         bwaFP_lowAlnLen_allQ += bwaFP_lowAlnLen; 
    //         bwaFP_editPos_allQ += bwaFP_editPos;
    //         cmp << "bwaFP_editsExceed: " << bwaFP_editsExceed << ", bwaFP_lowAlnLen: " << bwaFP_lowAlnLen << ", bwaFP_editPos: " << bwaFP_editPos << endl;
    //         bwaBest_parmikOutperfomed_allQ += bwaBest_parmikOutperfomed;
    //         bwaBest_bwaOutperfomed_allQ += bwaBest_bwaOutperfomed;
    //         bwaBest_bwaEqualParmik_allQ += bwaBest_bwaEqualParmik;
    //         cmp << "bwaBest_parmikOutperfomed: " << bwaBest_parmikOutperfomed << ", bwaBest_bwaOutperfomed: " << bwaBest_bwaOutperfomed << ", bwaBest_bwaEqualParmik: " << bwaBest_bwaEqualParmik << endl;
    //         alnPmBWAHisto.push_back(make_pair(pmMatchSize, bwaMatchSize)); 
    //         //PM and BWA alignments per query
    //         alnPerQ << queryInd << " " << pmAlignments.container_.count(queryInd) << " " << bwaTP << endl;
    //     }
    //     Utilities<uint32_t> util;   
    //     cmp << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Overall Comparison Results>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    //     cmp << left << setw(80) << "# Of queries : " << queryCount << endl;
    //     cmp << left << setw(80) << "# Of queries that PARMIK found match : " << queriesPMFoundMatch << endl;
    //     cmp << left << setw(80) << "# Of queries that BWA found match (TP) : " << queriesBWAFoundMatch << endl;
        
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of queries that BWA didn't find TN : " << totalBwaTN << endl;
        
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of queries that BWA didn't find FN (Total) : " << bwaFN + bwaFN_noCriteriaFittedMatches << endl;
    //     cmp << left << setw(80) << "# Of queries that BWA didn't find any match FN : " << bwaFN << endl;
    //     pair<uint32_t, uint32_t> avgBwaFN_parmik_alnLen_allQ = util.calculateStatistics2(bwaFN_parmik_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgBwaFN_parmik_ed_allQ = util.calculateStatistics2(bwaFN_parmik_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BWA FN: " << avgBwaFN_parmik_alnLen_allQ.first << ", " << avgBwaFN_parmik_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BWA FN: " << avgBwaFN_parmik_ed_allQ.first << ", " << avgBwaFN_parmik_ed_allQ.second << endl;
    //     cmp << left << setw(80) << "# Of queries that BWA didn't find based on criteria FN : " << bwaFN_noCriteriaFittedMatches << endl;
    //     pair<uint32_t, uint32_t> avgbwaFN_noCriteria_parmik_alnLen_allQ = util.calculateStatistics2(bwaFN_noCriteria_parmik_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaFN_noCriteria_parmik_ed_allQ = util.calculateStatistics2(bwaFN_noCriteria_parmik_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BWA FN based on criteria: " << avgbwaFN_noCriteria_parmik_alnLen_allQ.first << ", " << avgbwaFN_noCriteria_parmik_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BWA FN based on criteria: " << avgbwaFN_noCriteria_parmik_ed_allQ.first << ", " << avgbwaFN_noCriteria_parmik_ed_allQ.second << endl;
       
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of queries that PARMIK didn't find FN : " << numnberOfPM_FN << endl;
    //     pair<uint32_t, uint32_t> avgPmFN_bwa_alnLen_allQ = util.calculateStatistics2(pmFN_bwa_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgPmFN_bwa_ed_allQ = util.calculateStatistics2(pmFN_bwa_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where PARMIK FN: " << avgPmFN_bwa_alnLen_allQ.first << ", " << avgPmFN_bwa_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where PARMIK FN: " << avgPmFN_bwa_ed_allQ.first << ", " << avgPmFN_bwa_ed_allQ.second << endl;
       
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of readIDs found by PARMIK (total): " << pmTotalNumberOfReadIDs << endl;
    //     cmp << left << setw(80) << "# Of readIDs found by BWA (total): " << bwaTotalNumberOfReadIDs << endl;
    //     cmp << left << setw(80) << "# Of readIDs found by BWA containing N (total): " << numberOfBwaReadsContainingN << endl;
        
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of total BWA (TP): " << totalBwaTP << endl;
    //     pair<uint32_t, uint32_t> avgbwaTPPerQuery = util.calculateStatistics2(bwaTPPerQuery);
    //     cmp << left << setw(80) << "(avg, median) BWA's TP per query: " << avgbwaTPPerQuery.first << ", " << avgbwaTPPerQuery.second << endl;

    //     cmp << left << setw(80) << "# Of BWA (TP) where PARMIK didn't find match: " << bwaTP_noParmikMatches_allQ << endl;
    //     pair<uint32_t, uint32_t> avgBwaTP_noParmikMatches_alnLen_allQ = util.calculateStatistics2(bwaTP_noParmikMatches_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgBwaTP_noParmikMatches_ed_allQ = util.calculateStatistics2(bwaTP_noParmikMatches_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where PARMIK didn't find match: " << avgBwaTP_noParmikMatches_alnLen_allQ.first << ", " << avgBwaTP_noParmikMatches_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where PARMIK didn't find match: " << avgBwaTP_noParmikMatches_ed_allQ.first << ", " << avgBwaTP_noParmikMatches_ed_allQ.second << endl;
        
    //     cmp << left << setw(80) << "# Of BWA (TP) where PARMIK outperformed: " << bwaTP_parmikOutperfomed_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaTP_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(bwaTP_parmikOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaTP_parmikOutperfomed_ed_allQ = util.calculateStatistics2(bwaTP_parmikOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where PARMIK outperformed: " << avgbwaTP_parmikOutperfomed_alnLen_allQ.first << ", " << avgbwaTP_parmikOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where PARMIK outperformed: " << avgbwaTP_parmikOutperfomed_ed_allQ.first << ", " << avgbwaTP_parmikOutperfomed_ed_allQ.second << endl;
    //     pair<uint32_t, uint32_t> avgparmikTP_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikTP_parmikOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgparmikTP_parmikOutperfomed_ed_allQ = util.calculateStatistics2(parmikTP_parmikOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where PARMIK outperformed: " << avgparmikTP_parmikOutperfomed_alnLen_allQ.first << ", " << avgparmikTP_parmikOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where PARMIK outperformed: " << avgparmikTP_parmikOutperfomed_ed_allQ.first << ", " << avgparmikTP_parmikOutperfomed_ed_allQ.second << endl;

    //     cmp << left << setw(80) << "# Of BWA (TP) where BWA outperformed: " << bwaTP_bwaOutperfomed_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaTP_bwaOutperfomed_alnLen_allQ = util.calculateStatistics2(bwaTP_bwaOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaTP_bwaOutperfomed_ed_allQ = util.calculateStatistics2(bwaTP_bwaOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where BWA outperformed: " << avgbwaTP_bwaOutperfomed_alnLen_allQ.first << ", " << avgbwaTP_bwaOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where BWA outperformed: " << avgbwaTP_bwaOutperfomed_ed_allQ.first << ", " << avgbwaTP_bwaOutperfomed_ed_allQ.second << endl;
    //     pair<uint32_t, uint32_t> avgparmikTP_bwaOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikTP_bwaOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgparmikTP_bwaOutperfomed_ed_allQ = util.calculateStatistics2(parmikTP_bwaOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BWA outperformed: " << avgparmikTP_bwaOutperfomed_alnLen_allQ.first << ", " << avgparmikTP_bwaOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BWA outperformed: " << avgparmikTP_bwaOutperfomed_ed_allQ.first << ", " << avgparmikTP_bwaOutperfomed_ed_allQ.second << endl;

    //     cmp << left << setw(80) << "# Of BWA (TP) where they perfromed equaly: " << bwaTP_bwaEqualParmik_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaTP_bwaEqualParmik_alnLen_allQ = util.calculateStatistics2(bwaTP_bwaEqualParmik_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaTP_bwaEqualParmik_ed_allQ = util.calculateStatistics2(bwaTP_bwaEqualParmik_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where they perfromed equaly: " << avgbwaTP_bwaEqualParmik_alnLen_allQ.first << ", " << avgbwaTP_bwaEqualParmik_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where they perfromed equaly: " << avgbwaTP_bwaEqualParmik_ed_allQ.first << ", " << avgbwaTP_bwaEqualParmik_ed_allQ.second << endl;
    //     pair<uint32_t, uint32_t> avgparmikTP_bwaEqualParmik_alnLen_allQ = util.calculateStatistics2(parmikTP_bwaEqualParmik_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgparmikTP_bwaEqualParmik_ed_allQ = util.calculateStatistics2(parmikTP_bwaEqualParmik_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where they perfromed equaly: " << avgparmikTP_bwaEqualParmik_alnLen_allQ.first << ", " << avgparmikTP_bwaEqualParmik_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where they perfromed equaly: " << avgparmikTP_bwaEqualParmik_ed_allQ.first << ", " << avgparmikTP_bwaEqualParmik_ed_allQ.second << endl;
        
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of BWA best where PARMIK outperformed: " << bwaBest_parmikOutperfomed_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaBest_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(bwaBest_parmikOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaBest_parmikOutperfomed_ed_allQ = util.calculateStatistics2(bwaBest_parmikOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's best alignment length where PARMIK outperformed: " << avgbwaBest_parmikOutperfomed_alnLen_allQ.first << ", " << avgbwaBest_parmikOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's best edit distance where PARMIK outperformed: " << avgbwaBest_parmikOutperfomed_ed_allQ.first << ", " << avgbwaBest_parmikOutperfomed_ed_allQ.second << endl;
    //     pair<uint32_t, uint32_t> avgparmikBest_parmikOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikBest_parmikOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgparmikBest_parmikOutperfomed_ed_allQ = util.calculateStatistics2(parmikBest_parmikOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's best alignment length where PARMIK outperformed: " << avgparmikBest_parmikOutperfomed_alnLen_allQ.first << ", " << avgparmikBest_parmikOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's best edit distance where PARMIK outperformed: " << avgparmikBest_parmikOutperfomed_ed_allQ.first << ", " << avgparmikBest_parmikOutperfomed_ed_allQ.second << endl;

    //     cmp << left << setw(80) << "# Of BWA best where BWA outperformed: " << bwaBest_bwaOutperfomed_allQ  << endl;
    //     pair<uint32_t, uint32_t> avgbwaBest_bwaOutperfomed_alnLen_allQ = util.calculateStatistics2(bwaBest_bwaOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaBest_bwaOutperfomed_ed_allQ = util.calculateStatistics2(bwaBest_bwaOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's best alignment length where BWA outperformed: " << avgbwaBest_bwaOutperfomed_alnLen_allQ.first << ", " << avgbwaBest_bwaOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's best edit distance where BWA outperformed: " << avgbwaBest_bwaOutperfomed_ed_allQ.first << ", " << avgbwaBest_bwaOutperfomed_ed_allQ.second << endl;
    //     pair<uint32_t, uint32_t> avgparmikBest_bwaOutperfomed_alnLen_allQ = util.calculateStatistics2(parmikBest_bwaOutperfomed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgparmikBest_bwaOutperfomed_ed_allQ = util.calculateStatistics2(parmikBest_bwaOutperfomed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's best alignment length where BWA outperformed: " << avgparmikBest_bwaOutperfomed_alnLen_allQ.first << ", " << avgparmikBest_bwaOutperfomed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's best edit distance where BWA outperformed: " << avgparmikBest_bwaOutperfomed_ed_allQ.first << ", " << avgparmikBest_bwaOutperfomed_ed_allQ.second << endl;

    //     cmp << left << setw(80) << "# Of BWA best where they perfromed equaly: " << bwaBest_bwaEqualParmik_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaBest_bwaEqualParmik_alnLen_allQ = util.calculateStatistics2(bwaBest_bwaEqualParmik_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaBest_bwaEqualParmik_ed_allQ = util.calculateStatistics2(bwaBest_bwaEqualParmik_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's best alignment length where they perfromed equaly: " << avgbwaBest_bwaEqualParmik_alnLen_allQ.first << ", " << avgbwaBest_bwaEqualParmik_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's best edit distance where they perfromed equaly: " << avgbwaBest_bwaEqualParmik_ed_allQ.first << ", " << avgbwaBest_bwaEqualParmik_ed_allQ.second << endl;
    //     pair<uint32_t, uint32_t> avgparmikBest_bwaEqualParmik_alnLen_allQ = util.calculateStatistics2(parmikBest_bwaEqualParmik_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgparmikBest_bwaEqualParmik_ed_allQ = util.calculateStatistics2(parmikBest_bwaEqualParmik_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's best alignment length where they perfromed equaly: " << avgparmikBest_bwaEqualParmik_alnLen_allQ.first << ", " << avgparmikBest_bwaEqualParmik_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's best edit distance where they perfromed equaly: " << avgparmikBest_bwaEqualParmik_ed_allQ.first << ", " << avgparmikBest_bwaEqualParmik_ed_allQ.second << endl;
        
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of total BWA (FP): " << totalBwaFP << endl;
    //     pair<uint32_t, uint32_t> avgbwaFPPerQuery = util.calculateStatistics2(bwaFPPerQuery);
    //     cmp << left << setw(80) << "(avg, median) BWA's FP per query: " << avgbwaFPPerQuery.first << ", " << avgbwaFPPerQuery.second << endl;

    //     cmp << left << setw(80) << "# Of BWA (FP) that exceed max edit distance: " << bwaFP_editsExceed_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaFP_editsExceed_alnLen_allQ = util.calculateStatistics2(bwaFP_editsExceed_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaFP_editsExceed_ed_allQ = util.calculateStatistics2(bwaFP_editsExceed_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where exceed max edit distance: " << avgbwaFP_editsExceed_alnLen_allQ.first << ", " << avgbwaFP_editsExceed_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where exceed max edit distance: " << avgbwaFP_editsExceed_ed_allQ.first << ", " << avgbwaFP_editsExceed_ed_allQ.second << endl;
        
    //     cmp << left << setw(80) << "# Of BWA (FP) that alignment length < R: " << bwaFP_lowAlnLen_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaFP_lowAlnLen_alnLen_allQ = util.calculateStatistics2(bwaFP_lowAlnLen_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaFP_lowAlnLen_ed_allQ = util.calculateStatistics2(bwaFP_lowAlnLen_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where alignment length < R: " << avgbwaFP_lowAlnLen_alnLen_allQ.first << ", " << avgbwaFP_lowAlnLen_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where alignment length < R: " << avgbwaFP_lowAlnLen_ed_allQ.first << ", " << avgbwaFP_lowAlnLen_ed_allQ.second << endl;
        
    //     cmp << left << setw(80) << "# Of BWA (FP) that edit pos breaks criteria: " << bwaFP_editPos_allQ << endl;
    //     pair<uint32_t, uint32_t> avgbwaFP_editPos_alnLen_allQ = util.calculateStatistics2(bwaFP_editPos_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgbwaFP_editPos_ed_allQ = util.calculateStatistics2(bwaFP_editPos_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) BWA's alignment length where edit pos breaks criteria: " << avgbwaFP_editPos_alnLen_allQ.first << ", " << avgbwaFP_editPos_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's edit distance where edit pos breaks criteria: " << avgbwaFP_editPos_ed_allQ.first << ", " << avgbwaFP_editPos_ed_allQ.second << endl;

    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of PARMIK (TP) where BWA didn't find match: " << parmikTP_noBwaMatches_allQ << endl;
    //     pair<uint32_t, uint32_t> avgparmikTP_noBwaMatches_alnLen_allQ = util.calculateStatistics2(parmikTP_noBwaMatches_alnLen_allQ);
    //     pair<uint32_t, uint32_t> avgparmikTP_noBwaMatches_ed_allQ = util.calculateStatistics2(parmikTP_noBwaMatches_ed_allQ);
    //     cmp << left << setw(80) << "(avg, median) PARMIK's alignment length where BWA didn't find match: " << avgparmikTP_noBwaMatches_alnLen_allQ.first << ", " << avgparmikTP_noBwaMatches_alnLen_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) PARMIK's edit distance where BWA didn't find match: " << avgparmikTP_noBwaMatches_ed_allQ.first << ", " << avgparmikTP_noBwaMatches_ed_allQ.second << endl;
        
    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     cmp << left << setw(80) << "# Of queries contains N: " << numberOfQueryContainN << endl;

    //     cmp << "------------------------------------------------------------------------------------------" << endl;
    //     pair<uint32_t, uint32_t> pmReadPerQuery_allQ  = util.calculateStatistics2(pmReadPerQuerySet);
    //     pair<uint32_t, uint32_t> bwaReadPerQuery_allQ  = util.calculateStatistics2(bwaReadPerQuerySet);
    //     cmp << left << setw(80) << "(avg, median) PARMIKS's Matches per Query: " << pmReadPerQuery_allQ.first << ", " << pmReadPerQuery_allQ.second << endl;
    //     cmp << left << setw(80) << "(avg, median) BWA's Matches per Query: " << bwaReadPerQuery_allQ.first << ", " << bwaReadPerQuery_allQ.second << endl;

    //     cmp.close();
    //     alnPerQ.close();
    // }
};

#endif