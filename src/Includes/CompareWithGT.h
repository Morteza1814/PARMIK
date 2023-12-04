#ifndef COMAPREWITHGT_H
#define COMAPREWITHGT_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>

using namespace std;

class CompareWithGroundTruth {
public: 

    void comparePmWithGroundTruth(IndexContainer<uint32_t, LevAlign> &gtAlignments, const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress, IndexContainer<uint32_t, LevAlign>& pmAlignments, vector<std::pair<uint32_t, uint32_t>>& alnPmGTHisto, const uint32_t queryCount, string alnPerQueryFileAddress)
    {
        ofstream cmp(comparisonResultsFileAddress);
        ofstream alnPerQ(alnPerQueryFileAddress);
        SamReader sam(cfg.otherToolOutputFileAddress);
        uint64_t numberOfEqualalignment = 0, numberOfGTbetter = 0, numberOfPMbetter = 0, 
        numberOfGTbetterExceedMaxEdits = 0, numberOfGTbetterWithLowMatchSize = 0,  
        onlyGTFoundMatchForQuery = 0, pmReadIdNotFoundInGT = 0, gtReadIdNotFoundInPM = 0, 
        onlyPmFoundMatchForQuery = 0, nonePmGTFoundAlignmentForQuery = 0, numberOfGTbetterNotObserveOurCriteria = 0, 
        numberOfQueryContainN = 0, queriesGTFoundMatch = 0, pmTotalNumberOfReadIDs = 0, gtTotalNumberOfReadIDs = 0;
        set<uint32_t> pmReadPerQuerySet, gtReadPerQuerySet;
        uint32_t pmQueriesFound = 0;
        for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            bool gtFound = false, pmFound = false;
            uint32_t gtMatchSize = 0, pmMatchSize = 0;
            string query = queries[queryInd];
            if(query.find('N') != string::npos || query.find('n') != string::npos)
            {
                cmp << "query contains N!" << endl;
                numberOfQueryContainN++; 
                alnPerQ << queryInd << " 0 0"<< endl;
                continue;
            }
            //PM and GT alignments per query
            alnPerQ << queryInd << " " << pmAlignments.container_.count(queryInd) << " " << gtAlignments.container_.count(queryInd) << endl;
            LevAlign gtAlignment;
            auto gtRrange = gtAlignments.getRange(queryInd);
            size_t gtReadPerQuery = distance(gtRrange.first, gtRrange.second);
            gtTotalNumberOfReadIDs += gtReadPerQuery;
            gtReadPerQuerySet.insert(gtReadPerQuery);
            for (auto it = gtRrange.first; it != gtRrange.second; it++) 
            {
                gtFound = true;
                LevAlign aln = it->second;
                if (aln.numberOfMatches + aln.numberOfInDel > gtAlignment.numberOfMatches + gtAlignment.numberOfInDel)
                {
                    gtAlignment = aln;
                } else if (aln.numberOfMatches + aln.numberOfInDel == gtAlignment.numberOfMatches + gtAlignment.numberOfInDel)
                {
                    if (aln.numberOfInDel + aln.numberOfSub < gtAlignment.numberOfInDel + gtAlignment.numberOfSub) // InDel has the same wight as substitution
                    {
                        gtAlignment = aln;
                    }
                }
            }
            string gtRead = reads[gtAlignment.readID];
            auto pmRange = pmAlignments.getRange(queryInd);
            size_t pmReadPerQuery = distance(pmRange.first, pmRange.second);
            pmTotalNumberOfReadIDs += pmReadPerQuery;
            pmReadPerQuerySet.insert(pmReadPerQuery);
            LevAlign pmAlignment;
            for (auto it = pmRange.first; it != pmRange.second; it++) 
            {
                pmFound = true;
                LevAlign aln = it->second;
                if (aln.numberOfMatches + aln.numberOfInDel > pmAlignment.numberOfMatches + pmAlignment.numberOfInDel)
                {
                    pmAlignment = aln;
                } else if (aln.numberOfMatches + aln.numberOfInDel == pmAlignment.numberOfMatches + pmAlignment.numberOfInDel)
                {
                    if (aln.numberOfInDel + aln.numberOfSub < pmAlignment.numberOfInDel + pmAlignment.numberOfSub) // InDel has the same wight as substitution
                    {
                        pmAlignment = aln;
                    }
                }
            }
            if (pmFound)
            {
                pmQueriesFound++;
                pmMatchSize = sam.countMatches(pmAlignment.cigar);
            }
            if (gtFound)
            {
                queriesGTFoundMatch++;
                gtMatchSize = sam.countMatches(gtAlignment.cigar);
            }
            alnPmGTHisto.push_back(std::make_pair(pmMatchSize, gtMatchSize));
            if (pmFound && !gtFound)
            {
                cmp << "PM outperformed and GT did not found" << endl; 
                numberOfPMbetter++;
                onlyPmFoundMatchForQuery++;
                cmp << "<<<<<<<<<<<PM alignment>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[pmAlignment.readID] << endl;
                // cmp << "q : " << pmAlignment.alignedQuery << endl;
                // cmp << "r : " << pmAlignment.alignedRead << endl;
                // cmp << "E : " << pmAlignment.editDistanceTypes << endl;
                cmp << "CIGAR : " << pmAlignment.cigar << ", subs : " << pmAlignment.numberOfSub << ", query id : " << queryInd << ", read id : " << pmAlignment.readID << ", flag : " << pmAlignment.flag << ", partial match size : " << pmAlignment.partialMatchSize << ", edits : " << pmAlignment.editDistance << ", R per Q : " << pmReadPerQuery << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            } else if (!pmFound && gtFound)
            {
                cmp << "GT outperformed and PM did not found"; 
                numberOfGTbetter++;
                onlyGTFoundMatchForQuery++;
                // if (gtAlignment.Mismatches + gtAlignment.InDels > cfg.editDistance)
                // {
                //     cmp << "with more edits > " << cfg.editDistance << endl;
                //     numberOfGTbetterExceedMaxEdits++;
                // } else if(gtMatchSize < cfg.regionSize) {
                //     cmp << "with low match size : " << gtMatchSize << endl;
                //     numberOfGTbetterWithLowMatchSize++;
                // } else {
                //     if (!checkBlastEditPositions(gtAlignment, query, gtRead, cfg)){
                //         numberOfGTbetterNotObserveOurCriteria++;
                //         cmp << "the read was supposed to be discarded based on our criteria (edit pos)" << endl;
                //     }
                // }
                cmp << "<<<<<<<<<<<GT alignment>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << gtRead << endl;
                cmp << "CIGAR : " << gtAlignment.cigar << ", subs : " << gtAlignment.numberOfSub << ", query id : " << queryInd << ", read id : " << gtAlignment.readID << ", flag : " << gtAlignment.flag << ", partial match size : " << gtAlignment.partialMatchSize << ", edits : " << gtAlignment.editDistance << ", R per Q : " << gtReadPerQuery << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            
            } else if (pmFound && gtFound) 
            {
                // int bwaInDelSize = bwaSam.countInsertions(gtAlignment.cigar) + bwaSam.countDeletions(gtAlignment.cigar);
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
               
                if (gtMatchSize + gtAlignment.numberOfInDel < pmMatchSize + pmAlignment.numberOfInDel)
                {
                    cmp << "PM outperformed "; 
                    numberOfPMbetter++;
                } else if (gtMatchSize + gtAlignment.numberOfInDel > pmMatchSize + pmAlignment.numberOfInDel)
                {
                    cmp << "GT outperformed";
                    numberOfGTbetter++;
                    // if (gtAlignment.Mismatches + gtAlignment.InDels > cfg.editDistance)
                    // {
                    //     cmp << "with more edits > " << cfg.editDistance << endl;
                    //     numberOfGTbetterExceedMaxEdits++;
                    // } else if(gtMatchSize < cfg.regionSize) {
                    //     cmp << "with low match size : " << gtMatchSize << endl;
                    //     numberOfGTbetterWithLowMatchSize++;
                    // } else {
                    //     if (!checkBlastEditPositions(gtAlignment, query, gtRead, cfg)){
                    //         numberOfGTbetterNotObserveOurCriteria++;
                    //         cmp << "the read was supposed to be discarded based on our criteria (edit pos)" << endl;
                    //     }
                    // }
                } else
                {
                    if (pmAlignment.numberOfSub + pmAlignment.numberOfInDel < gtAlignment.numberOfSub + gtAlignment.numberOfInDel)
                    {
                        cmp << "PM outperformed with fewer edits" << endl;
                        numberOfPMbetter++;
                    } else if (pmAlignment.numberOfSub + pmAlignment.numberOfInDel > gtAlignment.numberOfSub + gtAlignment.numberOfInDel)
                    {
                        cmp << "GT outperformed with fewer edits" << endl;
                        numberOfGTbetter++;
                        // if (gtAlignment.numberOfSub + gtAlignment.numberOfInDel > cfg.editDistance)
                        // {
                        //     numberOfGTbetterExceedMaxEdits++;
                        // }  
                    } else
                    {
                        cmp << "GT and PM performed equal" << endl;   
                        numberOfEqualalignment++;
                    }
                }
                //check for the read id
                bool found = false;
                for (auto it = pmRange.first; it != pmRange.second; it++) 
                {
                    LevAlign aln = it->second;
                    found = false;
                    for (auto itt = gtRrange.first; itt != gtRrange.second; itt++) 
                    {
                        LevAlign gt = itt->second;
                        if((uint32_t) aln.readID == (uint32_t)  gt.readID)
                        {
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        pmReadIdNotFoundInGT++;
                        // cmp << "pmReadIdNotFoundInGT = Q :" << gtAlignment.queryId << ", R : " << gtAlignment.readId << endl;
                    }
                }
                for (auto itt = gtRrange.first; itt != gtRrange.second; itt++) 
                {
                    LevAlign gt = itt->second;
                    string gtRead = reads[gt.readID];
                    found = false;
                    if(gtRead.find('N') == string::npos && gtRead.find('n') == string::npos)// gt read does not contain 'N'
                    {
                        for (auto it = pmRange.first; it != pmRange.second; it++) 
                        {
                            LevAlign aln = it->second;
                            if((uint32_t) aln.readID == (uint32_t) gt.readID)
                            {
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                        {
                            gtReadIdNotFoundInPM++;
                            // cmp << "pmReadIdNotFoundInGT = Q :" << gtAlignment.queryId << ", R : " << gtAlignment.readId << endl;
                        }
                    }
                }
                cmp << "<<<<<<<<<<<PM alignment>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[pmAlignment.readID] << endl;
                // cmp << "q : " << pmAlignment.alignedQuery << endl;
                // cmp << "r : " << pmAlignment.alignedRead << endl;
                // cmp << "E : " << pmAlignment.editDistanceTypes << endl;
                cmp << "CIGAR : " << pmAlignment.cigar << ", subs : " << pmAlignment.numberOfSub << ", query id : " << queryInd << ", read id : " << pmAlignment.readID << ", flag : " << pmAlignment.flag << ", partial match size : " << pmAlignment.partialMatchSize << ", edits : " << pmAlignment.editDistance << ", R per Q : " << pmReadPerQuery << endl;
                cmp << "<<<<<<<<<<<<<GT alignment>>>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << gtRead << endl;
                cmp << "CIGAR : " << gtAlignment.cigar << ", subs : " << gtAlignment.numberOfSub << ", query id : " << queryInd << ", read id : " << gtAlignment.readID << ", flag : " << gtAlignment.flag << ", partial match size : " << gtAlignment.partialMatchSize << ", edits : " << gtAlignment.editDistance << ", R per Q : " << gtReadPerQuery << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";  
            } else {
                cmp << "query " << queryInd << " not found in the PM results nor in GT" << endl;
                nonePmGTFoundAlignmentForQuery++;
            }
        }
        Utilities<uint32_t> util;   
        tuple<uint32_t, uint32_t, uint32_t> pmReadPerQueryTuple = util.calculateStatistics(pmReadPerQuerySet);
        tuple<uint32_t, uint32_t, uint32_t> gtReadPerQueryTuple = util.calculateStatistics(gtReadPerQuerySet);
        cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Overall Comparison Results>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cmp << left << setw(80) << "# Of queries : " << queryCount << endl;
        cmp << left << setw(80) << "# Of queries that PARMIK found match : " << pmQueriesFound << endl;
        cmp << left << setw(80) << "# Of queries that GT found match : " << queriesGTFoundMatch << endl;
        cmp << left << setw(80) << "# Of Matched Hits between PM and GT : " << numberOfEqualalignment << endl;
        cmp << left << setw(80) << "# Of readIDs found by PARMIK (total): " << pmTotalNumberOfReadIDs << endl;
        cmp << left << setw(80) << "# Of readIDs found by GT (total): " << gtTotalNumberOfReadIDs << endl;
        cmp << left << setw(80) << "# Of PM outperformed : " << numberOfPMbetter << endl;
        cmp << left << setw(80) << "# Of GT outperformed : " << numberOfGTbetter << endl;
        cmp << left << setw(80) << "# Of GT outperformed that Exceed Max Edits : " << numberOfGTbetterExceedMaxEdits << endl;
        cmp << left << setw(80) << "# Of GT outperformed that matchsize < (R - E): " << numberOfGTbetterWithLowMatchSize << endl;
        cmp << left << setw(80) << "# Of GT outperformed that has edits in minExact-mers of region : " << numberOfGTbetterNotObserveOurCriteria << endl;
        cmp << left << setw(80) << "# Of queries that GT found matches, not PM : " << onlyGTFoundMatchForQuery << endl;
        cmp << left << setw(80) << "# Of queries that PM found matches, not GT : " << onlyPmFoundMatchForQuery << endl;
        cmp << left << setw(80) << "# Of readIDs that PM found but GT did not : " << pmReadIdNotFoundInGT << endl;
        cmp << left << setw(80) << "# Of readIDs that GT found but PM did not : " << gtReadIdNotFoundInPM << endl;
        cmp << left << setw(80) << "# Of queries that neigther GT nor PM found any match : " << nonePmGTFoundAlignmentForQuery << endl;
        cmp << left << setw(80) << "# Of queries contains N: " << numberOfQueryContainN << endl;
        cmp << left << setw(80) << "# of Read Per Query found by PM" << "average: " <<  get<0>(pmReadPerQueryTuple) << ", median: " <<  get<1>(pmReadPerQueryTuple) << ", mean: " << get<2>(pmReadPerQueryTuple)<< endl;
        cmp << left << setw(80) << "# of Read Per Query found by GT" << "average: " <<  get<0>(gtReadPerQueryTuple) << ", median: " <<  get<1>(gtReadPerQueryTuple) << ", mean: " << get<2>(gtReadPerQueryTuple)<< endl;
        cmp.close();
        alnPerQ.close();
    }
};

#endif