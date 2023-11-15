#ifndef COMAPREWITHBWA_H
#define COMAPREWITHBWA_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>

using namespace std;

class ComparatorWithBWA {
public: 
    pair<int, int> countBwaClippedBases(const string& cigarString) {
        int softClippedStart = 0;
        int softClippedEnd = 0;
        int currentCount = 0;
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
        int currentCount = 0, cumCount = 0;
        for (char c : cigar) {
            if (isdigit(c)) {
                currentCount = currentCount * 10 + (c - '0');
            } else {
                currentCount++;
                cumCount += currentCount;
                if (c == 'I' || c=='D')
                {    
                    indels.push_back(cumCount + frontClip);
                    // cout << "indel in : " << cumCount + frontClip << endl;
                }
                currentCount = 0;
            }
        }
        return indels;
    }

    bool checkBwaAlignmentClips(const string& cigarString, int regionSize, int minExactMatchSize)
    {
        bool frontDismissed = false;//the alignment does not have min exact size front of query in it
        bool backDismissed = false;//the alignment does not have min exact size back of query in it
        pair<int, int> clips = countBwaClippedBases(cigarString);
        // cout << "front clips: " << clips.first <<  ", back clips: " << clips.second << endl;
        if (clips.first > regionSize-minExactMatchSize)
        {
            frontDismissed = true;
            // cout << "frontDismissed\n";
        }
        if (clips.second > regionSize-minExactMatchSize)
        {
            backDismissed = true;
            // cout << "backDismissed\n";
        }
        if (frontDismissed & backDismissed)
            return false;
        return true;
    }

    bool checkBwaAlignmentBasedOnOurCriteria(int regionSize, int minExactMatchSize, int contigSize, string cigar, string md, ofstream &cmp)
    {
        //check if the clips covered the front and back region
        if(checkBwaAlignmentClips(cigar, regionSize, minExactMatchSize))
        {
            //check if clips + edits covered the front and back region
            pair<int, int> clips = countBwaClippedBases(cigar);
            vector<int> edits = editPositionsInBwaAlignment(md, clips.first);
            vector<int>Indels = extractInDelLoc(cigar, clips.first);
            bool firstKmerInFrontRegionDismissed = false, lastKmerInFrontRegionDismissed = false,
            firstKmerInBackRegionDismissed = false, lastKmerInBacktRegionDismissed = false;
            for(auto v : edits)
            {
                // cout << "loc : " << v << endl;
                if((v >= 0 && v <= minExactMatchSize) || clips.first > 0)
                {
                    firstKmerInFrontRegionDismissed = true;
                    // cout << "firstKmerInFrontRegionDismissed" << endl;
                }
                if((v >= regionSize - minExactMatchSize && v <= regionSize) || (clips.first > regionSize - minExactMatchSize))
                {
                    lastKmerInFrontRegionDismissed = true;
                    // cout << "lastKmerInFrontRegionDismissed" << endl;
                }
                if((v >= contigSize - regionSize && v <= contigSize - regionSize + minExactMatchSize)  || (clips.second > regionSize - minExactMatchSize))
                {
                    firstKmerInBackRegionDismissed = true;
                    // cout << "firstKmerInBackRegionDismissed" << endl;
                }
                if((v >= contigSize - minExactMatchSize && v <= contigSize) || clips.second > 0)
                {
                    lastKmerInBacktRegionDismissed = true;
                    // cout << "lastKmerInBacktRegionDismissed" << endl;
                }
            }
            for(auto v : Indels)
            {
                // cout << "indel loc : " << v << endl;
                if(v >= 0 && v <= minExactMatchSize)
                {
                    firstKmerInFrontRegionDismissed = true;
                    // cout << "firstKmerInFrontRegionDismissed" << endl;
                }
                if(v >= regionSize - minExactMatchSize && v <= regionSize)
                {
                    lastKmerInFrontRegionDismissed = true;
                    // cout << "lastKmerInFrontRegionDismissed" << endl;
                }
                if((v >= contigSize - regionSize && v <= contigSize - regionSize + minExactMatchSize)  || (clips.second > regionSize - minExactMatchSize))
                {
                    firstKmerInBackRegionDismissed = true;
                    // cout << "firstKmerInBackRegionDismissed" << endl;
                }
                if((v >= contigSize - minExactMatchSize && v <= contigSize) || clips.second > 0)
                {
                    lastKmerInBacktRegionDismissed = true;
                    // cout << "lastKmerInBacktRegionDismissed" << endl;
                }
            }
            if(firstKmerInFrontRegionDismissed && lastKmerInFrontRegionDismissed && firstKmerInBackRegionDismissed && lastKmerInBacktRegionDismissed)
            {
                cmp << "the read was supposed to be discarded based on our criteria (edit pos)" << endl;
                return false;
            }
        } else 
        {
            cmp << "the read was supposed to be discarded based on our criteria (clips)" << endl;
            return false;
        }
        return true;
    }

    void comparePmWithBWA(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress, IndexContainer<uint32_t, LevAlign>& pmAlignments, vector<std::pair<uint32_t, uint32_t>>& alnPmBwaHisto, const uint32_t queryCount, string alnPerQueryFileAddress)
    {
        ofstream cmp(comparisonResultsFileAddress);
        ofstream alnPerQ(alnPerQueryFileAddress);
        SamReader bwaSam(cfg.otherToolOutputFileAddress);
        vector<SamReader::Sam> bwaAlignments = bwaSam.parseFile(queryCount);
        uint32_t numberOfEqualalignment = 0, numberOfBWAbetter = 0, numberOfPMbetter = 0, 
        numberOfBWAbetterExceedMaxEdits = 0, numberOfBwaClippedAlignments = 0, 
        numberOfBWAbetterWithLowMatchSize = 0, numberOfBWAbetterNotObserveOurCriteria = 0, 
        onlyBwaFoundMatchForQuery = 0, bwaReadIdFoundInPM = 0, bwaReadIdNotFoundInPM = 0, 
        onlyPmFoundMatchForQuery = 0, nonePmBWAFoundAlignmentForQuery = 0, bwaClipCount = 0, 
        numberOfQueryContainN = 0;
        LevAlign pmAlignment;
        uint32_t pmQueriesFound = 0;
        SamReader::Sam bwaAlignment;
        for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            bool bwaFound = false, pmFound = false;
            uint32_t bwaMatchSize = 0, pmMatchSize = 0;
            string query = queries[queryInd];
            if(query.find('N') != string::npos || query.find('n') != string::npos)
            {
                cmp << "query contains N!" << endl;
                numberOfQueryContainN++; 
                alnPerQ << queryInd << " 0 0"<< endl;
                continue;
            }
            //PM and BWA alignments per query
            alnPerQ << queryInd << " " << pmAlignments.container_.count(queryInd);
            bool bwaFoundReadForQuery = false;
            for(auto it = bwaAlignments.begin(); it != bwaAlignments.end(); it++)
            {
                if(it->queryId == queryInd)
                {
                    bwaFoundReadForQuery = true;
                    break;
                }
            }
            alnPerQ << ((bwaFoundReadForQuery==true) ? " 1" : " 0") << endl;

            for (const SamReader::Sam& aln : bwaAlignments) 
            {
                if(aln.queryId == queryInd)
                {
                    bwaFound = true;
                    bwaAlignment = aln;
                    bwaMatchSize = bwaSam.countMatches(bwaAlignment.cigar);
                    bwaClipCount = bwaSam.countClips(bwaAlignment.cigar);
                    if (bwaClipCount > 0)
                    {
                        numberOfBwaClippedAlignments++;
                    }
                    break;
                }
            }
            auto range = pmAlignments.getRange(queryInd);
            LevAlign pmAlignment;
            for (auto it = range.first; it != range.second; it++) 
            {
                pmFound = true;
                LevAlign aln = it->second;
                if (aln.numberOfMatches > pmAlignment.numberOfMatches)
                {
                    pmAlignment = aln;
                } else if (aln.numberOfMatches == pmAlignment.numberOfMatches)
                {
                    if (aln.numberOfInDel < pmAlignment.numberOfInDel && aln.numberOfSub <= pmAlignment.numberOfSub + 1) // this is a preference on sub over InDel
                    {
                        pmAlignment = aln;
                    } else if (aln.numberOfInDel == pmAlignment.numberOfInDel && aln.numberOfSub < pmAlignment.numberOfSub)
                    {
                        pmAlignment = aln;
                    }
                }
            }
            if (pmFound)
            {
                pmQueriesFound++;
                pmMatchSize = bwaSam.countMatches(pmAlignment.cigar);
            }
            alnPmBwaHisto.push_back(std::make_pair(pmMatchSize, bwaMatchSize));
            if (pmFound && !bwaFound)
            {
                cmp << "PM outperformed and BWA did not found" << endl; 
                numberOfPMbetter++;
                onlyPmFoundMatchForQuery++;
                cmp << "<<<<<<<<<<<PM alignment>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[pmAlignment.readID] << endl;
                // cmp << "q : " << pmAlignment.alignedQuery << endl;
                // cmp << "r : " << pmAlignment.alignedRead << endl;
                // cmp << "E : " << pmAlignment.editDistanceTypes << endl;
                cmp << "CIGAR : " << pmAlignment.cigar << ", subs : " << pmAlignment.numberOfSub << ", query id : " << queryInd << ", read id : " << pmAlignment.readID << ", flag : " << pmAlignment.flag << ", partial match size : " << pmAlignment.partialMatchSize << ", edits : " << pmAlignment.editDistance << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            } else if (!pmFound && bwaFound)
            {
                cmp << "BWA outperformed and PM did not found"; 
                numberOfBWAbetter++;
                onlyBwaFoundMatchForQuery++;
                if (bwaAlignment.editDistance > cfg.editDistance)
                {
                    cmp << "with more edits > " << cfg.editDistance << endl;
                    numberOfBWAbetterExceedMaxEdits++;
                } else if(bwaMatchSize - bwaAlignment.editDistance < cfg.regionSize - cfg.editDistance) {
                    cmp << "with low match size : " << bwaMatchSize << endl;
                    numberOfBWAbetterWithLowMatchSize++;
                } else
                {
                    if(!checkBwaAlignmentBasedOnOurCriteria(cfg.regionSize, cfg.minExactMatchLen, cfg.contigSize, bwaAlignment.cigar, bwaAlignment.mismatchPositions, cmp))
                    {
                        numberOfBWAbetterNotObserveOurCriteria++;
                    }
                }
                cmp << "<<<<<<<<<<<<<BWA alignment>>>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[bwaAlignment.readId] << endl;
                cmp << "CIGAR : " << bwaAlignment.cigar << ", subs : " << bwaAlignment.editDistance << ", query id : " << bwaAlignment.queryId << ", read id : " << bwaAlignment.readId << ", flag : " << bwaAlignment.flag << ", start pos : " << bwaAlignment.pos << ", mismatchPositions : " << bwaAlignment.mismatchPositions << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";  
            } else if (pmFound && bwaFound) 
            {
                // int bwaInDelSize = bwaSam.countInsertions(bwaAlignment.cigar) + bwaSam.countDeletions(bwaAlignment.cigar);
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
               
                if (bwaMatchSize < pmMatchSize)
                {
                    cmp << "PM outperformed "; 
                    numberOfPMbetter++;
                    cmp << "and # of BWA clips are: " << bwaClipCount << endl;
                } else if (bwaMatchSize > pmMatchSize)
                {
                    cmp << "BWA outperformed";
                    numberOfBWAbetter++;
                    if (bwaAlignment.editDistance > cfg.editDistance)
                    {
                        cmp << "with more edits > " << cfg.editDistance << endl;
                        numberOfBWAbetterExceedMaxEdits++;
                    } else if(bwaMatchSize - bwaAlignment.editDistance < cfg.regionSize - cfg.editDistance) {
                        cmp << "with low match size : " << bwaMatchSize << endl;
                        numberOfBWAbetterWithLowMatchSize++;
                    } else
                    {
                        if(!checkBwaAlignmentBasedOnOurCriteria(cfg.regionSize, cfg.minExactMatchLen, cfg.contigSize, bwaAlignment.cigar, bwaAlignment.mismatchPositions, cmp))
                        {
                            numberOfBWAbetterNotObserveOurCriteria++;
                        }
                    }
                } else
                {
                    if (pmAlignment.numberOfSub < bwaAlignment.editDistance)
                    {
                        cmp << "PM outperformed with less subs" << endl;
                        numberOfPMbetter++;
                    } else if (pmAlignment.numberOfSub > bwaAlignment.editDistance)
                    {
                        cmp << "BWA outperformed with less subs" << endl;
                        numberOfBWAbetter++;
                        if (bwaAlignment.editDistance > cfg.editDistance)
                        {
                            numberOfBWAbetterExceedMaxEdits++;
                        }  
                    } else
                    {
                        cmp << "BWA and PM performed equal" << endl;   
                        numberOfEqualalignment++;
                    }
                }
                //check for the read id
                auto range = pmAlignments.getRange(bwaAlignment.queryId);
                bool found = false;
                for (auto itt=range.first; itt != range.second; itt++)
                {
                    LevAlign aln = itt->second;
                    if((uint32_t) aln.readID == bwaAlignment.readId)
                    {
                        bwaReadIdFoundInPM++;
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    bwaReadIdNotFoundInPM++;
                    cmp << "bwaReadIdNotFoundInPM = Q :" << bwaAlignment.queryId << ", R : " << bwaAlignment.readId << endl;
                }
                cmp << "<<<<<<<<<<<PM alignment>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[pmAlignment.readID] << endl;
                // cmp << "q : " << pmAlignment.alignedQuery << endl;
                // cmp << "r : " << pmAlignment.alignedRead << endl;
                // cmp << "E : " << pmAlignment.editDistanceTypes << endl;
                cmp << "CIGAR : " << pmAlignment.cigar << ", subs : " << pmAlignment.numberOfSub << ", query id : " << queryInd << ", read id : " << pmAlignment.readID << ", flag : " << pmAlignment.flag << ", partial match size : " << pmAlignment.partialMatchSize << ", edits : " << pmAlignment.editDistance << endl;
                cmp << "<<<<<<<<<<<<<BWA alignment>>>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[bwaAlignment.readId] << endl;
                cmp << "CIGAR : " << bwaAlignment.cigar << ", subs : " << bwaAlignment.editDistance << ", query id : " << bwaAlignment.queryId << ", read id : " << bwaAlignment.readId << ", flag : " << bwaAlignment.flag << ", start pos : " << bwaAlignment.pos << ", mismatchPositions : " << bwaAlignment.mismatchPositions << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            } else {
                cmp << "query " << queryInd << " not found in the PM results nor in BWA" << endl;
                nonePmBWAFoundAlignmentForQuery++;
            }
        }
        cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Overall Comparison Results>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cmp << left << setw(80) << "# Of queries : " << queryCount << endl;
        cmp << left << setw(80) << "# Of queries that PARMIK found match : " << pmQueriesFound << endl;
        cmp << left << setw(80) << "# Of queries that BWA found match : " << bwaAlignments.size() << endl;
        cmp << left << setw(80) << "# Of Matched Hits between PM and BWA : " << numberOfEqualalignment << endl;
        cmp << left << setw(80) << "# Of PM outperformed : " << numberOfPMbetter << endl;
        cmp << left << setw(80) << "# Of BWA outperformed : " << numberOfBWAbetter << endl;
        cmp << left << setw(80) << "# Of BWA outperformed that Exceed Max Edits : " << numberOfBWAbetterExceedMaxEdits << endl;
        cmp << left << setw(80) << "# Of BWA outperformed that matchsize < (R - E): " << numberOfBWAbetterWithLowMatchSize << endl;
        cmp << left << setw(80) << "# Of BWA outperformed that has edits or clips in minExact-mers of region : " << numberOfBWAbetterNotObserveOurCriteria << endl;
        cmp << left << setw(80) << "# Of BWA clipped alignments : " << numberOfBwaClippedAlignments << endl;
        cmp << left << setw(80) << "# Of queries that BWA found matches, not PM : " << onlyBwaFoundMatchForQuery << endl;
        cmp << left << setw(80) << "# Of queries that PM found matches, not BWA : " << onlyPmFoundMatchForQuery << endl;
        cmp << left << setw(80) << "# Of readIDs that BWA found but PM did not : " << bwaReadIdNotFoundInPM << endl;
        cmp << left << setw(80) << "# Of queries that neigther BWA nor PM found any match : " << nonePmBWAFoundAlignmentForQuery << endl;
        cmp << left << setw(80) << "# Of queries contains N: " << numberOfQueryContainN << endl;
        cmp.close();
        alnPerQ.close();
    }
};

#endif