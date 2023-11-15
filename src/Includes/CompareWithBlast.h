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

    void comparePmWithBlast(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, string comparisonResultsFileAddress, IndexContainer<uint32_t, LevAlign>& pmAlignments, vector<std::pair<uint32_t, uint32_t>>& alnPmBLASTHisto, const uint32_t queryCount, string alnPerQueryFileAddress)
    {
        ofstream cmp(comparisonResultsFileAddress);
        ofstream alnPerQ(alnPerQueryFileAddress);
        SamReader sam(cfg.otherToolOutputFileAddress);
        BlastReader blastReader(cfg.otherToolOutputFileAddress);
        IndexContainer<uint32_t, BlastReader::Blast> blastAlignments;
        blastReader.parseFile(queryCount, blastAlignments);
        uint32_t numberOfEqualalignment = 0, numberOfBLASTbetter = 0, numberOfPMbetter = 0, 
        numberOfBLASTbetterExceedMaxEdits = 0, numberOfBLASTbetterWithLowMatchSize = 0,  
        onlyblastFoundMatchForQuery = 0, pmReadIdNotFoundInBLAST = 0, blastReadIdNotFoundInPM = 0, 
        onlyPmFoundMatchForQuery = 0, nonePmblastFoundAlignmentForQuery = 0,  
        numberOfQueryContainN = 0;
        LevAlign pmAlignment;
        uint32_t pmQueriesFound = 0;
        BlastReader::Blast blastAlignment;
        for(uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            bool blastFound = false, pmFound = false;
            uint32_t blastMatchSize = 0, pmMatchSize = 0;
            string query = queries[queryInd];
            if(query.find('N') != string::npos || query.find('n') != string::npos)
            {
                cmp << "query contains N!" << endl;
                numberOfQueryContainN++; 
                alnPerQ << queryInd << " 0 0"<< endl;
                continue;
            }
            //PM and BLAST alignments per query
            alnPerQ << queryInd << " " << pmAlignments.container_.count(queryInd) << " " << blastAlignments.container_.count(queryInd);

            BlastReader::Blast blastAlignment;
            auto blRrange = blastAlignments.getRange(queryInd);
            for (auto it = blRrange.first; it != blRrange.second; it++) 
            {
                blastFound = true;
                BlastReader::Blast aln = it->second;
                if (aln.AlignmentLength > blastAlignment.AlignmentLength)
                {
                    blastAlignment = aln;
                } else if (aln.AlignmentLength == blastAlignment.AlignmentLength)
                {
                    if (aln.AlignmentLength == blastAlignment.AlignmentLength && aln.Mismatches < blastAlignment.Mismatches)
                    {
                        blastAlignment = aln;
                    }
                }
            }
            auto pmRange = pmAlignments.getRange(queryInd);
            LevAlign pmAlignment;
            for (auto it = pmRange.first; it != pmRange.second; it++) 
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
                pmMatchSize = sam.countMatches(pmAlignment.cigar);
            }
            if (blastFound)
            {
                blastMatchSize = blastAlignment.AlignmentLength;
            }
            alnPmBLASTHisto.push_back(std::make_pair(pmMatchSize, blastMatchSize));
            if (pmFound && !blastFound)
            {
                cmp << "PM outperformed and BLAST did not found" << endl; 
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
            } else if (!pmFound && blastFound)
            {
                cmp << "BLAST outperformed and PM did not found"; 
                numberOfBLASTbetter++;
                onlyblastFoundMatchForQuery++;
                if (blastAlignment.Mismatches > cfg.editDistance)
                {
                    cmp << "with more edits > " << cfg.editDistance << endl;
                    numberOfBLASTbetterExceedMaxEdits++;
                } else if(blastMatchSize - blastAlignment.Mismatches < cfg.regionSize - cfg.editDistance) {
                    cmp << "with low match size : " << blastMatchSize << endl;
                    numberOfBLASTbetterWithLowMatchSize++;
                }
                cmp << "<<<<<<<<<<<<<BLAST alignment>>>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[blastAlignment.readId] << endl;
                cmp << "AlignmentLength : " << blastAlignment.AlignmentLength << ", subs : " << blastAlignment.Mismatches << ", query id : " << blastAlignment.queryId << ", read id : " << blastAlignment.readId << ", flag : " << blastAlignment.flag << ", query start pos : " << blastAlignment.queryS << ", read start pos : " << blastAlignment.readS << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";  
            } else if (pmFound && blastFound) 
            {
                // int bwaInDelSize = bwaSam.countInsertions(blastAlignment.cigar) + bwaSam.countDeletions(blastAlignment.cigar);
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
               
                if (blastMatchSize < pmMatchSize)
                {
                    cmp << "PM outperformed "; 
                    numberOfPMbetter++;
                } else if (blastMatchSize > pmMatchSize)
                {
                    cmp << "BLAST outperformed";
                    numberOfBLASTbetter++;
                    if (blastAlignment.Mismatches > cfg.editDistance)
                    {
                        cmp << "with more edits > " << cfg.editDistance << endl;
                        numberOfBLASTbetterExceedMaxEdits++;
                    } else if(blastMatchSize - blastAlignment.Mismatches < cfg.regionSize - cfg.editDistance) {
                        cmp << "with low match size : " << blastMatchSize << endl;
                        numberOfBLASTbetterWithLowMatchSize++;
                    }
                } else
                {
                    if (pmAlignment.numberOfSub < blastAlignment.Mismatches)
                    {
                        cmp << "PM outperformed with less subs" << endl;
                        numberOfPMbetter++;
                    } else if (pmAlignment.numberOfSub > blastAlignment.Mismatches)
                    {
                        cmp << "BLAST outperformed with less subs" << endl;
                        numberOfBLASTbetter++;
                        if (blastAlignment.Mismatches > cfg.editDistance)
                        {
                            numberOfBLASTbetterExceedMaxEdits++;
                        }  
                    } else
                    {
                        cmp << "BLAST and PM performed equal" << endl;   
                        numberOfEqualalignment++;
                    }
                }
                //check for the read id
                bool found = false;
                for (auto it = pmRange.first; it != pmRange.second; it++) 
                {
                    LevAlign aln = it->second;
                    for (auto itt = blRrange.first; itt != blRrange.second; itt++) 
                    {
                        BlastReader::Blast blaln = itt->second;
                        if((uint32_t) aln.readID == blaln.readId)
                        {
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        pmReadIdNotFoundInBLAST++;
                        // cmp << "pmReadIdNotFoundInBLAST = Q :" << blastAlignment.queryId << ", R : " << blastAlignment.readId << endl;
                    }
                }

                for (auto itt = blRrange.first; itt != blRrange.second; itt++) 
                {
                    BlastReader::Blast blaln = itt->second;
                    for (auto it = pmRange.first; it != pmRange.second; it++) 
                    {
                        LevAlign aln = it->second;
                        if((uint32_t) aln.readID == blaln.readId)
                        {
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        blastReadIdNotFoundInPM++;
                        // cmp << "pmReadIdNotFoundInBLAST = Q :" << blastAlignment.queryId << ", R : " << blastAlignment.readId << endl;
                    }
                }
                cmp << "<<<<<<<<<<<PM alignment>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[pmAlignment.readID] << endl;
                // cmp << "q : " << pmAlignment.alignedQuery << endl;
                // cmp << "r : " << pmAlignment.alignedRead << endl;
                // cmp << "E : " << pmAlignment.editDistanceTypes << endl;
                cmp << "CIGAR : " << pmAlignment.cigar << ", subs : " << pmAlignment.numberOfSub << ", query id : " << queryInd << ", read id : " << pmAlignment.readID << ", flag : " << pmAlignment.flag << ", partial match size : " << pmAlignment.partialMatchSize << ", edits : " << pmAlignment.editDistance << endl;
                cmp << "<<<<<<<<<<<<<BLAST alignment>>>>>>>>>>>>>" << endl;
                cmp << "Q : " << query << endl;
                cmp << "R : " << reads[blastAlignment.readId] << endl;
                cmp << "AlignmentLength : " << blastAlignment.AlignmentLength << ", subs : " << blastAlignment.Mismatches << ", query id : " << blastAlignment.queryId << ", read id : " << blastAlignment.readId << ", flag : " << blastAlignment.flag << ", query start pos : " << blastAlignment.queryS << ", read start pos : " << blastAlignment.readS << endl;
                cmp << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";  
            } else {
                cmp << "query " << queryInd << " not found in the PM results nor in BLAST" << endl;
                nonePmblastFoundAlignmentForQuery++;
            }
        }
        cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Overall Comparison Results>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cmp << left << setw(80) << "# Of queries : " << queryCount << endl;
        cmp << left << setw(80) << "# Of queries that PARMIK found match : " << pmQueriesFound << endl;
        cmp << left << setw(80) << "# Of queries that BLAST found match : " << blastAlignments.size() << endl;
        cmp << left << setw(80) << "# Of Matched Hits between PM and BLAST : " << numberOfEqualalignment << endl;
        cmp << left << setw(80) << "# Of PM outperformed : " << numberOfPMbetter << endl;
        cmp << left << setw(80) << "# Of BLAST outperformed : " << numberOfBLASTbetter << endl;
        cmp << left << setw(80) << "# Of BLAST outperformed that Exceed Max Edits : " << numberOfBLASTbetterExceedMaxEdits << endl;
        cmp << left << setw(80) << "# Of BLAST outperformed that matchsize < (R - E): " << numberOfBLASTbetterWithLowMatchSize << endl;
        cmp << left << setw(80) << "# Of queries that BLAST found matches, not PM : " << onlyblastFoundMatchForQuery << endl;
        cmp << left << setw(80) << "# Of queries that PM found matches, not BLAST : " << onlyPmFoundMatchForQuery << endl;
        cmp << left << setw(80) << "# Of readIDs that PM found but BLAST did not : " << pmReadIdNotFoundInBLAST << endl;
        cmp << left << setw(80) << "# Of readIDs that BLAST found but PM did not : " << blastReadIdNotFoundInPM << endl;
        cmp << left << setw(80) << "# Of queries that neigther BLAST nor PM found any match : " << nonePmblastFoundAlignmentForQuery << endl;
        cmp << left << setw(80) << "# Of queries contains N: " << numberOfQueryContainN << endl;
        cmp.close();
        alnPerQ.close();
    }
};

#endif