#ifndef EVALUATE_SECOND_CHANCE_H
#define EVALUATE_SECOND_CHANCE_H

#include <iostream>
#include <vector>
#include "SamReader.h"
#include "Alignment.h"
#include "Utils.h"
#include "tsl/robin_map.h"
#include "omp.h"

using namespace std;

class EvaluateSecondChance {
public:
    void evaluateSecondChance(string parmikAlnFileAddressWithSC, string parmikAlnFileAddressWithoutSC, const uint32_t queryCount, string outputDir, string queryFileAddress){
        cout << "<<<<<<<<<<<<<<<<<<<<Evaluate Second Chance>>>>>>>>>>>>>>>>>>>>" << endl;
        Utilities<uint32_t> util;
        vector<Alignment> parmikAlignmentsWithSC;
        SamReader parmikSamSC(parmikAlnFileAddressWithSC);
        parmikSamSC.parseFile(queryCount, parmikAlignmentsWithSC, false);
        cout << "Number of alignments with second chance: " << parmikAlignmentsWithSC.size() << endl;
        vector<Alignment> parmikAlignmentsWithoutSC;
        SamReader parmikSamNSC(parmikAlnFileAddressWithoutSC);
        parmikSamNSC.parseFile(queryCount, parmikAlignmentsWithoutSC, false);
        cout << "Number of alignments without second chance: " << parmikAlignmentsWithoutSC.size() << endl;
        tsl::robin_map <uint32_t, string> reads, queries;
        uint32_t baseAddress = 0;
        uint32_t qc = util.readContigsFromFile(queryFileAddress, queryCount, queries, baseAddress);
        ofstream out(outputDir + "/SC_vs_NoSC.txt");
        ofstream scmissedaln(outputDir + "/SC_MissedAlignments.txt");
        tsl::robin_map<int32_t, uint32_t> scAlignmentSizeImprovements;
        #pragma omp parallel for schedule(dynamic, 1)
        for (uint32_t queryInd = 0; queryInd < qc; queryInd++){
            if (queryInd % 1000 == 0) {
                cout << "Processing Query " << queryInd << endl;
            }
            tsl::robin_map<uint32_t, Alignment> parmikAlignmentsWithSCMap, parmikAlignmentsWithoutSCMap;
            for (Alignment it : parmikAlignmentsWithSC) {
                if (it.queryID == (int)queryInd){
                    if (parmikAlignmentsWithSCMap.find(it.readID) == parmikAlignmentsWithSCMap.end()){
                        parmikAlignmentsWithSCMap[it.readID] = it;
                    } else {
                        if (it.partialMatchSize > parmikAlignmentsWithSCMap[it.readID].partialMatchSize){
                            parmikAlignmentsWithSCMap[it.readID] = it;
                        }
                    }
                }
            }

            for (Alignment it : parmikAlignmentsWithoutSC) {
                if (it.queryID == (int)queryInd){
                    if (parmikAlignmentsWithoutSCMap.find(it.readID) == parmikAlignmentsWithoutSCMap.end()){
                        parmikAlignmentsWithoutSCMap[it.readID] = it;
                    } else {
                        if (it.partialMatchSize > parmikAlignmentsWithoutSCMap[it.readID].partialMatchSize){
                            parmikAlignmentsWithoutSCMap[it.readID] = it;
                        }
                    }
                }
            }
            std::map<int32_t, int> localScAlignmentSizeImprovements;
            for (auto it : parmikAlignmentsWithoutSCMap){
                auto itt = parmikAlignmentsWithSCMap.find(it.first);
                if (itt != parmikAlignmentsWithSCMap.end()){
                    #pragma omp critical
                    {
                        out << it.second.partialMatchSize << " " << itt->second.partialMatchSize << endl;
                    }
                    int32_t scAlignmentSizeImprovement = itt->second.partialMatchSize - it.second.partialMatchSize;
                    #pragma omp critical
                    {
                        if (scAlignmentSizeImprovement < 0){
                            scmissedaln << "SCQID: " << it.second.queryID << ", SCRID: " << it.second.readID 
                            << ", noSCQID: " << itt->second.queryID << ", noSCRID: " << itt->second.readID 
                            << ", SC partialMatchSize: " << itt->second.partialMatchSize << ", noSC partialMatchSize: " << it.second.partialMatchSize
                            << ", SC cigar: " << itt->second.cigar << ", noSC cigar: " << it.second.cigar  << ", scAlignmentSizeImprovement: " << scAlignmentSizeImprovement << endl;
                        }
                    }
                    if (scAlignmentSizeImprovements.find (scAlignmentSizeImprovement) != scAlignmentSizeImprovements.end()){
                        localScAlignmentSizeImprovements[scAlignmentSizeImprovement]++;
                    } else{
                        localScAlignmentSizeImprovements[scAlignmentSizeImprovement] = 1;
                    }
                }
            }
            #pragma omp critical
            {
                for (auto &entry : localScAlignmentSizeImprovements) {
                    scAlignmentSizeImprovements[entry.first] += entry.second;
                }
            }
        }
        cout << "---------------------scAlignmentSizeImprovements--------------------------" << endl;
        for (auto it : scAlignmentSizeImprovements){
            cout << it.first << " " << it.second << endl;
        }
        out.close();
    }
};

#endif