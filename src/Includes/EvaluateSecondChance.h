#ifndef EVALUATE_SECOND_CHANCE_H
#define EVALUATE_SECOND_CHANCE_H

#include <iostream>
#include <vector>
#include "SamReader.h"
#include "Alignment.h"
#include "Utils.h"
#include "tsl/robin_map.h"

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
        ofstream out(outputDir + "SC_vs_NoSC.txt");
        tsl::robin_map<int32_t, uint32_t> scAlignmentSizeImprovents;
        for (uint32_t queryInd = 0; queryInd < qc; queryInd++){
            tsl::robin_map<uint32_t, Alignment> parmikAlignmentsWithSCMap, parmikAlignmentsWithoutSCMap;
            for (Alignment it : parmikAlignmentsWithSC) {
                if (it.queryID == (int)queryInd){
                    parmikAlignmentsWithSCMap[it.readID] = it;
                }
            }

            for (Alignment it : parmikAlignmentsWithoutSC) {
                if (it.queryID == (int)queryInd){
                    parmikAlignmentsWithoutSCMap[it.readID] = it;
                }
            }

            for (auto it : parmikAlignmentsWithoutSCMap){
                auto itt = parmikAlignmentsWithSCMap.find(it.first);
                if (itt != parmikAlignmentsWithSCMap.end()){
                    out << it.second.partialMatchSize << " " << itt->second.partialMatchSize << endl;
                    int32_t scAlignmentSizeImprovent = itt->second.partialMatchSize - it.second.partialMatchSize;
                    if (scAlignmentSizeImprovents.find (scAlignmentSizeImprovent) != scAlignmentSizeImprovents.end()){
                        scAlignmentSizeImprovents[scAlignmentSizeImprovent]++;
                    } else{
                        scAlignmentSizeImprovents[scAlignmentSizeImprovent] = 1;
                    }
                }
            }
        }
        cout << "---------------------scAlignmentSizeImprovents--------------------------" << endl;
        for (auto it : scAlignmentSizeImprovents){
            cout << it.first << " " << it.second << endl;
        }
        out.close();
    }
};

#endif