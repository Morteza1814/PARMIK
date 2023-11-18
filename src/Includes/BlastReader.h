#ifndef BLASTREADER_H
#define BLASTREADER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include "Utils.h"

using namespace std;

class BlastReader {
public:
    struct Blast {
        uint32_t queryId = 0;
        uint32_t readId = 0;
        uint32_t queryS = 0;
        uint32_t queryE = 0;
        uint32_t readS = 0;
        uint32_t readE = 0;
        uint32_t AlignmentLength = 0;
        uint32_t Mismatches = 0;
        uint32_t flag = 0;
        uint32_t InDels = 0;
    };

    BlastReader(const string& filename) : filename_(filename) {}

    void parseFile(uint32_t lastQueryID, IndexContainer<uint32_t, Blast>& alignments) {
        
        ifstream blastFile(filename_);
        Utilities<uint32_t> utl;
        if (!blastFile) {
            cerr << "Error opening file: " << filename_ << endl;
            return;
        }

        string line;
        uint32_t queryCount = 0;
        while (getline(blastFile, line)) {
            if (line.empty() || line == "") {
                continue; // Skip empty lines
            }
            Blast blast;
            istringstream iss(line);
            string queryID, readID, queryS, queryE, readS, readE, AlignmentLength, Mismatches, flag;
            iss >> queryID >>  readID >>  queryS >>  queryE >>  readS >>  readE >>  AlignmentLength >>  Mismatches >>  flag;

            if (queryID == "") blast.queryId = 0; else blast.queryId = utl.extractContigId(queryID);
            if (readID == "") blast.readId = 0; else blast.readId = utl.extractContigId(readID);
            if (queryS == "") blast.queryS = 0; else blast.queryS = atoi(queryS.c_str());
            if (queryE == "") blast.queryE = 0; else blast.queryE = atoi(queryE.c_str());
            if (readS == "") blast.readS = 0; else blast.readS = atoi(readS.c_str());
            if (readE == "") blast.readE = 0; else blast.readE = atoi(readE.c_str());
            if (AlignmentLength == "") blast.AlignmentLength = 0; else blast.AlignmentLength = atoi(AlignmentLength.c_str());
            if (Mismatches == "") blast.Mismatches = 0; else blast.Mismatches = atoi(Mismatches.c_str());
            if (flag == "") blast.flag = 0; else if (flag == "plus") blast.flag = 0; else blast.flag = 16;

            if(blast.queryId > lastQueryID)
                break;

            alignments.put(blast.queryId, blast);
            queryCount++;
        }
    }

    bool isReverseComplemented(const Blast& blast) const {
        return (blast.flag & 0x10) != 0; // Check if the 5th bit (0x10) is set
    }

private:
    string filename_;
};

#endif 