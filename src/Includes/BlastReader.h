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
#include "Alignment.h"

using namespace std;

class BlastReader {
public:

    BlastReader(const string& filename) : filename_(filename) {}

    uint32_t getIndels(string read, string query) {
        uint32_t inDelCount = 0;
        if (read == "" || query == "") {
            // cout << "q : " << query << endl;
            // cout << "r : " << read << endl;
            cerr << "Error: Sequences are empty." << endl;
            return 0;
        }
        // Check if the sequences have the same length
        if (read.length() != query.length()) {
            // cout << "q : " << query << endl;
            // cout << "r : " << read << endl;
            cerr << "Error: Sequences have different lengths." << endl;
            return 0;
        }
        // Iterate through the sequences and compare corresponding characters
        for (size_t i = 0; i < read.length(); ++i) {
            if (read[i] != query[i]) {
                if (read[i] == '-' || query[i] == '-') {
                    inDelCount++;
                }
            } 
        }
 
        return inDelCount;
    }

    void parseFile(uint32_t lastQueryID, IndexContainer<uint32_t, Alignment>& alignments) {
        
        ifstream blastFile(filename_);
        Utilities<uint32_t> utl;
        if (!blastFile) {
            cerr << "Error opening file: " << filename_ << endl;
            return;
        }

        string line;
        uint32_t queryCount = 0;
        while (getline(blastFile, line)) {
            if (line.empty() || line == "" || line.find('#') != string::npos) {
                continue; // Skip empty lines
            }
            Alignment blast;
            istringstream iss(line);
            string queryID, readID, queryS, queryE, readS, readE, AlignmentLength, Mismatches, flag, queryAligned, readAligned;
            iss >> queryID >>  readID >>  queryS >>  queryE >>  readS >>  readE >>  AlignmentLength >>  Mismatches >>  flag >> queryAligned >> readAligned;

            if (queryID == "") blast.queryID = 0; else blast.queryID = utl.extractContigId(queryID);
            if (readID == "") blast.readID = 0; else blast.readID = utl.extractContigId(readID);
            if (queryS == "") blast.queryRegionStartPos = 0; else blast.queryRegionStartPos = stoi(queryS.c_str());
            if (queryE == "") blast.queryRegionEndPos = 0; else blast.queryRegionEndPos = stoi(queryE.c_str());
            if (readS == "") blast.readRegionStartPos = 0; else blast.readRegionStartPos = stoi(readS.c_str());
            if (readE == "") blast.readRegionEndPos = 0; else blast.readRegionEndPos = stoi(readE.c_str());
            if (AlignmentLength == "") blast.partialMatchSize = 0; else blast.partialMatchSize = stoi(AlignmentLength.c_str());
            if (Mismatches == "") blast.substitutions = 0; else blast.substitutions = stoi(Mismatches.c_str());
            if (flag == "") blast.flag = 0; else if (flag == "plus") blast.flag = 0; else blast.flag = 16;
            if (queryAligned == "") blast.alignedQuery = ""; else blast.alignedQuery = queryAligned;
            if (readAligned == "") blast.alignedRead = ""; else blast.alignedRead = readAligned;

            blast.inDels = getIndels(readAligned, queryAligned);

            if(blast.queryID > (int)lastQueryID)
                break;

            alignments.put(blast.queryID, blast);
            queryCount++;
        }
    }

private:
    string filename_;
};

#endif 