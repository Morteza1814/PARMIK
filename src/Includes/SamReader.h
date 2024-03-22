#ifndef SAMREADER_H
#define SAMREADER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include "Utils.h"

using namespace std;

class SamReader {
public:
    struct Sam {
        uint32_t queryId;
        int flag;
        uint32_t readId;
        int pos;
        int mapQ;
        string cigar;
        string rNext;
        string pNext;
        int tLen;
        string seq;
        string qual;
        unsigned int editDistance; // NM flag
        string mismatchPositions; // MD flag
        int alignmentScore; // AS flag
    };

    struct EditCounts {
        uint32_t matchCount = 0;
        uint32_t subCount = 0;
        uint32_t insertCount = 0;
        uint32_t deleteCount = 0;
        uint32_t clipCount = 0;
    };

    SamReader(const string& filename) : filename_(filename) {}

    EditCounts extractEditTypes(const string& cigar) const {
        EditCounts counts;
        regex cigarRegex("([0-9]+)([MIDNSHPX=])");
        sregex_iterator iter(cigar.begin(), cigar.end(), cigarRegex);
        sregex_iterator end;

        for (; iter != end; ++iter) {
            char op = (*iter)[2].str()[0];
            cout << "jaaaaaaaa" << endl;
            int length = stoi((*iter)[1].str());
            switch (op) {
                case 'X':
                    counts.subCount += length;
                    break;
                case 'M':
                case '=':
                    counts.matchCount += length;
                    break;
                case 'I':
                    counts.insertCount += length;
                    break;
                case 'D':
                    counts.deleteCount += length;
                    break;
                case 'S':
                case 'H':
                    counts.clipCount += length;
                    break;
                default:
                    break;
            }
        }

        return counts;
    }

    uint32_t countMatches(const string& cigar) const {
        return extractEditTypes(cigar).matchCount;
    }

    uint32_t countSubstitutions(const string& cigar) const {
        return extractEditTypes(cigar).subCount;
    }

    uint32_t countInsertions(const string& cigar) const {
        return extractEditTypes(cigar).insertCount;
    }

    uint32_t countDeletions(const string& cigar) const {
        return extractEditTypes(cigar).deleteCount;
    }

    uint32_t countClips(const string& cigar) const {
        return extractEditTypes(cigar).clipCount;
    }

    vector<Sam> parseFile(uint32_t lastQueryID) {
        vector<Sam> alignments;
        ifstream samFile(filename_);
        Utilities<uint32_t> utl;
        if (!samFile) {
            cerr << "Error opening file: " << filename_ << endl;
            return alignments;
        }

        string line;
        uint32_t queryCount = 0;
        while (getline(samFile, line)) {
            if (line.empty() || line[0] == '@') {
                continue; // Skip header lines
            }
            Sam sam;
            istringstream iss(line);
            string queryID, readID, flag, pos, mapQ, cigar, rNext, pNext, tLen, seq, qual;
            iss >> queryID >>  flag >> readID >>  pos >>  mapQ >> cigar >>   rNext >>  pNext >>  tLen >>  seq >>  qual;
            
            size_t found = readID.find("*");
            if (found != string::npos)
            {
                queryCount++;
                continue;
            }
            cout << "daaaaaaaaar" << endl;
            cout << "flag: " << flag << " pos:" << pos << " map: " <<  mapQ << "tLen:" << tLen << endl;
            if (flag.find("*") != string::npos) sam.flag = 0; else sam.flag = stoi(flag.c_str());
            if (pos.find("*") != string::npos) sam.pos = 0; else sam.pos = stoi(pos.c_str());
            if (mapQ.find("*") != string::npos) sam.mapQ = 0; else sam.mapQ = stoi(mapQ.c_str());
            if (cigar.find("*") != string::npos) sam.cigar = ""; else sam.cigar = cigar;
            if (rNext.find("*") != string::npos) sam.rNext = ""; else sam.rNext = rNext;
            if (pNext.find("*") != string::npos) sam.pNext = ""; else sam.pNext = pNext;
            if (tLen.find("*") != string::npos) sam.tLen = 0; else sam.tLen = stoi(tLen.c_str());
            if (seq.find("*") != string::npos) sam.seq = ""; else sam.seq = seq;
            if (qual.find("*") != string::npos) sam.qual = ""; else sam.qual = qual;

            sam.queryId = utl.extractContigId(queryID);
            sam.readId = utl.extractContigId(readID);
            if(sam.queryId > lastQueryID)
                break;
            // Parse other flags as needed

            // Find the CIGAR field
            // for (int i = 0; i < 9; ++i) {
            //     iss >> sam.cigar;
            // }

            string token;
            while (iss >> token) {
                if (token.find("NM:") == 0) {
                    sam.editDistance = stoi(token.substr(5));
                } else if (token.find("MD:") == 0) {
                    sam.mismatchPositions = token.substr(5);
                } else if (token.find("AS:") == 0) {
                    sam.alignmentScore = stoi(token.substr(5));
                }
            }
            alignments.push_back(sam);
            queryCount++;
        }

        return alignments;
    }

    bool isReverseComplemented(const Sam& sam) const {
        return (sam.flag & 0x10) != 0; // Check if the 5th bit (0x10) is set
    }

private:
    string filename_;
};

#endif 