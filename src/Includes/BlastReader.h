#ifndef BLASTREADER_H
#define BLASTREADER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include <limits>
#include <stdexcept>
#include "Utils.h"
#include "Alignment.h"

using namespace std;

class BlastReader {
    static uint32_t number(const string& value, uint32_t maximum) {
        if (value.empty() || value.find_first_not_of("0123456789") != string::npos)
            throw runtime_error("expected an unsigned integer: " + value);
        const auto result = stoull(value);
        if (result > maximum)
            throw runtime_error("integer outside supported range: " + value);
        return static_cast<uint32_t>(result);
    }

    static uint32_t contigId(const string& name) {
        const auto dot = name.find_last_of('.');
        return number(name.substr(dot == string::npos ? 0 : dot + 1),
                      numeric_limits<int>::max());
    }

    static uint32_t coordinate(const string& value) {
        // The input is one-based; Alignment stores zero-based 8-bit positions.
        const auto position = number(value, numeric_limits<uint8_t>::max() + 1u);
        if (position == 0) throw runtime_error("alignment coordinates must be positive");
        return position - 1;
    }

public:

    BlastReader(const string& filename) : filename_(filename) {}

    uint32_t getIndels(const string& read, const string& query) const {
        uint32_t inDelCount = 0;
        if (read.empty() || query.empty())
            throw runtime_error("aligned sequences must not be empty");
        if (read.length() != query.length())
            throw runtime_error("aligned sequences have different lengths");
        // Iterate through the sequences and compare corresponding characters
        for (size_t i = 0; i < read.length(); ++i) {
            if (read[i] == '-' && query[i] == '-')
                throw runtime_error("alignment column contains two gaps");
            if (read[i] == '-' || query[i] == '-') ++inDelCount;
        }
 
        return inDelCount;
    }

    void parseFile(uint32_t lastQueryID, vector<Alignment>& alignments) {
        ifstream blastFile(filename_);
        if (!blastFile) throw runtime_error("Cannot open BLAST alignments: " + filename_);

        string line;
        size_t lineNumber = 0;
        while (getline(blastFile, line)) {
            ++lineNumber;
            const auto first = line.find_first_not_of(" \t\r");
            if (first == string::npos || line[first] == '#') continue;
            try {
                istringstream iss(line);
                string queryID, readID, queryS, queryE, readS, readE, alignmentLength,
                       mismatches, strand, queryAligned, readAligned;
                if (!(iss >> queryID >> readID >> queryS >> queryE >> readS >> readE
                          >> alignmentLength >> mismatches >> strand >> queryAligned >> readAligned))
                    throw runtime_error("expected 11 BLAST fields including aligned sequences");
                const auto queryId = contigId(queryID);
                // BLAST records need not be sorted by query ID.
                if (queryId > lastQueryID) continue;
                const auto length = number(alignmentLength, numeric_limits<uint8_t>::max());
                const auto substitutions = number(mismatches, length);
                if (length == 0 || queryAligned.size() != length || readAligned.size() != length)
                    throw runtime_error("alignment length disagrees with aligned sequences");
                const auto indels = getIndels(readAligned, queryAligned);
                if (substitutions + indels > length)
                    throw runtime_error("mismatches and gaps exceed alignment length");
                if (strand != "plus" && strand != "minus")
                    throw runtime_error("expected plus or minus BLAST strand");

                Alignment blast;
                blast.queryID = queryId;
                blast.readID = contigId(readID);
                // Preserve endpoint order for reverse-strand alignments.
                blast.queryRegionStartPos = coordinate(queryS);
                blast.queryRegionEndPos = coordinate(queryE);
                blast.readRegionStartPos = coordinate(readS);
                blast.readRegionEndPos = coordinate(readE);
                blast.partialMatchSize = length;
                blast.substitutions = substitutions;
                blast.inDels = indels;
                blast.matches = length - substitutions - indels;
                blast.editDistance = substitutions + indels;
                blast.flag = strand == "plus" ? 0 : 16;
                blast.alignedQuery = queryAligned;
                blast.alignedRead = readAligned;
                alignments.push_back(blast);
            } catch (const exception& error) {
                throw runtime_error(filename_ + ":" + to_string(lineNumber) + ": " + error.what());
            }
        }
        if (blastFile.bad()) throw runtime_error("Error reading BLAST alignments: " + filename_);
    }

private:
    string filename_;
};

#endif
