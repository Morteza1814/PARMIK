#ifndef SAMREADER_H
#define SAMREADER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <tuple>
#include "Utils.h"
#include "Alignment.h"

using namespace std;

class SamReader {
public:
    struct Sam {
        uint32_t queryId = 0;
        int flag = 0;
        uint32_t readId = 0;
        int pos = 0;
        int mapQ = 0;
        string cigar;
        string rNext;
        string pNext;
        int tLen = 0;
        string seq;
        string qual;
        unsigned int editDistance = 0; // NM flag
        string mismatchPositions; // MD flag
        int alignmentScore = 0; // AS flag
    };

    struct EditCounts {
        uint32_t matchCount = 0;
        uint32_t subCount = 0;
        uint32_t insertCount = 0;
        uint32_t deleteCount = 0;
        uint32_t clipCount = 0;
        uint32_t ambiguousMatchCount = 0; // M columns, unlike explicit = columns
    };

    SamReader(const string& filename) : filename_(filename) {}

    EditCounts extractEditTypes(const string& cigar) const {
        if (cigar.empty() || cigar == "*")
            throw runtime_error("mapped record has no CIGAR");
        EditCounts counts;
        const auto add = [](uint32_t& total, uint32_t count) {
            if (uint64_t(total) + count > numeric_limits<uint32_t>::max())
                throw runtime_error("CIGAR count exceeds supported range");
            total += count;
        };
        for (size_t i = 0; i < cigar.size();) {
            const auto start = i;
            while (i < cigar.size() && cigar[i] >= '0' && cigar[i] <= '9') ++i;
            if (start == i || i == cigar.size())
                throw runtime_error("invalid CIGAR: " + cigar);
            const auto length = number(cigar.substr(start, i - start),
                                       numeric_limits<uint32_t>::max());
            if (length == 0) throw runtime_error("zero-length CIGAR operation");
            switch (cigar[i++]) {
                case 'M':
                    add(counts.ambiguousMatchCount, length);
                    add(counts.matchCount, length);
                    break;
                case '=': add(counts.matchCount, length); break;
                case 'X': case 'x': add(counts.subCount, length); break;
                case 'I': case 'i': add(counts.insertCount, length); break;
                case 'D': case 'd': add(counts.deleteCount, length); break;
                case 'S': case 'H': add(counts.clipCount, length); break;
                default: throw runtime_error("unsupported CIGAR operation: " + cigar);
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

    void parseFile(uint32_t lastQueryID, vector<Alignment>& alignments,
                   bool compressedFormat = false, bool isBWA = false) {
        readFile(lastQueryID, compressedFormat, isBWA, [&](const Alignment& aln) {
            alignments.push_back(aln);
        });
    }

    void parseFileBaseline(uint32_t lastQueryID, vector<Alignmen_Baseline>& alignments,
                           bool compressedFormat = false, bool isBWA = false) {
        readFile(lastQueryID, compressedFormat, isBWA, [&](const Alignment& source) {
            Alignmen_Baseline aln;
            aln.queryID = source.queryID;
            aln.readID = source.readID;
            aln.cigar = source.cigar;
            aln.matches = source.matches;
            aln.substitutions = source.substitutions;
            aln.inDels = source.inDels;
            aln.flag = source.flag;
            alignments.push_back(aln);
        });
    }

    std::tuple<int, int, int, int> countCigarOperations(const std::string& cigar) {
        const auto counts = extractEditTypes(cigar);
        for (auto count : {counts.insertCount, counts.deleteCount,
                           counts.matchCount, counts.subCount}) {
            if (count > static_cast<uint32_t>(numeric_limits<int>::max()))
                throw runtime_error("CIGAR count exceeds supported range");
        }
        return make_tuple(static_cast<int>(counts.insertCount),
                          static_cast<int>(counts.deleteCount),
                          static_cast<int>(counts.matchCount),
                          static_cast<int>(counts.subCount));
    }

    bool isReverseComplemented(const Sam& sam) const {
        return (sam.flag & 0x10) != 0; // Check if the 5th bit (0x10) is set
    }

private:
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

    template<typename Emit>
    void readFile(uint32_t lastQueryID, bool compressedFormat, bool isBWA, Emit emit) {
        ifstream input(filename_);
        if (!input) throw runtime_error("Cannot open SAM input: " + filename_);
        string line;
        size_t lineNumber = 0;
        while (getline(input, line)) {
            ++lineNumber;
            if (line.empty() || line[0] == '@' ||
                line.find_first_not_of(" \t\r") == string::npos) continue;
            try {
                istringstream row(line);
                string fields[11];
                const unsigned fieldCount = compressedFormat ? 5 : 11;
                for (unsigned i = 0; i < fieldCount; ++i) {
                    if (!(row >> fields[i]))
                        throw runtime_error("expected " + to_string(fieldCount) + " SAM fields");
                }
                const auto flag = number(fields[1], 65535);
                if (flag & 0x4) continue; // Unmapped records can retain reference coordinates.
                const auto& cigar = fields[compressedFormat ? 4 : 5];
                if (fields[2] == "*" || cigar == "*")
                    throw runtime_error("mapped record has no reference or CIGAR");
                const auto queryId = contigId(fields[0]);
                // The API takes an inclusive maximum; input need not be query-sorted.
                if (queryId > lastQueryID) continue;
                const auto readId = contigId(fields[2]);
                const auto counts = extractEditTypes(cigar);
                const uint64_t indels = uint64_t(counts.insertCount) + counts.deleteCount;
                const uint64_t length = uint64_t(counts.matchCount) + counts.subCount + indels;
                if (length == 0 || length > numeric_limits<uint8_t>::max())
                    throw runtime_error("alignment must contain 1 to 255 columns");

                bool hasNM = false;
                uint32_t edits = 0;
                string mismatchPositions, token;
                while (row >> token) {
                    if (token.rfind("NM:", 0) == 0) {
                        if (hasNM || token.rfind("NM:i:", 0) != 0)
                            throw runtime_error("invalid or duplicate NM tag");
                        edits = number(token.substr(5), static_cast<uint32_t>(length));
                        hasNM = true;
                    } else if (token.rfind("MD:", 0) == 0) {
                        if (token.rfind("MD:Z:", 0) != 0)
                            throw runtime_error("invalid MD tag");
                        mismatchPositions = token.substr(5);
                    }
                }
                if (isBWA && counts.ambiguousMatchCount && !hasNM)
                    throw runtime_error("M CIGAR requires NM:i edit counts for Bowtie/BWA");
                uint32_t substitutions = counts.subCount;
                if (isBWA && hasNM) {
                    // NM includes inserted/deleted bases, not only substitutions.
                    if (edits < indels + counts.subCount ||
                        edits > indels + counts.subCount + counts.ambiguousMatchCount)
                        throw runtime_error("NM disagrees with CIGAR");
                    substitutions = edits - indels;
                }
                // Legacy PARMIK/baseline writers use a substitutions-only NM
                // tag; their explicit CIGAR remains authoritative for counts.
                Alignment aln;
                aln.queryID = queryId;
                aln.readID = readId;
                aln.cigar = cigar;
                aln.matches = length - indels - substitutions;
                aln.substitutions = substitutions;
                aln.inDels = indels;
                aln.editDistance = substitutions + indels;
                aln.partialMatchSize = length;
                aln.flag = flag & 0x10; // Alignment stores strand, not all SAM flag bits.
                aln.mismatchPositions = mismatchPositions;
                emit(aln);
            } catch (const exception& error) {
                throw runtime_error(filename_ + ":" + to_string(lineNumber) + ": " + error.what());
            }
        }
        if (input.bad()) throw runtime_error("Error reading SAM input: " + filename_);
    }

    string filename_;
};

#endif
