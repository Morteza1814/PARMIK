#ifndef MINIMAP2READER_H
#define MINIMAP2READER_H

#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "SamReader.h"

// Read existing minimap2 SAM records without changing legacy SAM interpretation.
class Minimap2Reader {
    std::string filename_;

    static uint32_t number(const std::string& value, uint32_t maximum) {
        if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
            throw std::runtime_error("expected an unsigned integer: " + value);
        const auto result = std::stoull(value);
        if (result > maximum)
            throw std::runtime_error("integer outside supported range: " + value);
        return static_cast<uint32_t>(result);
    }

    static uint32_t contigId(const std::string& name) {
        const auto dot = name.find_last_of('.');
        return number(name.substr(dot == std::string::npos ? 0 : dot + 1),
                      std::numeric_limits<int>::max());
    }

public:
    explicit Minimap2Reader(const std::string& filename) : filename_(filename) {}

    void parseFile(uint32_t queryCount, std::vector<Alignment>& alignments) const {
        std::ifstream input(filename_);
        if (!input)
            throw std::runtime_error("Cannot open minimap2 SAM: " + filename_);
        std::string line;
        size_t lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            if (line.empty() || line[0] == '@') continue;
            try {
                std::istringstream row(line);
                std::string fields[11];
                for (auto& field : fields) {
                    if (!(row >> field))
                        throw std::runtime_error("expected 11 SAM fields (PAF is unsupported)");
                }
                const auto flag = number(fields[1], 65535);
                if (flag & 0x4) continue; // Unmapped records can retain reference coordinates.
                if (fields[2] == "*" || fields[5] == "*")
                    throw std::runtime_error("mapped record has no reference or CIGAR");
                const auto queryId = contigId(fields[0]);
                // SAM need not be sorted by query ID.
                if (queryId >= queryCount) continue;
                const auto readId = contigId(fields[2]);

                uint32_t ambiguous = 0, exact = 0, substitutions = 0, indels = 0;
                const auto& cigar = fields[5];
                for (size_t i = 0; i < cigar.size();) {
                    const auto start = i;
                    while (i < cigar.size() && cigar[i] >= '0' && cigar[i] <= '9') ++i;
                    if (start == i || i == cigar.size())
                        throw std::runtime_error("invalid CIGAR: " + cigar);
                    const auto count = number(cigar.substr(start, i - start),
                                              std::numeric_limits<int>::max());
                    if (count == 0) throw std::runtime_error("zero-length CIGAR operation");
                    switch (cigar[i++]) {
                        case 'M': ambiguous += count; break;
                        case '=': exact += count; break;
                        case 'X': substitutions += count; break;
                        case 'I': case 'D': indels += count; break;
                        case 'S': case 'H': break;
                        default: throw std::runtime_error("unsupported CIGAR operation: " + cigar);
                    }
                    // The existing Alignment representation uses 8-bit length fields.
                    if (uint64_t(ambiguous) + exact + substitutions + indels > 255)
                        throw std::runtime_error("alignment exceeds PARMIK's 255-column limit");
                }
                const auto length = ambiguous + exact + substitutions + indels;
                if (length == 0) throw std::runtime_error("empty alignment");

                bool hasNM = false;
                uint32_t edits = 0;
                std::string tag;
                while (row >> tag) {
                    if (tag.rfind("NM:", 0) == 0) {
                        if (hasNM || tag.rfind("NM:i:", 0) != 0)
                            throw std::runtime_error("invalid or duplicate NM tag");
                        edits = number(tag.substr(5), length);
                        hasNM = true;
                    }
                }
                if (ambiguous && !hasNM)
                    throw std::runtime_error("M CIGAR requires NM:i edit counts");
                if (hasNM) {
                    if (edits < indels + substitutions ||
                        edits > indels + substitutions + ambiguous)
                        throw std::runtime_error("NM disagrees with CIGAR");
                    substitutions = edits - indels;
                }
                Alignment aln;
                aln.queryID = queryId;
                aln.readID = readId;
                aln.cigar = cigar;
                aln.matches = length - indels - substitutions;
                aln.substitutions = substitutions;
                aln.inDels = indels;
                aln.editDistance = substitutions + indels;
                aln.partialMatchSize = length;
                // Alignment stores strand only; include all mapped secondary and
                // supplementary records as separate candidates, without joining them.
                aln.flag = flag & 0x10;
                alignments.push_back(aln);
            } catch (const std::exception& error) {
                throw std::runtime_error(filename_ + ":" + std::to_string(lineNumber) +
                                         ": " + error.what());
            }
        }
        if (input.bad()) throw std::runtime_error("Error reading minimap2 SAM: " + filename_);
    }
};

#endif
