#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <set>
#include <vector>
#include "tsl/robin_map.h"

using namespace std;

template<typename T>
class Utilities {
public:
    tuple<T, T, T> calculateStatistics(const multiset<T>& data) 
    {
        //initial check
        if(data.size() == 0) {
            return make_tuple(0, 0, 0);
        }
        // Calculate sum
        T sum = 0;
        for (T num : data) {
            sum += num;
        }
        // Calculate average
        T average = (data.size() > 0) ? sum / data.size() : 0;

        // Calculate the index of the middle element
        int size = data.size();
        int middleIndex = size / 2;

        // If the size is odd, return the middle element
        T median;
        if (size % 2 != 0) {
            auto it = data.begin();
            advance(it, middleIndex);
            median = *it;
        } 
        // If the size is even, return the average of the two middle elements
        else {
            auto mid = data.begin();
            auto mid_1 = data.begin();
            advance(mid, middleIndex);
            advance(mid_1, middleIndex-1);
            median = ((*mid_1) + (*mid)) / 2.0;
        }

        // Return the results as a tuple
        return make_tuple(average, median, sum);
    }

    pair<T, T> calculateStatistics2(const multiset<T>& data) 
    {
        //initial check
        if(data.size() == 0) {
            return make_pair(0, 0);
        }
        // Calculate sum
        T sum = 0;
        for (T num : data) {
            sum += num;
        }
        // Calculate average
        T average = (data.size() > 0) ? sum / data.size() : 0;

        // Calculate the index of the middle element
        int size = data.size();
        int middleIndex = size / 2;

         // If the size is odd, return the middle element
        T median;
        if (size % 2 != 0) {
            auto it = data.begin();
            advance(it, middleIndex);
            median = *it;
        } 
        // If the size is even, return the average of the two middle elements
        else {
            auto mid = data.begin();
            auto mid_1 = data.begin();
            advance(mid, middleIndex);
            advance(mid_1, middleIndex-1);
            median = ((*mid_1) + (*mid)) / 2.0;
        }

        // Return the results as a tuple
        return make_pair(average, median);
    }

    string reverseComplement(const string& genome) {
        string complement = genome;
        transform(complement.begin(), complement.end(), complement.begin(), [](char c) {
            switch (c) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default: return c;
            }
        });

        reverse(complement.begin(), complement.end());
        return complement;
    }

    tsl::robin_map <uint32_t, string> reverseComplementMapValues(const tsl::robin_map <uint32_t, string>& genomeMap) {
        tsl::robin_map <uint32_t, string> revMap;
        for (auto& pair : genomeMap) {
            string rev = reverseComplement(pair.second);
            revMap[pair.first] = rev;
        }
        return revMap;
    }

    T extractContigId(const string& line) 
    {
        string contigID;
        size_t dotPos = line.find_last_of('.');
        if (dotPos != string::npos) {
            // Add 1 to dotPos to skip the '.' character
            contigID = line.substr(dotPos + 1).c_str();
        } else if (line.find('.') == string::npos) {
            contigID = line.c_str();
        } else {
            cout << line << endl;
            throw runtime_error("Invalid contig ID format : " + line);
        }
        try {
            T val = stoi(contigID);
            return val; 
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument error: " << e.what() << std::endl;
            // Handle the error, maybe return a default value or rethrow
            // depending on your application's logic
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range error: " << e.what() << std::endl;
            // Handle the error, maybe return a default value or rethrow
            // depending on your application's logic
        }
        return 0;
    }

    uint32_t readContigsFromFile(string fileAddress, T numberOfEntriesToRead, tsl::robin_map <uint32_t, string>& queries, uint32_t& baseIndex)
    {
        ifstream file(fileAddress);
        if (!file.is_open()) {
            cerr << "fileAddress : " << fileAddress << endl;
            throw runtime_error("Failed to open the query file.");
        }
        string line;
        T queryCount = 0;
        T currentContigId;
        uint32_t numberOfReadsWithN = 0;
        while (getline(file, line) && queryCount < numberOfEntriesToRead) {
            if (line.empty() || line == "") 
            {
                continue;
            }
            if (line[0] == '>') {
                // Line contains the contig ID
                currentContigId = extractContigId(line);
                if (queryCount == 0) {
                    baseIndex = currentContigId;
                }
            } else 
            {
                queryCount++;
                if (line.find('n') != string::npos || line.find('N') != string::npos)
                {
                    numberOfReadsWithN++;
                }
                // cout << "currentContigId : " << currentContigId << ", line : " << line << endl;
                queries[currentContigId] = line;
                // cout << line << endl;
            }
        }
        cout << "numberOfReadsWithN : " << numberOfReadsWithN << endl;
        return queryCount;
    }

    void parseCigar(const string& cigar, uint8_t& matches, uint8_t& substitutions, uint8_t& inDels, vector<uint16_t> &editLocations) {
        substitutions = inDels = 0;
        int currentPos = 0; // Current position in the read
        int len = cigar.length();
        matches = 0; substitutions = 0; inDels = 0;
        editLocations.clear();
        for (int i = 0; i < len; ++i) {
            // Parse the numeric part of CIGAR operation
            int num = 0;
            while (i < len && isdigit(cigar[i])) {
                num = num * 10 + (cigar[i] - '0');
                ++i;
            }
            // Extract the CIGAR operation
            char op = cigar[i];
            // Perform actions based on CIGAR operation
            switch (op) {
                case '=':
                    // Match or mismatch (substitution)
                    currentPos += num;
                    matches += num;
                    break;
                case 'X':
                    // Substitution
                    substitutions += num;
                    for (int j = 0; j < num; ++j) {
                        editLocations.push_back(currentPos + j);
                    }
                    currentPos += num;
                    break;
                case 'I':
                case 'D':
                    // Insertion
                    inDels += num;
                    for (int j = 0; j < num; ++j) {
                        editLocations.push_back(currentPos + j);
                    }
                    currentPos += num;
                    break;
                case 'S':
                    // Soft clipping
                    // currentPos += num;
                    break;
                default:
                    // Unsupported CIGAR operation
                    cerr << "Unsupported CIGAR operation: " << op << endl;
                    break;
            }
        }
    }

};

#endif