#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <set>
#include "tsl/robin_map.h"

using namespace std;

template<typename T>
class Utilities {
public:
    tuple<T, T, T> calculateStatistics(const set<T>& data) 
    {
        // Calculate sum
        T sum = 0;
        for (T num : data) {
            sum += num;
        }
        // Calculate average
        T average = (data.size() > 0) ? sum / data.size() : 0;

        // Calculate median
        typename set<T>::const_iterator it = data.begin();
        advance(it, data.size() / 2);
        T median = *it;

        // Calculate mean
        T mean = (data.size() > 0) ? sum / data.size() : 0;

        // Return the results as a tuple
        return make_tuple(average, median, mean);
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
        size_t dotPos = line.find_last_of('.');
        if (dotPos != string::npos) {
            // Add 1 to dotPos to skip the '.' character
            return atoi(line.substr(dotPos + 1).c_str());
        } else if (line.find('.') == string::npos)
        {
            return atoi(line.c_str());
        } else {
            cout << line << endl;
            throw runtime_error("Invalid contig ID format : " + line);
        }
    }

    uint32_t readContigsFromFile(string fileAddress, T numberOfEntriesToRead, tsl::robin_map <uint32_t, string>& queries)
    {
        ifstream file(fileAddress);
        if (!file) {
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

};

#endif