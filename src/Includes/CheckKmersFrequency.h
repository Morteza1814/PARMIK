#ifndef CHECKKMERSFREQUENCY_H
#define CHECKKMERSFREQUENCY_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include "NucleotideEncoder.h"
#include "Utils.h"
#include "tsl/robin_map.h"

using namespace tsl;
using namespace std;

template <typename kmerT>
class CheckKmerFrequency
{
private:
    NucleotideEncoder<kmerT> encoder_;
    size_t k_;
    char frontRegion = 'F';
    char backRegion = 'B';
    size_t regionSize;
    size_t minExactMatchLen;

public:
    CheckKmerFrequency(size_t k, size_t R, uint32_t m) : encoder_(k), k_(k), regionSize(R), minExactMatchLen(m)
    {
    }

    vector<kmerT> extractReadKmers(const string &sequence)
    {
        vector<kmerT> kmers;
        if (sequence.size() < k_)
        {
            return kmers;
        }

        for (size_t i = 0; i <= sequence.size() - k_; i += 1)
        {
            string kmerString = sequence.substr(i, k_);
            kmerT encodedKmer = encoder_.encode(kmerString);
            kmers.push_back(encodedKmer);
        }

        return kmers;
    }

    void extractKmersFromRegions(const string &sequence, vector<kmerT> &regionKmers, char region)
    {
        if (sequence.size() < k_)
        {
            cout << "sequence : " << sequence << endl;
            throw invalid_argument("less than k sequence size");
        }
        if (region == frontRegion)
        {
            for (size_t i = 0; i <= regionSize - k_; i++)
            {
                string kmerString = sequence.substr(i, k_);
                kmerT encodedKmer = encoder_.encode(kmerString);
                regionKmers.push_back(encodedKmer);
                // cout << "kmer = " << kmerString << endl;
            }
        }
        else if (region == backRegion)
        {
            for (size_t i = (sequence.size() - regionSize); i <= sequence.size() - k_; i++)
            {
                string kmerString = sequence.substr(i, k_);
                // cout << "kmer = " << kmerString << endl;
                kmerT encodedKmer = encoder_.encode(kmerString);
                regionKmers.push_back(encodedKmer);
            }
        }
        else
        {
            cerr << "wrong region char!!" << endl;
        }
    }

    void checkKmerFreq(string missedMatchesFileAddress, string kmerFrequLocFileAddress, string queryFileAddress, string readFileAddress, uint32_t numbQueries, uint32_t numReads)
    {
        tsl::robin_map <uint32_t, string> reads, queries;
        Utilities<uint32_t> util;
        util.readContigsFromFile(queryFileAddress, numbQueries, queries);
        size_t readCount = util.readContigsFromFile(readFileAddress, numReads, reads);
        
        cout << "create index started\n";
        //create the index
        //create the partial matching inverted read index
        InvertedIndexBuilder<uint64_t, uint32_t> builder(k_, k_-1);
        // Build the inverted index
        IndexContainer<uint64_t, uint32_t> invertedIndex = builder.build(readFileAddress, readCount, 0);
        cout << "create index finished\n";

        ifstream inputFile(missedMatchesFileAddress);
        ofstream outputFile(kmerFrequLocFileAddress);
        
        if (!inputFile.is_open() || !outputFile.is_open()) {
            cerr << "Error opening file." << endl;
            return;
        }

        unordered_multimap<uint32_t, pair<uint32_t, uint32_t>> dataMultimap;

        string line;
        while (getline(inputFile, line)) {
            istringstream iss(line);
            uint32_t firstColumn, secondColumn, thirdColumn;
            
            if ((iss >> firstColumn) && (iss >> secondColumn) && (iss >> thirdColumn)) {
                // Store the values in your data structure (here using unordered_multimap)
                dataMultimap.insert({firstColumn, make_pair(secondColumn, thirdColumn)});
            } else {
                cerr << "Error parsing line: " << line << endl;
            }
        }
        inputFile.close();

        for (const auto& pair : dataMultimap) 
        {
            auto queryId = pair.first;
            auto readId = pair.second.first;
            auto strand = pair.second.second;
            string query = queries[queryId]; 
            string read = reads[readId]; 
            if (strand == 16)
                reverse(query.begin(), query.end());
            vector<kmerT> readKmers = extractReadKmers(read);
            //check front region
            vector<kmerT> frontKmers;
            extractKmersFromRegions(query, frontKmers, frontRegion);
            outputFile << "Q: " << queryId << "\t R:" << readId << "\n front region: " << endl;
            uint16_t consecutiveMatchCount = 0;
            for (unsigned int i = 0; i < frontKmers.size(); i++)
            {
                for (unsigned int j = 0; j < readKmers.size(); j++)
                {
                    if (frontKmers[i] == readKmers[j])
                    {
                        consecutiveMatchCount++;
                        uint64_t freq =  invertedIndex.container_.count(frontKmers[i]);
                        outputFile <<"(" << i << ", " << j << ", "<< freq <<"), " << endl;
                        if(consecutiveMatchCount >= minExactMatchLen - k_ + 1)
                            outputFile << ", [min exact match = " << consecutiveMatchCount << "] ";
                    } else 
                    {
                        consecutiveMatchCount = 0;
                    }
                }
                outputFile << endl;
            }
            
            //check back region
            vector<kmerT> backKmers;
            extractKmersFromRegions(query, backKmers, backRegion);
            outputFile << "\n back region: " << endl;
            consecutiveMatchCount = 0;
            for (unsigned int i = 0; i < backKmers.size(); i++)
            {
                for (unsigned int j = 0; j < readKmers.size(); j++)
                {
                    if (backKmers[i] == readKmers[j])
                    {
                        uint64_t freq =  invertedIndex.container_.count(backKmers[i]);
                        outputFile <<"(" << i << ", " << j << ", "<< freq <<"), " << endl;
                        if(consecutiveMatchCount >= minExactMatchLen - k_ + 1)
                            outputFile << ", [min exact match = " << consecutiveMatchCount << "] ";
                    } else
                    {
                        consecutiveMatchCount = 0;
                    }
                }
                outputFile << endl;
            }
            outputFile << "-------------------------------------------------------------------" << endl;
        }
    }
};
#endif