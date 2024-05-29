#ifndef EXPKMEREVAL_H
#define EXPKMEREVAL_H

#include <iostream>
#include <vector>
#include "SamReader.h"
#include "Alignment.h"
#include "Utils.h"
#include "Container.h"
#include "NucleotideEncoder.h"
#include "tsl/robin_map.h"

using namespace std;

class ExpensiveKmerFNEvalutor {
public:

    // Function to read data from the file and fill the containers with integers
    void readDataFromFile(const std::string& filename, std::set<uint32_t>& queryIds, std::set<uint32_t>& readIds, std::multiset<uint32_t>& inexpKmersSet, std::multiset<uint32_t>& expKmersSet) {
        std::ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string currentContigIdStr, readIdStr, inexpKmers, expKmers;

            if (!(iss >> currentContigIdStr >> readIdStr >> inexpKmers >> expKmers)) {
                std::cerr << "Error parsing line: " << line << std::endl;
                continue; // Skip malformed lines
            }

            try {
                uint32_t currentContigId = std::stoi(currentContigIdStr);
                uint32_t readId = std::stoi(readIdStr);
                uint32_t inexpKmersInt = std::stoi(inexpKmers);
                uint32_t expKmersInt = std::stoi(expKmers);
                inexpKmersSet.insert(inexpKmersInt);
                expKmersSet.insert(expKmersInt);
                queryIds.insert(currentContigId);
                readIds.insert(readId);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid integer conversion: " << line << std::endl;
                continue; // Skip lines with invalid integer conversion
            } catch (const std::out_of_range& e) {
                std::cerr << "Integer out of range: " << line << std::endl;
                continue; // Skip lines with out of range integer
            }
        }

        infile.close();
    }

    vector<uint32_t> extractKmers(const string& sequence, uint32_t k_) {
        NucleotideEncoder<uint32_t> encoder_(k_);
        vector<uint32_t> kmers;
        if (sequence.size() < k_) {
            return kmers;
        }

        for (size_t i = 0; i <= sequence.size() - k_; i++) {
            string kmerString = sequence.substr(i, k_);
            uint32_t encodedKmer = encoder_.encode(kmerString);
            kmers.push_back(encodedKmer);
        }

        return kmers;
    }

    void expensiveKmerFNEvalutor(string parmikAlnFileAddress, const uint32_t queryCount, const uint32_t readsCount, string expKmersFNFileAddress, string offlineExpensiveIndexAddress, string readDatabaseAddress, uint16_t k, string outputDir){
        string outputFileAddress = outputDir + "/expensiveKmerFNEvalutor.txt";
        ofstream outputFile(outputFileAddress);
        Utilities<uint32_t> util;
        tsl::robin_map <uint32_t, string> reads;
        uint32_t readBaseIndex = 0;
        util.readContigsFromFile(readDatabaseAddress, readsCount, reads, readBaseIndex);
        Container<uint32_t, uint32_t> expensiveKmers;
        expensiveKmers.deserialize(offlineExpensiveIndexAddress);
        vector<Alignment> parmikAlignments;
        SamReader parmikSam(parmikAlnFileAddress);
        parmikSam.parseFile(queryCount, parmikAlignments, false);
        std::set<uint32_t> fnQueryIds;
        std::set<uint32_t> fnReadIds;
        std::multiset<uint32_t> inexpKmersSet;
        std::multiset<uint32_t> expKmersSet;
        readDataFromFile(expKmersFNFileAddress, fnQueryIds, fnReadIds, inexpKmersSet, expKmersSet);
        outputFile << "fnQueryIds.size(): " << fnQueryIds.size() << endl;
        outputFile << "fnReadIds.size(): " << fnReadIds.size() << endl;
        std::set<uint32_t> tpQueryIds;
        std::set<uint32_t> tpReadIds;
        for (const auto& alignment : parmikAlignments) {
            uint32_t queryId = alignment.queryID;
            uint32_t readId = alignment.readID;
            tpQueryIds.insert(queryId);
            tpReadIds.insert(readId);
        }
        outputFile << "tpQueryIds.size(): " << tpQueryIds.size() << endl;
        outputFile << "tpReadIds.size(): " << tpReadIds.size() << endl;
        uint32_t fnQueryiesNotFoundInTP = 0;
        for (const auto& queryId : fnQueryIds){
            if (tpQueryIds.find(queryId) == tpQueryIds.end()){
                outputFile << "fn QueryId: " << queryId << " is not present in tpQueryIds file." << endl;
                fnQueryiesNotFoundInTP++;
            }
        }
        outputFile << "fnQueryiesNotFoundInTP: " << fnQueryiesNotFoundInTP << endl;
        uint32_t fnReadsNotFoundInTP = 0;
        for (const auto& readId : fnReadIds){
            if (tpReadIds.find(readId) == tpReadIds.end()){
                outputFile << "fn [" << readId;
                fnReadsNotFoundInTP++;
                string read = reads[readId];
                outputFile << "]: " << read << endl;
                vector<uint32_t> kmers = extractKmers(read, k);
                outputFile << "expensive kmer: ";
                for (const auto& kmer : kmers){
                    unordered_set<uint32_t> expKmerSet = expensiveKmers.get(kmer);
                    
                    if (expKmerSet.size() > 0){
                        NucleotideEncoder<uint32_t> decoder_(k);
                        outputFile << decoder_.decode(kmer) << ", ";
                    }
                }
                outputFile << "\n------------------------------------" << endl;
            }
        }
        outputFile << "fnReadsNotFoundInTP: " << fnReadsNotFoundInTP << endl;
        pair<uint32_t, uint32_t> inexpKmersSetPair = util.calculateStatistics2(inexpKmersSet);
        pair<uint32_t, uint32_t> expKmersSetPair = util.calculateStatistics2(expKmersSet);
        outputFile << left << setw(80) << "(avg, median) No. of inexpensive kmers in FN reads: " << inexpKmersSetPair.first << ", " << inexpKmersSetPair.second << endl;
        outputFile << left << setw(80) << "(avg, median) No. of expensive kmers in FN reads: " << expKmersSetPair.first << ", " << expKmersSetPair.second << endl;
    }

};

#endif