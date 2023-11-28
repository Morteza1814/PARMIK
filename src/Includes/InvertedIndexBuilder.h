#ifndef INVERTEDINDEXBUILDER_H
#define INVERTEDINDEXBUILDER_H

#include "NucleotideEncoder.h"
#include "IndexContainer.h"
#include "IndexFile.h"
#include "Utils.h"
#include "tsl/robin_map.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

template <typename keyT, typename valT>
class InvertedIndexBuilder {
public:
    InvertedIndexBuilder(size_t k, size_t overlapSize) : k_(k), overlapSize_(overlapSize), encoder_(k) {
        bp_jump_size_ = k_ - overlapSize_;
        cout << "bp_jump_size : " << bp_jump_size_ << endl;
    }

    vector<keyT> extractKmers(const string& sequence)
    {
        vector<keyT> kmers;
        if (sequence.size() < k_) {
            return kmers;
        }

        for (size_t i = 0; i <= sequence.size() - k_; i += bp_jump_size_) {
            string kmerString = sequence.substr(i, k_);
            keyT encodedKmer = encoder_.encode(kmerString);
            kmers.push_back(encodedKmer);
        }

        return kmers;
    }

    IndexContainer<keyT, valT> build(const string& fastqFile, uint32_t contigCountToRead, bool isIndexOffline)
    {
        ifstream file(fastqFile);
        if (!file) {
            cout << "fastqFile address: " << fastqFile << endl;
            throw runtime_error("Failed to open the FASTQ file.");
        }
        IndexContainer<keyT, valT> invertedIndex;
        string line;
        uint32_t currentContigId = 0;
        string currentContigSequence;
        uint32_t contigCount = 0;
        uint32_t kmerIndCount = 0;
        uint32_t readsContainingN = 0;
        Utilities<uint32_t> utl;
        cout << "<<<<<<<<<<<<<<<<<<<<<<<< Index Construction Started!!>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        clock_t start_time = clock();
        // if (!isIndexOffline)
        // {
            while (getline(file, line) && contigCount <= contigCountToRead) {
                if (line.empty()) {
                    // Skip empty lines
                    continue;
                }

                if (line[0] == '>') {
                    // Line contains the contig ID
                    currentContigId = utl.extractContigId(line);
                } else {
                    // read should not contain N
                    if (line.find('n') != string::npos || line.find('N') != string::npos)
                    {
                        readsContainingN++;
                        continue;
                    }
                    // Line contains the contig sequence
                    currentContigSequence = line;

                    // Extract k-mers from the contig sequence
                    vector<keyT> kmers = extractKmers(currentContigSequence);

                    // Add k-mers to the inverted index
                    for (const auto& kmer : kmers) {
                        invertedIndex.put(kmer, currentContigId);
                    }
                    contigCount++;
                    kmerIndCount += kmers.size();
                    if ((currentContigId % 100000) == 0)
                    {
                        cout << "The " << (currentContigId / 100000) << " 100k contigs have been processed!" << endl;
                    } 
                }
            }
            printf("Index Construction took: %.2f seconds\n", (double)(clock() - start_time)/CLOCKS_PER_SEC);
            cout << left << setw(40) << "Number of kmers: " << kmerIndCount << endl;
            cout << left << setw(40) << "Number of contigs: " << contigCount << endl;
            cout << left << setw(40) << "Size of inverted index: " << invertedIndex.size() << endl;
            cout << left << setw(40) << "Number of Reads Containing N: " << readsContainingN << endl;
            cout <<  "<<<<<<<<<<<<<<<<<<<<<<<Index Construction Ended!!>>>>>>>>>>>>>>>>>>>>>>> " << endl;
            file.close();
            //store the index
            // ife.saveToFile(invertedIndex);
        // }else
        // {
        //     ife.loadFromFile(invertedIndex);
        //     cout << "invertedIndex loaded size: " << invertedIndex.size() << endl;
        // }
        return invertedIndex;
    }

    void printIndex(IndexContainer<keyT, valT> index)
    {
        for (const auto& pair : index.getContainer()) 
        {
            cout << "Key: " << encoder_.decode(pair.first) << ", Value: " << pair.second << endl;
        }
    }

private:
    size_t k_;
    size_t overlapSize_;
    size_t bp_jump_size_;
    NucleotideEncoder<keyT> encoder_;
};

#endif // INVERTEDINDEXBUILDER_H
