#ifndef BASELINEPARTIALMATCHER_H
#define BASELINEPARTIALMATCHER_H
#include <fstream>
#include <iostream>
#include <map>
#include <chrono>
#include <set>
#include "NucleotideEncoder.h"
#include "InvertedIndexBuilder.h"
#include "KmersFrequencyCounter.h"
#include "IndexContainer.h"
#include "Utils.h"
#include "IndexFile.h"

template<typename queryIndT, typename kmerT, typename readIndT>
class BaseLinePartialMatcher {
private:
    NucleotideEncoder<kmerT> encoder_;
    size_t k_;
    char frontRegion = 'F';
    char backRegion = 'B';
    size_t regionSize;
public:
    BaseLinePartialMatcher(size_t k, size_t R) : encoder_(k), k_(k), regionSize(R) {}

    struct QueryBaseLineSeeds {
        IndexContainer<kmerT, readIndT> frontSeeds;
        IndexContainer<kmerT, readIndT> backSeeds;
    };

    void extractKmersFromRegions(const string& sequence, vector<kmerT> &frontKmers, vector<kmerT> &backKmers) {
        if (sequence.size() < k_) {
            throw invalid_argument("less than k sequence size");
        }

        for (size_t i = 0; i <= sequence.size() - k_; i++) {
            if (i<=regionSize-k_ || (i >= (sequence.size() - regionSize)))
            {
                string kmerString = sequence.substr(i, k_);
                kmerT encodedKmer = encoder_.encode(kmerString);
                
                if (i<=regionSize-k_)
                {
                    frontKmers.push_back(encodedKmer);
                } else if (i >= (sequence.size() - regionSize))
                {
                    backKmers.push_back(encodedKmer);
                }
            }
        }
    }

    QueryBaseLineSeeds collectQueryRegionSeeds(IndexContainer<kmerT, readIndT> baslineInvertedIndex, vector<kmerT> frontKmers, vector<kmerT> backKmers)
    {
        QueryBaseLineSeeds qblm;
        for (auto it = frontKmers.begin(); it != frontKmers.end(); it++) 
        {
            auto range = baslineInvertedIndex.getRange(*it);
            qblm.frontSeeds.insertRange(range.first, range.second);
        }
        for (auto it = backKmers.begin(); it != backKmers.end(); it++) 
        {
            auto range = baslineInvertedIndex.getRange(*it);
            qblm.backSeeds.insertRange(range.first, range.second);
        }
        return qblm;
    }

    void countReadsSimilarToQ(QueryBaseLineSeeds qbls, map<readIndT, uint32_t> &frontSimilarReads, map<readIndT, uint32_t> &backSimilarReads)
    {
        for (auto it = qbls.frontSeeds.begin(); it != qbls.frontSeeds.end(); it++) 
        {
            if (frontSimilarReads.find(it->second) != frontSimilarReads.end())
            {
                frontSimilarReads[it->second]++;
            }
            else
            {
                frontSimilarReads[it->second]=1;
            }
        }
        for (auto it = qbls.backSeeds.begin(); it != qbls.backSeeds.end(); it++) 
        {
            if (backSimilarReads.find(it->second) != backSimilarReads.end())
            {
                backSimilarReads[it->second]++;
            }
            else
            {
                backSimilarReads[it->second]=1;
            }
        }
    }

    void seedFinder(IndexContainer<kmerT, readIndT> baselineInvertedIndex, size_t queryCountToRead, const string& queryFileAddress, 
     IndexContainer<queryIndT, readIndT> &queriesFrontSeeds, IndexContainer<queryIndT, readIndT> &queriesBackSeeds)
    {
        ifstream file(queryFileAddress);
        if (!file) {
            throw runtime_error("Failed to open the query file.");
        }
        string line;
        string currentContigSequence;
        uint32_t currentContigId = 0;
        queryIndT queryCount = 0;
        long long extractKmersFromRegions_totalExeTime = 0, collectQueryRegionSeeds_totalExeTime = 0, collectReadsCountsSimilarToQ_totalExeTime = 0;
        multiset<uint32_t> queriesFrontNumSimilarReads, queriesBackNumSimilarReads;
        cout <<  "<<<<<<<<<<<<<<<<< BaseLine Seed Finding Started!! >>>>>>>>>>>>>>>>>>" << endl;
        clock_t start_time = clock();
        Utilities<queryIndT> utl;
        while (getline(file, line) && queryCount < queryCountToRead) {
            if (line.empty()) {
                // Skip empty lines
                continue;
            }
            // numberOfCheapKmersInQ = 0;
            if (line[0] == '>') {
                // Line contains the contig ID
                currentContigId = utl.extractContigId(line);
            } else {
                // Line contains the contig sequence
                currentContigSequence = line;
                //Extract k-mers from the regions
                vector<kmerT> frontKmers, backKmers;
                auto start = chrono::high_resolution_clock::now();
                extractKmersFromRegions(currentContigSequence, frontKmers, backKmers); 
                auto end = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::seconds>(end - start);
                extractKmersFromRegions_totalExeTime += duration.count();
                //collect reads containing seeds
                start = chrono::high_resolution_clock::now();
                QueryBaseLineSeeds qbls = collectQueryRegionSeeds(baselineInvertedIndex, frontKmers, backKmers);
                end = chrono::high_resolution_clock::now();
                duration = chrono::duration_cast<chrono::seconds>(end - start);
                collectQueryRegionSeeds_totalExeTime += duration.count();
                //count reads containing seeds per query
                start = chrono::high_resolution_clock::now();
                map<readIndT, uint32_t> frontSimilarReads, backSimilarReads;
                countReadsSimilarToQ(qbls, frontSimilarReads, backSimilarReads);
                end = chrono::high_resolution_clock::now();
                duration = chrono::duration_cast<chrono::seconds>(end - start);
                collectReadsCountsSimilarToQ_totalExeTime += duration.count();
                queriesFrontNumSimilarReads.insert(frontSimilarReads.size());
                queriesBackNumSimilarReads.insert(backSimilarReads.size());
                // preparing the seeds for output
                for(auto it = frontSimilarReads.begin(); it != frontSimilarReads.end(); it++)
                {
                    queriesFrontSeeds.put(currentContigId, it->first);
                }
                for(auto it = backSimilarReads.begin(); it != backSimilarReads.end(); it++)
                {
                    queriesBackSeeds.put(currentContigId, it->first);
                }

                queryCount++;
            }
        }
        Utilities<uint32_t> util;
        tuple<uint32_t, uint32_t, uint32_t> queriesFrontNumSimilarReadsTuple = util.calculateStatistics(queriesFrontNumSimilarReads);
        tuple<uint32_t, uint32_t, uint32_t> queriesBackNumSimilarReadsTuple = util.calculateStatistics(queriesBackNumSimilarReads);
        printf("Number of reads with at least 1 similar seed per query up to q(%d) => Front [average: %d, median: %d, mean: %d] - Back [average: %d, median: %d, mean: %d]\n", queryCount, get<0>(queriesFrontNumSimilarReadsTuple), get<1>(queriesFrontNumSimilarReadsTuple), get<2>(queriesFrontNumSimilarReadsTuple), get<0>(queriesBackNumSimilarReadsTuple), get<1>(queriesBackNumSimilarReadsTuple), get<2>(queriesBackNumSimilarReadsTuple));
        printf("BaseLine Seed Finding took: %.2f seconds\n", (double)(clock() - start_time)/CLOCKS_PER_SEC);
        cout << "------------BaseLine Seed Finding timing details-------------" << endl;
        cout << left << setw(40) << "extractKmersFromRegions_totalExeTime = " << extractKmersFromRegions_totalExeTime << endl;
        cout << left << setw(40) << "collectQueryRegionSeeds_totalExeTime = " << collectQueryRegionSeeds_totalExeTime << endl;
        cout << left << setw(40) << "collectReadsCountsSimilarToQ_totalExeTime = " << collectReadsCountsSimilarToQ_totalExeTime << endl;
        cout << "--------------------------------------------------------" << endl;
        cout <<  "<<<<<<<<<<<<<<<<< BaseLine Seed Finding Ended!! >>>>>>>>>>>>>>>>>>" << endl;
    }

    void constructBaseline(const string& readFileAddress, uint32_t readCount, const string& queryFileAddress, uint32_t queryCount, uint16_t minExactMatchSize,
     uint16_t cheapKmerThreshold, IndexContainer<queryIndT, readIndT> &queriesFrontSeeds, IndexContainer<queryIndT, readIndT> &queriesBackSeeds, bool isIndexOffline, string offlineIndexAddress)
    {
        string baseline_R = offlineIndexAddress + "baseline/baseline_R.json";
        string queriesFrontSeedsAddress = offlineIndexAddress + "baseline/Q" + to_string(queryCount) +"/baselinefrontmatches.json";
        string queriesBackSeedsAddress = offlineIndexAddress + "baseline/Q" + to_string(queryCount) +"/baselinebackmatches.json";
        if(!isIndexOffline)
        {
            //first create baseline inverted read index
            cout << "=====================constructing the baseline system=========================" << endl;
            //create the partial matching inverted read index
            InvertedIndexBuilder<kmerT, readIndT>baselineBuilder(minExactMatchSize, minExactMatchSize-1);
            // Build the read DB inverted index
            IndexContainer<kmerT, readIndT> baselineInvertedIndex = baselineBuilder.build(readFileAddress, readCount, isIndexOffline, baseline_R);
            // Print the Size of inverted index
            cout << "baselineInvertedIndex orig size: " << baselineInvertedIndex.size() << endl;
            //count the k-mer frequencies in the baseline
            // KmersFrequencyCounter<kmerT, readIndT> blkfc(cheapKmerThreshold);
            // blkfc.collectCheapKmers(baselineInvertedIndex);
            cout << "==============================================================================" << endl;
            seedFinder(baselineInvertedIndex, queryCount, queryFileAddress, queriesFrontSeeds, queriesBackSeeds);
            IndexFile<queryIndT, readIndT> ifbl_front(queriesFrontSeedsAddress);
            IndexFile<queryIndT, readIndT> ifbl_back(queriesBackSeedsAddress);
            ifbl_front.saveToFile(queriesFrontSeeds);
            ifbl_back.saveToFile(queriesBackSeeds);
        } else
        {
            // IndexContainer<kmerT, readIndT> baselineInvertedIndexloaded;
            // IndexFile<kmerT, readIndT> ifbl(baseline_R);
            // ifbl.loadFromFile(baselineInvertedIndexloaded);
            // cout << "baselineInvertedIndex loaded size: " << baselineInvertedIndexloaded.size() << endl;
            IndexFile<queryIndT, readIndT> ifbl_front(queriesFrontSeedsAddress);
            IndexFile<queryIndT, readIndT> ifbl_back(queriesBackSeedsAddress);
            ifbl_front.loadFromFile(queriesFrontSeeds);
            cout << "queriesFrontSeeds loaded size: " << queriesFrontSeeds.size() << endl;
            ifbl_back.loadFromFile(queriesBackSeeds);
            cout << "queriesBackSeeds loaded size: " << queriesBackSeeds.size() << endl;
        }
    }
};

#endif