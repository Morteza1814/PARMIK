#ifndef CHEAPKMERPARTIALMATCHER_H
#define CHEAPKMERPARTIALMATCHER_H

#include <fstream>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <map>
#include <set>
#include <chrono>
#include "NucleotideEncoder.h"
#include "InvertedIndexBuilder.h"
#include "Utils.h"
#include "Container.h"

using namespace std;

template<typename queryIndT, typename kmerT, typename readIndT>
class CheapKmerPartialMatcher {
private:
    NucleotideEncoder<kmerT> encoder_;
    size_t k_;
    size_t contigSize;
    uint32_t minNumberOfCheapSeeds;
    // uint32_t *numberOfKmersMatchedInQuery;
    bool isVeboseLog_ = false;
public:
    CheapKmerPartialMatcher(size_t k, size_t c, uint32_t minCheaps, bool isVeboseLog) : encoder_(k), k_(k), contigSize(c), minNumberOfCheapSeeds(minCheaps), isVeboseLog_(isVeboseLog) 
    {
        // numberOfKmersMatchedInQuery = new uint32_t[contigSize-k_+2];
        // for (size_t i = 0; i <= contigSize-k+1; i++)
        // {
        //     numberOfKmersMatchedInQuery[i] = 0;
        // }
    }

    ~CheapKmerPartialMatcher()
    {
        // delete [] numberOfKmersMatchedInQuery;
    }

    void printArrays() 
    {
        // for (size_t i = 0; i <= contigSize-k_; i++) {
        //     cout << "numberOfKmersMatchedInQuery[" << i << "]: " << numberOfKmersMatchedInQuery[i] << endl;
        // }
    }

    string printKeys(const map<readIndT, uint32_t>& map) {
        ostringstream out;
        for (const auto& entry : map) {
            out << entry.first << ", ";
        }
        return out.str();
    }

    void cheapSeedFilter(Container<kmerT, readIndT>& cheapKmers, tsl::robin_map <uint32_t, string>& queries, 
    Container<queryIndT, readIndT>& minThCheapSeeds)
    {
        Utilities<uint32_t> utilint;
        Utilities<double> utildouble;
        // uint64_t currentContigId;
        // uint64_t numberOfQueriesWithMinCheapAtLeast = 0, numberOfCheapSeedsInQ = 0;
        string currentContigSequence;
        queryIndT currentContigId;
        cout <<  "<<<<<<<<<<<<<<<<< Pre-filtering based on Cheap Kmers Started!! >>>>>>>>>>>>>>>>>>" << endl;
        clock_t start_time = clock();
        multiset<uint32_t> queriesNumCheapKmers; 
        // set<uint32_t> queriesnumberOfCheapSeedsInMatchedRead;
        // set<uint32_t> queriesNumReadsWithAtLeastOneSimilarCheapKmer;
        multiset<uint32_t> queriesMinThCheapSeeds;
        multiset<double> queriesPMtime;
        queryIndT queryCount = 0;
        long long extractKmersFromRegions_totalExeTime = 0, collectCheapSeeds_totalExeTime = 0, findMatchedReadWithHihestSeeds_totalExeTime = 0;
        // long long collectReadsWithMinThCheapSeedsNumber_totalExeTime = 0;
        for(auto itq = queries.begin(); itq != queries.end(); itq++) {
            currentContigSequence = itq->second;
            if (currentContigSequence.find('n') != string::npos || currentContigSequence.find('N') != string::npos)
                continue;
            currentContigId = itq->first;
            if (isVeboseLog_) printf("===================query [%d]===================\n", currentContigId);
            auto pm_start = chrono::high_resolution_clock::now();
            // Line contains the contig sequence

            //Extract k-mers from the regions
            auto start = chrono::high_resolution_clock::now();
            vector<kmerT> queryRegionKmers = extractKmers(currentContigSequence); 
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            extractKmersFromRegions_totalExeTime += duration.count();

            //only collect the cheap kmers from the kmers in the regions of query
            // tsl::robin_map <readIndT, uint32_t> cheapSeedReads;
            start = chrono::high_resolution_clock::now();
            pair<uint64_t, uint64_t> cheapSeeds = cheapKmers.collectCheapSeeds(queryRegionKmers, minThCheapSeeds, minNumberOfCheapSeeds, currentContigId);
            queriesNumCheapKmers.insert(cheapSeeds.first);
            end = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            collectCheapSeeds_totalExeTime += duration.count();
            cout << "query [" << currentContigId << "] contains " << cheapSeeds.first << " cheap seed k-mers and " << cheapSeeds.second << " cheap seeds (reads candidates) and it took: " << static_cast<double>(duration.count()) / 1'000'000.0 << "ms" << endl;
            // if (isVeboseLog_) printf("# of reads with at least 1 cheap seed => [%ld]\n", cheapSeedReads.size());
            //collect Reads With Min Threshold number of Cheap Seeds
            // IndexContainer<queryIndT, readIndT> minThCheapSeedsOfQuery;
            // start = chrono::high_resolution_clock::now();
            // collectReadsWithMinThCheapSeedsNumber(cheapSeedReads, minThCheapSeedsOfQuery, currentContigId);
            // end = chrono::high_resolution_clock::now();
            // duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            // collectReadsWithMinThCheapSeedsNumber_totalExeTime += duration.count();
            // if (isVeboseLog_) printf("# of reads with at least min threshold cheap seed => [%ld]\n", minThCheapSeedsOfQuery.size());
            // minThCheapSeeds.insertRange(minThCheapSeedsOfQuery.begin(), minThCheapSeedsOfQuery.end());
            //collect the read with max number of similar k-mers for Front and Back
            // uint32_t maxValue = 0;
            // uint32_t maxKey = 0;
            // start = chrono::high_resolution_clock::now();
            // for (const auto& pair : cheapSeedReads) {
            //     if (pair.second > maxValue) {
            //         maxValue = pair.second;
            //         maxKey = pair.first;
            //     }
            // }
            // numberOfKmersMatchedInQuery[maxValue]++;        
            //collect the read that has max number of matched k-mers between front and back
            // map<queryIndT, readIndT> kmerMatchedResults;
            // kmerMatchedResults[currentContigId] = maxKey;
            // end = chrono::high_resolution_clock::now();
            // duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            // findMatchedReadWithHihestSeeds_totalExeTime += duration.count();
            auto pm_end = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::nanoseconds>(pm_end - pm_start);
            double queriesPMtime_seconds = static_cast<double>(duration.count()) / 1'000'000'000.0;
            queriesPMtime.insert((double)(queriesPMtime_seconds));
            if (isVeboseLog_) printf("match for this query took: %.3f seconds\n", (double)(queriesPMtime_seconds));
            queryCount++;
            // queriesnumberOfCheapSeedsInMatchedRead.insert(maxValue);
            // queriesNumReadsWithAtLeastOneSimilarCheapKmer.insert(cheapSeedReads.size());
            queriesMinThCheapSeeds.insert(cheapSeeds.second);
        }
        printf("Pre-filtering based on Cheap Kmer took: %.2f seconds\n", (double)(clock() - start_time)/CLOCKS_PER_SEC);
        cout << "------------more partial matching statistics-------------" << endl;
        // report statistics
        tuple<uint32_t, uint32_t, uint32_t> queriesNumberOfCheapKmersTuple = utilint.calculateStatistics(queriesNumCheapKmers);
        printf("Number of cheap seeds per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(queriesNumberOfCheapKmersTuple), get<1>(queriesNumberOfCheapKmersTuple), get<2>(queriesNumberOfCheapKmersTuple));
        
        // tuple<uint32_t, uint32_t, uint32_t> queriesnumberOfCheapSeedsInMatchedReadTuple = utilint.calculateStatistics(queriesnumberOfCheapSeedsInMatchedRead);
        // printf("Number of cheap seeds in matched read per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(queriesnumberOfCheapSeedsInMatchedReadTuple), get<1>(queriesnumberOfCheapSeedsInMatchedReadTuple), get<2>(queriesnumberOfCheapSeedsInMatchedReadTuple));

        // tuple<uint32_t, uint32_t, uint32_t> queriesNumReadsWithAtLeastOneSimilarCheapKmerTuple = utilint.calculateStatistics(queriesNumReadsWithAtLeastOneSimilarCheapKmer);
        // printf("Number of reads with at least 1 similar cheap seed per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(queriesNumReadsWithAtLeastOneSimilarCheapKmerTuple), get<1>(queriesNumReadsWithAtLeastOneSimilarCheapKmerTuple), get<2>(queriesNumReadsWithAtLeastOneSimilarCheapKmerTuple));

        tuple<uint32_t, uint32_t, uint32_t> queriesMinThCheapSeedsTuple = utilint.calculateStatistics(queriesMinThCheapSeeds);
        printf("Number of matches with at least min threshold (%d) seeds per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", minNumberOfCheapSeeds, queryCount, get<0>(queriesMinThCheapSeedsTuple), get<1>(queriesMinThCheapSeedsTuple), get<2>(queriesMinThCheapSeedsTuple));
        tuple<double, double, double> queriesPMtimeTuple = utildouble.calculateStatistics(queriesPMtime);
        printf("Query time for queries up to q(%d) => [average: %f, median: %f, sum: %f]\n", queryCount, get<0>(queriesPMtimeTuple), get<1>(queriesPMtimeTuple), get<2>(queriesPMtimeTuple));                       
        cout << "--------------------------------------------------------" << endl;
        cout << "------------partial matching timing details-------------" << endl;
        double extractKmersFromRegions_totalExeTime_second = static_cast<double>(extractKmersFromRegions_totalExeTime) / 1'000'000'000.0;
        double collectCheapSeeds_totalExeTime_second = static_cast<double>(collectCheapSeeds_totalExeTime) / 1'000'000'000.0;
        double findMatchedReadWithHihestSeeds_totalExeTime_second = static_cast<double>(findMatchedReadWithHihestSeeds_totalExeTime) / 1'000'000'000.0;
        double extractKmersFromRegions_totalExeTime_millisecond = static_cast<double>(extractKmersFromRegions_totalExeTime) / 1'000'000.0;
        double collectCheapSeeds_totalExeTime_millisecond = static_cast<double>(collectCheapSeeds_totalExeTime) / 1'000'000.0;
        double findMatchedReadWithHihestSeeds_totalExeTime_millisecond = static_cast<double>(findMatchedReadWithHihestSeeds_totalExeTime) / 1'000'000.0;
        cout << left << setw(40) << "extractKmersFromRegions_totalExeTime = " << extractKmersFromRegions_totalExeTime << " nanoseconds, " << extractKmersFromRegions_totalExeTime_second << " seconds, " << extractKmersFromRegions_totalExeTime_millisecond << " milliseconds" << endl;
        cout << left << setw(40) << "collectCheapSeeds_totalExeTime = " << collectCheapSeeds_totalExeTime << " nanoseconds, " << collectCheapSeeds_totalExeTime_second << " seconds, " << collectCheapSeeds_totalExeTime_millisecond << " milliseconds" << endl;
        cout << left << setw(40) << "findMatchedReadWithHihestSeeds_totalExeTime = " << findMatchedReadWithHihestSeeds_totalExeTime << " nanoseconds, " << findMatchedReadWithHihestSeeds_totalExeTime_second << " seconds, " << findMatchedReadWithHihestSeeds_totalExeTime_millisecond << " milliseconds" << endl;
        cout << "--------------------------------------------------------" << endl;
        cout <<  "<<<<<<<<<<<<<<<<< Partial Matcher based on Cheap Kmers  Ended!! >>>>>>>>>>>>>>>>>>" << endl;
    }
    
    vector<kmerT> extractKmers(const string& sequence) {
        vector<kmerT> kmers;
        if (sequence.size() < k_) {
            return kmers;
        }

        for (size_t i = 0; i <= sequence.size() - k_; i++) {
            string kmerString = sequence.substr(i, k_);
            kmerT encodedKmer = encoder_.encode(kmerString);
            kmers.push_back(encodedKmer);
        }

        return kmers;
    }

    void collectReadsWithMinThCheapSeedsNumber( tsl::robin_map <readIndT, uint32_t> cheapSeedReads, IndexContainer<queryIndT, readIndT> &minThCheapSeedReads, queryIndT queryInd)
    {
        for (auto it = cheapSeedReads.begin(); it != cheapSeedReads.end(); it++)
        {
            if (it->second >= minNumberOfCheapSeeds)
            {
                minThCheapSeedReads.put(queryInd, it->first);
            }
        }
    }

};

#endif // CHEAPKMERPARTIALMATCHER_H
