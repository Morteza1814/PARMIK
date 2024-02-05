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
    size_t regionSize;
    uint32_t minNumberOfCheapSeeds;
    char frontRegion = 'F';
    char backRegion = 'B';
    uint32_t *frontNumberOfKmersMatchedInQuery;
    uint32_t *backNumberOfKmersMatchedInQuery;
    bool isVeboseLog_ = false;
public:
    CheapKmerPartialMatcher(size_t k, size_t R, uint32_t minCheaps, bool isVeboseLog) : encoder_(k), k_(k), regionSize(R), minNumberOfCheapSeeds(minCheaps), isVeboseLog_(isVeboseLog) 
    {
        frontNumberOfKmersMatchedInQuery = new uint32_t[regionSize-k_+2];
        backNumberOfKmersMatchedInQuery = new uint32_t[regionSize-k_+2];
        for (size_t i = 0; i <= regionSize-k+1; i++)
        {
            frontNumberOfKmersMatchedInQuery[i] = 0;
            backNumberOfKmersMatchedInQuery[i] = 0;
        }
    }

    ~CheapKmerPartialMatcher()
    {
        delete [] frontNumberOfKmersMatchedInQuery;
        delete[] backNumberOfKmersMatchedInQuery;
    }

    void printArrays() 
    {
        for (size_t i = 0; i <= regionSize-k_; i++) {
            cout << "frontNumberOfKmersMatchedInQuery[" << i << "]: " << frontNumberOfKmersMatchedInQuery[i] << endl;
            cout << "backNumberOfKmersMatchedInQuery[" << i << "]: " << backNumberOfKmersMatchedInQuery[i] << endl;
        }
    }

    string printKeys(const map<readIndT, uint32_t>& map) {
        ostringstream out;
        for (const auto& entry : map) {
            out << entry.first << ", ";
        }
        return out.str();
    }

    void cheapSeedFilter(Container<kmerT, readIndT>& cheapKmers, tsl::robin_map <uint32_t, string>& queries, 
    IndexContainer<queryIndT, readIndT>& frontMinThCheapSeeds, IndexContainer<queryIndT, readIndT>& backMinThCheapSeeds)
    {
        Utilities<uint32_t> utilint;
        Utilities<double> utildouble;
        // uint64_t currentContigId;
        // uint64_t numberOfQueriesWithMinCheapAtLeast = 0, numberOfCheapSeedsInQ = 0;
        string currentContigSequence;
        queryIndT currentContigId;
        cout <<  "<<<<<<<<<<<<<<<<< Partial Matcher based on Cheap Kmers Started!! >>>>>>>>>>>>>>>>>>" << endl;
        clock_t start_time = clock();
        set<uint32_t> queriesFrontNumCheapKmers; 
        set<uint32_t> queriesBackNumCheapKmers;
        set<uint32_t> queriesnumberOfCheapSeedsInMatchedRead;
        set<uint32_t> queriesFrontNumReadsWithAtLeastOneSimilarCheapKmer;
        set<uint32_t> queriesBackNumReadsWithAtLeastOneSimilarCheapKmer;
        set<uint32_t> queriesFrontMinThCheapSeeds;
        set<uint32_t> queriesBackMinThCheapSeeds;
        set<double> queriesPMtime;
        queryIndT queryCount = 0;
        long long extractKmersFromRegions_totalExeTime = 0, collectCheapSeeds_totalExeTime = 0, findMatchedReadWithHihestSeeds_totalExeTime = 0,
        collectReadsWithMinThCheapSeedsNumber_totalExeTime = 0;
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
            IndexContainer<char, kmerT> queryRegionKmers = extractKmersFromRegions(currentContigSequence); 
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            extractKmersFromRegions_totalExeTime += duration.count();

            //only collect the cheap kmers from the kmers in the regions of query
            map<readIndT, uint32_t> frontCheapSeedReads;
            map<readIndT, uint32_t> backCheapSeedReads;
            start = chrono::high_resolution_clock::now();
            cheapKmers.collectCheapSeeds(queryRegionKmers, frontCheapSeedReads, backCheapSeedReads);
            end = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            collectCheapSeeds_totalExeTime += duration.count();
            if (isVeboseLog_) printf("# of reads with at least 1 cheap seed => Front [%ld] - Back [%ld]\n", frontCheapSeedReads.size(), backCheapSeedReads.size());
            //collect Reads With Min Threshold number of Cheap Seeds
            IndexContainer<queryIndT, readIndT> frontMinThCheapSeedsOfQuery;
            IndexContainer<queryIndT, readIndT> backMinThCheapSeedsOfQuery;
            start = chrono::high_resolution_clock::now();
            collectReadsWithMinThCheapSeedsNumber(frontCheapSeedReads, backCheapSeedReads, frontMinThCheapSeedsOfQuery, backMinThCheapSeedsOfQuery, currentContigId);
            end = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            collectReadsWithMinThCheapSeedsNumber_totalExeTime += duration.count();
            if (isVeboseLog_) printf("# of reads with at least min threshold cheap seed => Front [%ld] - Back [%ld]\n", frontMinThCheapSeedsOfQuery.size(), backMinThCheapSeedsOfQuery.size());
            frontMinThCheapSeeds.insertRange(frontMinThCheapSeedsOfQuery.begin(), frontMinThCheapSeedsOfQuery.end());
            backMinThCheapSeeds.insertRange(backMinThCheapSeedsOfQuery.begin(), backMinThCheapSeedsOfQuery.end());
            //collect the read with max number of similar k-mers for Front and Back
            uint32_t frontMaxValue = 0;
            uint32_t frontMaxKey = 0;
            start = chrono::high_resolution_clock::now();
            for (const auto& pair : frontCheapSeedReads) {
                if (pair.second > frontMaxValue) {
                    frontMaxValue = pair.second;
                    frontMaxKey = pair.first;
                }
            }
            frontNumberOfKmersMatchedInQuery[frontMaxValue]++;        
            uint32_t backMaxValue = 0;
            uint32_t backMaxKey = 0;
            for (const auto& pair : backCheapSeedReads) {
                if (pair.second > backMaxValue) {
                    backMaxValue = pair.second;
                    backMaxKey = pair.first;
                }
            }
            backNumberOfKmersMatchedInQuery[backMaxValue]++;
            //collect the read that has max number of matched k-mers between front and back
            map<queryIndT, readIndT> kmerMatchedResults;
            uint32_t maxVal;
            if (frontMaxValue >= backMaxValue)
            {
                kmerMatchedResults[currentContigId] = frontMaxKey;
                maxVal = frontMaxValue;
                if (isVeboseLog_) cout << "Q : " << currentContigId << " => R : " << frontMaxKey << endl;
            } else
            {
                kmerMatchedResults[currentContigId] = backMaxKey;
                maxVal = backMaxValue;
                if (isVeboseLog_) cout << "Q : " << currentContigId << " => R : " << backMaxKey << endl;
            }
            end = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            findMatchedReadWithHihestSeeds_totalExeTime += duration.count();
            auto pm_end = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::nanoseconds>(pm_end - pm_start);
            double queriesPMtime_seconds = static_cast<double>(duration.count()) / 1'000'000'000.0;
            queriesPMtime.insert((double)(queriesPMtime_seconds));
            if (isVeboseLog_) printf("match for this query took: %.3f seconds\n", (double)(queriesPMtime_seconds));
            queryCount++;
            queriesnumberOfCheapSeedsInMatchedRead.insert(maxVal);
            queriesFrontNumReadsWithAtLeastOneSimilarCheapKmer.insert(frontCheapSeedReads.size());
            queriesBackNumReadsWithAtLeastOneSimilarCheapKmer.insert(backCheapSeedReads.size());
            queriesFrontMinThCheapSeeds.insert(frontMinThCheapSeedsOfQuery.size());
            queriesBackMinThCheapSeeds.insert(backMinThCheapSeedsOfQuery.size());      
        }
        printf("Partial Mathing based on Cheap Kmer took: %.2f seconds\n", (double)(clock() - start_time)/CLOCKS_PER_SEC);
        cout << "------------more partial matching statistics-------------" << endl;
        // report statistics
        tuple<uint32_t, uint32_t, uint32_t> queriesNumberOfFrontCheapKmersTuple = utilint.calculateStatistics(queriesFrontNumCheapKmers);
        tuple<uint32_t, uint32_t, uint32_t> queriesNumberOfBackCheapKmersTuple = utilint.calculateStatistics(queriesBackNumCheapKmers);
        printf("Number of cheap seeds per query up to q (%d) => Front [average: %d, median: %d, sum: %d] - Back [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(queriesNumberOfFrontCheapKmersTuple), get<1>(queriesNumberOfFrontCheapKmersTuple), get<2>(queriesNumberOfFrontCheapKmersTuple), get<0>(queriesNumberOfBackCheapKmersTuple), get<1>(queriesNumberOfBackCheapKmersTuple), get<2>(queriesNumberOfBackCheapKmersTuple));
        
        tuple<uint32_t, uint32_t, uint32_t> queriesnumberOfCheapSeedsInMatchedReadTuple = utilint.calculateStatistics(queriesnumberOfCheapSeedsInMatchedRead);
        printf("Number of cheap seeds in matched read per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(queriesnumberOfCheapSeedsInMatchedReadTuple), get<1>(queriesnumberOfCheapSeedsInMatchedReadTuple), get<2>(queriesnumberOfCheapSeedsInMatchedReadTuple));

        tuple<uint32_t, uint32_t, uint32_t> queriesFrontNumReadsWithAtLeastOneSimilarCheapKmerTuple = utilint.calculateStatistics(queriesFrontNumReadsWithAtLeastOneSimilarCheapKmer);
        tuple<uint32_t, uint32_t, uint32_t> queriesBackNumReadsWithAtLeastOneSimilarCheapKmerTuple = utilint.calculateStatistics(queriesBackNumReadsWithAtLeastOneSimilarCheapKmer);
        printf("Number of reads with at least 1 similar cheap seed per query up to q (%d) => Front [average: %d, median: %d, sum: %d] - Back [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(queriesFrontNumReadsWithAtLeastOneSimilarCheapKmerTuple), get<1>(queriesFrontNumReadsWithAtLeastOneSimilarCheapKmerTuple), get<2>(queriesFrontNumReadsWithAtLeastOneSimilarCheapKmerTuple), get<0>(queriesBackNumReadsWithAtLeastOneSimilarCheapKmerTuple), get<1>(queriesBackNumReadsWithAtLeastOneSimilarCheapKmerTuple), get<2>(queriesBackNumReadsWithAtLeastOneSimilarCheapKmerTuple));

        tuple<uint32_t, uint32_t, uint32_t> queriesFrontMinThCheapSeedsTuple = utilint.calculateStatistics(queriesFrontMinThCheapSeeds);
        tuple<uint32_t, uint32_t, uint32_t> queriesBackMinThCheapSeedsTuple = utilint.calculateStatistics(queriesBackMinThCheapSeeds);
        printf("Number of matches with at least min threshold (%d) seeds per query up to q (%d) => Front [average: %d, median: %d, sum: %d] - Back [average: %d, median: %d, sum: %d]\n", minNumberOfCheapSeeds, queryCount, get<0>(queriesFrontMinThCheapSeedsTuple), get<1>(queriesFrontMinThCheapSeedsTuple), get<2>(queriesFrontMinThCheapSeedsTuple), get<0>(queriesBackMinThCheapSeedsTuple), get<1>(queriesBackMinThCheapSeedsTuple), get<2>(queriesBackMinThCheapSeedsTuple));
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

    IndexContainer<char, kmerT> extractKmersFromRegions(const string& sequence) {
        IndexContainer<char, kmerT> kmers;
        if (sequence.size() < k_) {
            return kmers;
        }

        for (size_t i = 0; i <= sequence.size() - k_; i++) {
            char key;
            if (i<=regionSize-k_)
            {
                key = frontRegion;
            } else if (i >= (sequence.size() - regionSize))
            {
                key = backRegion;
            } else
            {
                continue;
            }         
            
            string kmerString = sequence.substr(i, k_);
            kmerT encodedKmer = encoder_.encode(kmerString);
            kmers.put(key, encodedKmer);
        }

        return kmers;
    }

    void countCheapSeedReads(Container<kmerT, readIndT> frontCheapKmersReadList, Container<kmerT, readIndT> backCheapKmersReadList,
     map<readIndT, uint32_t> &frontSimilarReads, map<readIndT, uint32_t> &backSimilarReads)
    {
        for (auto it = frontCheapKmersReadList.begin(); it != frontCheapKmersReadList.end(); it++) 
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
        for (auto it = backCheapKmersReadList.begin(); it != backCheapKmersReadList.end(); it++) 
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

    void collectReadsWithMinThCheapSeedsNumber(map<readIndT, uint32_t> frontCheapSeedReads, map<readIndT, uint32_t> backCheapSeedReads,
    IndexContainer<queryIndT, readIndT> &frontMinThCheapSeedReads, IndexContainer<queryIndT, readIndT> &backMinThCheapSeedReads, queryIndT queryInd)
    {
        for (auto it = frontCheapSeedReads.begin(); it != frontCheapSeedReads.end(); it++)
        {
            if (it->second >= minNumberOfCheapSeeds)
            {
                frontMinThCheapSeedReads.put(queryInd, it->first);
            }
        }
        for (auto it = backCheapSeedReads.begin(); it != backCheapSeedReads.end(); it++)
        {
            if (it->second >= minNumberOfCheapSeeds)
            {
                backMinThCheapSeedReads.put(queryInd, it->first);
            }
        }
    }

};

#endif // CHEAPKMERPARTIALMATCHER_H
