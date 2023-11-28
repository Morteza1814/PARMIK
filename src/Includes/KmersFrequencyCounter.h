#ifndef KMERFREQUENCYCOUNTER_H
#define KMERFREQUENCYCOUNTER_H

#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include "IndexContainer.h"
#include "Container.h"
#include <iomanip>

template <typename keyT, typename valT>
class KmersFrequencyCounter {
private:
    size_t cheapKmersThreshold_; 
public:

    KmersFrequencyCounter(size_t cheapKmersThreshold){
        cheapKmersThreshold_ = cheapKmersThreshold;
    }

    void collectCheapKmers(Container<keyT, valT>& cheapKmers, IndexContainer<keyT, valT>& invertedIndex) {
        // Count the number of items with values less than 500 and more than 500
        size_t cheapKmersCount = 0;
        size_t expensiveKmersCount = 0; 
        IndexContainer<keyT, valT> expensiveKmers;
        cout << "<<<<<<<<<<<<<<<<<<<<<<<<Cheap Kmers Collection Started!!>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        cout << left << setw(30) << "cheap k-mer threshold: " << cheapKmersThreshold_ << endl;
        clock_t start_time = clock();
        size_t totalkmers = 0, uniqueKmers = 0, distinctKmers = 0;
        for (auto it = invertedIndex.begin(); it != invertedIndex.end(); ) {
            // Get the range of items with the current key  
                      
            auto range = invertedIndex.getRange(it->first);
            auto rangeBegin = range.first;
            auto rangeEnd = range.second;
            // Get the number of items for the current key
            size_t itemCount = distance(range.first, range.second);
            // count unique k-mers
            if (itemCount == 1)
            {
                uniqueKmers++;
            }
            // count distinct k-mers
            distinctKmers++;
            // count total k-mers
            totalkmers += itemCount;
            // count cheap k-mers
            if (itemCount < cheapKmersThreshold_ || cheapKmersThreshold_ == 0)
            {
                cheapKmersCount++;
                // cheapKmers.insertRange(rangeBegin, rangeEnd);
                for (auto itt=range.first; itt != range.second; itt++)
                {
                    cheapKmers.insertValue(itt->first, itt->second);

                }
                
            } else
            {
                expensiveKmersCount++;
                expensiveKmers.insertRange(rangeBegin, rangeEnd);
            }
            
            // Move the iterator to the next unique key
            it = range.second;
        }
        //get the most expensive k-mer
        uint32_t mostExpensiveKmer = 0;
        for (auto it = expensiveKmers.begin(); it != expensiveKmers.end(); ) {
            auto range = expensiveKmers.getRange(it->first);
            size_t itemCount = distance(range.first, range.second);
            if (mostExpensiveKmer < itemCount)
            {
                mostExpensiveKmer = itemCount;
            }
            it = range.second;
        }
        cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Collection Time: " << (double)(clock() - start_time)/CLOCKS_PER_SEC << " seconds" << endl;
        cout << left << setw(30) << "Total k-mers: " << totalkmers << endl;
        cout << left << setw(30) << "Distinct k-mers: " << distinctKmers << " [" << (distinctKmers*100) / totalkmers << "\% of totalkmers]" << endl;
        cout << left << setw(30) << "Unique k-mers: " << uniqueKmers << " [" << (uniqueKmers*100) / distinctKmers << "\% of distinctKmers]" << endl;
        cout << left << setw(30) << "Cheap k-mers: " << cheapKmersCount  << " [" << (cheapKmersCount*100) / distinctKmers << "\% of distinctKmers]" << endl;
        cout << left << setw(30) << "Expensive k-mers: " << expensiveKmersCount << " [" << (expensiveKmersCount*100) / distinctKmers << "\% of distinctKmers]" << endl;
        cout << left << setw(30) << "Most expensive k-mer has: " << mostExpensiveKmer << " read ptr\n";
        cout << "<<<<<<<<<<<<<<<<<<<<<<<<CollectCheapKmers Ended!!>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    }
};

#endif // KMERFREQUENCYCOUNTER_H