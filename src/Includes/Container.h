#ifndef CONTAINER_H
#define CONTAINER_H

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <map>
#include <chrono>
#include "tsl/robin_map.h"
#include "tsl/robin_set.h"

using namespace std;

template<typename Key, typename Value>
class Container {
public:
    using InnerContainerType = unordered_set<Value>;
    // using OuterContainerType = tsl::robin_map<Key, InnerContainerType>;
    using OuterContainerType = unordered_map<Key, InnerContainerType>;

    // Constructor
    Container() = default;

    // Insert values into the outer container
    void insert(const Key& key, const InnerContainerType& innerContainer) {
        container[key] = innerContainer;
    }

    // Insert a value into the inner container associated with the given key
    void insertValue(const Key& key, const Value& value) {
        container[key].insert(value);
    }

    // Select keys from the outer container
    InnerContainerType selectKeys(const vector<Key>& selectedKeys) {
        InnerContainerType selectedInnerContainer;
        for (const auto& key : selectedKeys) {
            const auto& innerContainer = container[key];
            selectedInnerContainer.insert(innerContainer.begin(), innerContainer.end());
        }
        return selectedInnerContainer;
    }

    // void getSize()
    // {
    //     size_t size = 0;
    //     for (const auto& entry : container) {
    //         // Size of the key (assuming int as the key type)
    //         size += sizeof(Key);

    //         // Size of the inner set (number of elements multiplied by size of Value)
    //         size += entry.second.size() * sizeof(Value);
    //     }

    //     // Convert to megabytes and gigabytes
    //     double sizeMB = static_cast<double>(size) / (1024 * 1024); // Convert to MB
    //     double sizeGB = sizeMB / 1024; // Convert to GB

    //     cout << "Size in bytes: " << size << " bytes" << endl;
    //     cout << "Size in megabytes: " << sizeMB << " MB" << endl;
    //     cout << "Size in gigabytes: " << sizeGB << " GB" << endl;

    // }

    // Count the occurrences of distinct values across the selected inner containers
    unordered_map<Value, int> countDistinctValues(const InnerContainerType& selectedInnerContainer) {
        unordered_map<Value, int> valueCounts;
        for (const auto& value : selectedInnerContainer) {
            valueCounts[value]++;
        }
        return valueCounts;
    }

    // Count the occurrences of distinct values across all inner containers
    map<Value, uint32_t> countDistinctValues() {
        map<Value, uint32_t> valueCounts;
        for (const auto& pair : container) {
            const auto& innerContainer = pair.second;
            for (const auto& value : innerContainer) {
                valueCounts[value]++;
            }
        }
        return valueCounts;
    }

    // Count the occurrences of distinct values for a key
    void collectCheapSeeds(IndexContainer<char, Key> kmers, map<Value, uint32_t> &frontReadCounts, map<Value, uint32_t> &backReadCounts) {
        char frontRegion = 'F';
        char backRegion = 'B';
        auto frontRange = kmers.getRange(frontRegion);
        uint32_t frontCheapSeeds = 0, backCheapSeeds = 0;
        for (auto it = frontRange.first; it != frontRange.second; it++) 
        {
            auto itt = container.find(it->second);
            if (itt != container.end()) 
            {
                frontCheapSeeds++;
                const auto& innerContainer = itt->second;
                for (const auto& value : innerContainer) {
                    frontReadCounts[value]++;
                }
            }
        }
        auto backRange = kmers.getRange(backRegion);
        for (auto it = backRange.first; it != backRange.second; it++) 
        {
            backCheapSeeds++;
            auto itt = container.find(it->second);
            if (itt != container.end()) 
            {
                const auto& innerContainer = itt->second;
                for (const auto& value : innerContainer) {
                    backReadCounts[value]++;
                }
            }
        }
        // printf("# of cheap seeds => Front [%d] - Back [%d]\n", frontCheapSeeds, backCheapSeeds);
    }

    void serializeSet(ofstream& ofs, const InnerContainerType& data) {
        size_t size = data.size();
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));

        for (const auto& element : data) {
            ofs.write(reinterpret_cast<const char*>(&element), sizeof(element));
        }
    }

    void deserializeSet(ifstream& ifs, InnerContainerType& data) {
        data.clear();

        size_t size;
        ifs.read(reinterpret_cast<char*>(&size), sizeof(size));

        for (size_t i = 0; i < size; ++i) {
            Value element;
            ifs.read(reinterpret_cast<char*>(&element), sizeof(element));
            data.insert(element);
        }
    }

    void serializeMap(string fileName) {
        ofstream ofs(fileName, ios::binary);
        if (!ofs.is_open()) {
            cerr << "Error opening the file for serialization cheap mers: " << fileName << endl;
            throw runtime_error("Error opening the file for serialization cheap mers");
        }
        size_t mapSize = container.size();
        ofs.write(reinterpret_cast<const char*>(&mapSize), sizeof(mapSize));

        for (const auto& kv : container) {
            ofs.write(reinterpret_cast<const char*>(&kv.first), sizeof(kv.first));
            serializeSet(ofs, kv.second);  // Serialize the unordered_set<Value>
        }
        cout << "ser cheap kmer size is : " << container.size() << endl;
        ofs.close();
    }

    void deserializeMap(string fileName) {
        ifstream ifs(fileName, ios::binary);
        if (!ifs.is_open()) {
            cerr << "Error opening the file for deserialization cheap mers: " << fileName << endl;
            throw runtime_error("Error opening the file for deserialization cheap mers");
        }

        container.clear();

        size_t mapSize;
        ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));

        for (size_t i = 0; i < mapSize; ++i) {
            Key key;
            ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
            deserializeSet(ifs, container[key]);  // Deserialize the unordered_set<Value>
        }
        cout << "deser cheap kmer size is : " << container.size() << endl;
        ifs.close();
    }


// private:
    OuterContainerType container;

};

#endif