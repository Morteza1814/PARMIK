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

    InnerContainerType get(const Key& key) {
        return container[key];
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
    pair<uint64_t, uint64_t> collectCheapSeeds(vector<Key> kmers, Container<Value, Value>& minThCheapSeeds, uint32_t minNumberOfCheapSeeds, Value queryInd) {
        uint64_t cheapSeeds = 0;
        uint64_t cheapSeedReads = 0;
        tsl::robin_map <Value, uint64_t> readCounts;
        for (auto it = kmers.begin(); it != kmers.end(); it++) 
        {
            auto itt = container.find(*it);
            if (itt != container.end()) 
            {
                cheapSeeds++;
                const auto& innerContainer = itt->second;
                for (const auto& value : innerContainer) {
                    auto ittt = readCounts.find(value);
                    if (ittt == readCounts.end()) {
                        readCounts[value] = 1;
                    } else {
                        if (ittt->second == minNumberOfCheapSeeds) {
                            minThCheapSeeds.insertValue(queryInd, value);
                            cheapSeedReads++;
                        }
                        readCounts[value]++;
                    }
                }
            }
        }
        return make_pair(cheapSeeds, cheapSeedReads);
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

        // Serialize the Container to a file
    void serialize(const string& filename) {

        ofstream outfile(filename, ios::binary);

        // Write size of outer map
        size_t outerSize = container.size();
        outfile.write((char*)&outerSize, sizeof(outerSize));

        // Loop through outer map
        for (auto& kv : container) {
            
            // Write key
            outfile.write((char*)&kv.first, sizeof(kv.first));
            
            // Get inner container
            auto& inner = kv.second;
            
            // Write size of inner container
            size_t innerSize = inner.size();
            outfile.write((char*)&innerSize, sizeof(innerSize));
            
            // Write each element of inner container
            for (auto& elem : inner) {
            outfile.write((char*)&elem, sizeof(elem)); 
            }
        }

        outfile.close();
    }


    // Deserialize the Container from a file 
    void deserialize(const string& filename) {
        // OuterContainerType container2;
        ifstream infile(filename, ios::binary); 

        // Read size of outer map
        size_t outerSize;
        infile.read((char*)&outerSize, sizeof(outerSize));

        // Clear existing data
        container.clear();

        // Read each entry
        for (size_t i = 0; i < outerSize; i++) {

            Key key;
            infile.read((char*)&key, sizeof(key));

            // Read size of inner container
            size_t innerSize;
            infile.read((char*)&innerSize, sizeof(innerSize));

            // Read each element
            InnerContainerType inner;
            for (size_t j = 0; j < innerSize; j++) {
            Value elem;
            infile.read((char*)&elem, sizeof(elem));
            inner.insert(elem);
            }

            // Insert into container
            container[key] = inner;
        }
        // cout << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<check serlizer>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        // cout << "deser cheap kmer size is : " << container2.size() << endl;
        // for (const auto& kv : container) {
        //     auto it = container2.find(kv.first);
        //     if (it == container2.end()) {
        //         cout << "Error: key not found in deserialized container" << endl;
        //     } else {
        //         // Check inner values
        //         const auto& originalValues = kv.second;
        //         const auto& deserializedValues = it->second;
                
        //         for (const auto& value : originalValues) {
        //             if (deserializedValues.find(value) == deserializedValues.end()) {
        //                 cout << "Error: value not found in deserialized container" << endl;
        //             } 
        //         }
        //     }
        // }
        // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<check serlizer>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n" << endl;
        infile.close();
    }

// private:
    OuterContainerType container;

};

#endif