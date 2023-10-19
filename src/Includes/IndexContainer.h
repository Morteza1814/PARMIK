#ifndef INDEXCONTAINER_H
#define INDEXCONTAINER_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>

using namespace std;

template <typename KeyType, typename ValueType>
class IndexContainer {
public: 
    // Type aliases
    using ContainerType = std::unordered_multimap<KeyType, ValueType>;
    using Iterator = typename ContainerType::iterator;
    using ConstIterator = typename ContainerType::const_iterator;
    // Constructor
    IndexContainer() = default;

    void put(const KeyType& key, const ValueType& value) {
        container_.insert(make_pair(key, value));
    }

    set<ValueType> get(const KeyType& key) {
        set<ValueType> result;
        auto range = container_.equal_range(key);
        for (auto it = range.first; it != range.second; it++) {
            result.insert(it->second);
        }
        return result;
    }

    // void getSize()
    // {
    //     // Calculate the size in bytes
    //     std::size_t sizeBytes = 0;

    //     for (const auto& entry : container_) {
    //         // Size of the key (assuming int as the key type)
    //         sizeBytes += sizeof(KeyType);
    //         // Size of the value (assuming int as the value type)
    //         sizeBytes += sizeof(ValueType);
    //     }

    //     // Convert to megabytes and gigabytes
    //     double sizeMB = static_cast<double>(sizeBytes) / (1024 * 1024); // Convert to MB
    //     double sizeGB = sizeMB / 1024; // Convert to GB

    //     std::cout << "Size in bytes: " << sizeBytes << " bytes" << std::endl;
    //     std::cout << "Size in megabytes: " << sizeMB << " MB" << std::endl;
    //     std::cout << "Size in gigabytes: " << sizeGB << " GB" << std::endl;
    // }

    // Get a range of values for a key
    std::pair<Iterator, Iterator> getRange(const KeyType& key) {
        return container_.equal_range(key);
    }

    // Insert a range of key-value pairs
    template<typename InputIterator>
    void insertRange(InputIterator rangeBegin, InputIterator rangeEnd) {
        container_.insert(rangeBegin, rangeEnd);
    }

    // Begin iterator
    Iterator begin() {
        return container_.begin();
    }

    // End iterator
    Iterator end() {
        return container_.end();
    }

    // Const begin iterator
    ConstIterator cbegin() const {
        return container_.cbegin();
    }

    // Const end iterator
    ConstIterator cend() const {
        return container_.cend();
    }

    // Clear the container
    void clear() {
        container_.clear();
    }

    size_t size() const {
        return container_.size();
    }

// private:
    unordered_multimap<KeyType, ValueType> container_;
};

#endif // INDEXCONTAINER_H
