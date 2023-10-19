#ifndef INDEXFILE_H
#define INDEXFILE_H
#include <iostream>
#include "nlohmann/json.hpp"
#include <fstream>

using json = nlohmann::json;

using namespace std;

template <typename KeyType, typename ValueType>
class IndexFile {
public: 
    IndexFile(string fileAddress)
    {
        // fstream f(fileAddress);
        cout << "json address : " <<  fileAddress << endl;
        // if (!f.is_open())
        // {
        //     cout << "Invalid file name!" << endl;
        // }
        // f.close(); 
        address = fileAddress;
    }
    void saveToFile(IndexContainer<KeyType, ValueType> &container)
    {
        ofstream o(address);
        if (!o.is_open())
        {
            throw runtime_error("Invalid save file name!");
        } 
        auto start = chrono::high_resolution_clock::now();
        json j(container.container_);
        o << j.dump();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(end - start);
        cout << "---------------\nwrite to file took : " << duration.count() << endl;
    }

    void loadFromFile(IndexContainer<KeyType, ValueType> &container)
    {
        ifstream i(address);
        if (!i.is_open())
        {
            throw runtime_error("Invalid load file name!");
        } 
        json j;
        auto start = chrono::high_resolution_clock::now();
        try
        {
            // json j = json::parse(i);
            i >> j;
        }
        catch (json::parse_error& e)
        {
            cout << e.what() << std::endl;
        }
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(end - start);
        cout << "loading json from file finished, json size: " << j.size() << " and took : " << duration.count() << endl;
        
        start = chrono::high_resolution_clock::now();
        for (const auto& entry : j.items()) {
            KeyType k = 0;
            ValueType v = 0;
            for (json::iterator it = entry.value().begin(); it != entry.value().end(); it++)
            {
                if (it == entry.value().begin())
                {
                    k = it->get<KeyType>();
                } else
                {
                    v = it->get<ValueType>();
                }
            }
            container.put(k, v);
        }
        
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::seconds>(end - start);
        cout << "read from json object took : " << duration.count() << endl;
    }
    private:
        string address;
};

#endif