#ifndef NUCLEOTIDEENCODER_H
#define NUCLEOTIDEENCODER_H

#include <string>
#include <cstdint>
#include <stdexcept>
#include <iostream>

using namespace std;

template <typename keyT>
class NucleotideEncoder {
public:
    NucleotideEncoder(size_t stringSize) : stringSize_(stringSize) {
        maxStringSize_ = (sizeof(keyT)*8)/2;
        cout << "k-mer max size : " << maxStringSize_ << ", stringSize = " << stringSize << endl;
        if (stringSize  > maxStringSize_) {
            throw invalid_argument("Invalid Nucleotide sequence size");
        }
    }

    keyT encode(const string& nucleotideString) {
        keyT encodedValue = 0;
        size_t numValues = min(nucleotideString.size(), stringSize_);

        for (size_t i = 0; i < numValues; i++) {
            keyT nucleotideValue = encodeNucleotide(nucleotideString[i]);
            encodedValue |= (nucleotideValue << (2 * i));
        }

        return encodedValue;
    }

    string decode(keyT encodedValue) {
        string nucleotideString;
        size_t numValues = min(stringSize_, maxStringSize_);

        for (size_t i = 0; i < numValues; i++) {
            uint16_t nucleotideValue = (encodedValue >> (2 * i)) & 3;
            char nucleotideChar = decodeNucleotide(nucleotideValue);
            nucleotideString += nucleotideChar;
        }

        return nucleotideString;
    }

    uint16_t encodeNucleotide(char nucleotide) {
        switch (nucleotide) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return 0; // or throw an exception for invalid nucleotides
        }
    }

    char decodeNucleotide(uint16_t nucleotideValue) {
        switch (nucleotideValue) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default: return 'A'; // or throw an exception for invalid nucleotide values
        }
    }

private:
    size_t maxStringSize_ = 32;
    size_t stringSize_;
};

#endif // NUCLEOTIDEENCODER_H
