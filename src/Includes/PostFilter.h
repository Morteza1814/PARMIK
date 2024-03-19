#ifndef POSTFILTER_H
#define POSTFILTER_H

#include <iostream>
#include "Aligner.h"

typedef struct Alignment;

using namespace std;

class PostFilter {
private:
    uint32_t regionSize;
    uint32_t allowedEditDistance;
    uint32_t contigSize;
    uint32_t minExactMatchLength;
    double identityPercentage;
public: 
    PostFilter(uint32_t R, uint32_t a, uint32_t c, uint32_t m, double i) : regionSize(R), allowedEditDistance(a), contigSize(c), minExactMatchLength(m), identityPercentage(i) {}

    bool hasMinConsecutiveMatches(const string& cigarStr) {
        uint32_t consecutiveMatchCount = 0;
    
        // Ensure that both strings have the same length for exact matching
        if (cigarStr.length() == 0) {
            cerr << "Error: Sequences have different lengths." << endl;
            return false;
        }
        // cout << "cigarStr: " << cigarStr << endl;
        // Iterate through each character in the strings and count consecutive exact matches
        uint32_t startPos = 0, numOfDeletionsSofare = 0;
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == '=') {
                consecutiveMatchCount++;
                if (consecutiveMatchCount >= minExactMatchLength) {
                    if ((startPos + minExactMatchLength <= regionSize) || ((startPos >= contigSize - regionSize) && (startPos + minExactMatchLength <= contigSize))) {
                        // cout << "11startPos: " << startPos << " consecutiveMatchCount: " << consecutiveMatchCount << endl;
                        return true;  // Found enough consecutive matches inth front or back region
                    }else{
                        startPos++;
                        consecutiveMatchCount--;
                        // cout << "22startPos: " << startPos  << " consecutiveMatchCount: " << consecutiveMatchCount << endl;
                    }
                }
            } else {
                consecutiveMatchCount = 0;  // Reset count if consecutive match is broken
                if (cigarStr[i] == 'D') {
                    numOfDeletionsSofare++;
                }
                startPos = i+1-numOfDeletionsSofare;
            }
        }

        // Return false if the number of consecutive exact matches is less than minConsecutiveMatch
        return false;
    }

    uint32_t getMatchesCount(string cigarStr) {
        uint32_t cnt = 0;
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == '=') {
                cnt++;
            }
        }
        return cnt;
    }

    bool checkIdentityPercentange(string cigarStr){
        uint32_t matches = getMatchesCount(cigarStr);
        uint32_t len = cigarStr.length();
        double identity = matches * 100 / len;
        if(identity < identityPercentage)
            return false;
        return true;
    }

    string convertCigarToStr(const string& cigar) {
        stringstream simplifiedCigar;
        int len = cigar.length();

        for (int i = 0; i < len; ++i) {
            // Parse the numeric part of CIGAR operation
            int num = 0;
            while (i < len && isdigit(cigar[i])) {
                num = num * 10 + (cigar[i] - '0');
                ++i;
            }

            // Extract the CIGAR operation
            char op = cigar[i];

            // Perform actions based on CIGAR operation
            switch (op) {
                case '=':
                    // Match or mismatch (substitution)
                    simplifiedCigar << string(num, '=');
                    break;
                case 'X':
                    // Match or mismatch (substitution)
                    simplifiedCigar << string(num, 'X');
                    break;
                case 'I':
                    // Insertion
                    simplifiedCigar << string(num, 'I');
                    break;
                case 'D':
                    // Deletion
                    simplifiedCigar << string(num, 'D');
                    break;
                case 'S':
                    // Soft clipping
                    simplifiedCigar << string(num, 'S');
                    break;
                default:
                    // Handle other CIGAR operations if needed
                    break;
            }
        }

        return simplifiedCigar.str();
    }

    string convertStrToCigar(const string cigarStr, uint32_t queryS, uint32_t queryE) {
        stringstream cigar;
        char prev = cigarStr[0];
        uint16_t num = 0;
        // cout << "queryS: " << queryS << ", queryE: " << queryE << endl;
        if (prev != 'S' && prev != 'H' && queryS > 0) {
            cigar << queryS << 'S';
        }
        for (size_t i = 1; i < cigarStr.length(); ++i) {
            char cur = cigarStr[i];
            num++;
            if (prev != cur)
            {
                if (prev == 'S' || prev == 'H') {
                    num += queryS;
                }
                cigar << num << prev;
                num = 0;
                prev = cur;
            }
        }
        num++;
        int remainedClips = contigSize - queryE - 1;
        if (remainedClips > 0 && (prev == 'S' || prev == 'H')) {
            num += remainedClips;
            remainedClips = 0;
        }
        cigar << num << prev;
        if (remainedClips > 0) {
            cigar << remainedClips << 'S';
        }
        return cigar.str();
    }

    bool checkAndUpdateBasedOnAlingmentCriteria(Alignment &aln)
    {
        // covert cigar to str
        string cigarStr = convertCigarToStr(aln.cigar);
        while (true){
            if(cigarStr.length() == 0 || cigarStr.length() < regionSize)
                return false;
            if(checkIdentityPercentange(cigarStr))
                return true;
            else {
                uint16_t lastEditLocation = max(cigarStr.find_last_of('X'), max(cigarStr.find_last_of('I'), cigarStr.find_last_of('D')));
                uint16_t firstEditLocation = min(cigarStr.find_first_of('X'), min(cigarStr.find_first_of('I'), cigarStr.find_first_of('D')));
                if (cigarStr.size() - lastEditLocation - 1 < firstEditLocation){
                    cigarStr = cigarStr.substr(0, lastEditLocation);
                    size_t insertionsBeforeEnd = 0, deletionsBeforeEnd = 0;
                    for (size_t i = 0; i < queryClips + lastEditLocation; i++) {
                        if (cigarStr[i] == 'I') {
                            insertionsBeforeEnd++;
                        } else if (cigarStr[i] == 'D') {
                            deletionsBeforeEnd++;
                        }
                    }
                    aln.queryRegionEndPos = aln.queryRegionStartPos + lastEditLocation - 1 - deletionsBeforeEnd;
                    aln.readRegionEndPos = aln.readRegionStartPos + lastEditLocation - 1 - insertionsBeforeEnd;
                    aln.cigar = convertStrToCigar(cigarStr, aln.queryRegionStartPos, aln.queryRegionEndPos);
                } else {
                    size_t insertionsBeforeStart = 0, deletionsBeforeStart = 0;
                    for (size_t i = firstEditLocation + 1; i < cigarStr.size(); i++) {
                        if (cigarStr[i] == 'I') {
                            insertionsBeforeStart++;
                        } else if (cigarStr[i] == 'D') {
                            deletionsBeforeStart++;
                        }
                    }
                    cigarStr = cigarStr.substr(firstEditLocation + 1);
                    aln.queryRegionStartPos = aln.queryRegionStartPos + firstEditLocation + 1 - deletionsBeforeStart;
                    aln.readRegionStartPos = aln.readRegionStartPos + firstEditLocation + 1 - insertionsBeforeStart;
                    aln.cigar = convertStrToCigar(cigarStr, aln.queryRegionStartPos, aln.queryRegionEndPos);
                }
            }
        }
    }

};

#endif