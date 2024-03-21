#ifndef POSTFILTER_H
#define POSTFILTER_H

#include <iostream>
#include "Alignment.h"

using namespace std;

class PostFilter {
private:
    uint32_t regionSize;
    uint32_t allowedEditDistance;
    uint32_t contigSize;
    double identityPercentage;
public: 
    PostFilter(uint32_t R, uint32_t a, uint32_t c, double i) : regionSize(R), allowedEditDistance(a), contigSize(c), identityPercentage(i) {}

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
        double identity =  (double) matches / (double) cigarStr.size();
        if(DEBUG_MODE) cout << "matches: " << matches << ", cigarStr.size : " << cigarStr.size() << ", identity: " << identity << endl;
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
        if (queryS > 0) {
            cigar << queryS << 'S';
        }
        for (size_t i = 1; i < cigarStr.length(); ++i) {
            char cur = cigarStr[i];
            num++;
            if (prev != cur)
            {
                cigar << num << prev;
                num = 0;
                prev = cur;
            }
        }
        num++;
        cigar << num << prev;
        int remainedClips = contigSize - queryE - 1;
        if (remainedClips > 0) {
            cigar << remainedClips << 'S';
        }
        return cigar.str();
    }

    uint32_t getClipCount(string cigarStr) {
        uint32_t cnt = 0;
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'S' || cigarStr[i] == 'H') {
                cnt++;
            }
        }
        return cnt;
    }

    string trimClips(string cigarStr){
        string trimmedCigarStr = "";
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'S' || cigarStr[i] == 'H') {
                continue;
            }
            trimmedCigarStr += cigarStr[i];
        }
        return trimmedCigarStr;
    }

    int getLastEditLocation (string cigarStr){
        for (size_t i = cigarStr.length()-1; i >= 0; i--) {
            if (cigarStr[i] == 'X' || cigarStr[i] == 'I' || cigarStr[i] == 'D') {
                return i;
            }
        }
        return -1;
    }

    int getFirstEditLocation (string cigarStr){
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'X' || cigarStr[i] == 'I' || cigarStr[i] == 'D') {
                return i;
            }
        }
        return -1;
    }

    bool checkAndUpdateBasedOnAlingmentCriteria(Alignment &aln)
    {
        Utilities<uint32_t> util;
        // covert cigar to str
        string cigarStr = convertCigarToStr(aln.cigar);
        cigarStr = trimClips(cigarStr);
        if(DEBUG_MODE) cout << "-----------------------------------------------------\n";
        while (true){
            if(DEBUG_MODE) cout << "cigar string: " << cigarStr << endl;
            if(cigarStr.length() == 0 || cigarStr.length() < regionSize)
                return false;
            if(checkIdentityPercentange(cigarStr))
                return true;
            else {
                int lastEditLocation = getLastEditLocation(cigarStr);
                int firstEditLocation = getFirstEditLocation(cigarStr);
                if(DEBUG_MODE) cout << "lastEditLocation: " << lastEditLocation << ", firstEditLocation: " << firstEditLocation << endl;
                if (lastEditLocation < 0 || firstEditLocation < 0) {
                    return false;
                }
                if (cigarStr.size() - (uint16_t)lastEditLocation - 1 < (uint16_t)firstEditLocation){
                    cigarStr = cigarStr.substr(0, lastEditLocation);
                    size_t insertionsBeforeEnd = 0, deletionsBeforeEnd = 0;
                    for (size_t i = 0; i < (uint16_t)lastEditLocation; i++) {
                        if (cigarStr[i] == 'I') {
                            insertionsBeforeEnd++;
                        } else if (cigarStr[i] == 'D') {
                            deletionsBeforeEnd++;
                        }
                    }
                    aln.partialMatchSize = cigarStr.size();
                    aln.queryRegionEndPos = aln.queryRegionStartPos + lastEditLocation - 1 - deletionsBeforeEnd;
                    aln.readRegionEndPos = aln.readRegionStartPos + lastEditLocation - 1 - insertionsBeforeEnd;
                    aln.cigar = convertStrToCigar(cigarStr, aln.queryRegionStartPos, aln.queryRegionEndPos);
                    util.parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, aln.editLocations);
                    aln.partialMatchSize = aln.matches + aln.substitutions + aln.inDels;
                    aln.editDistance = aln.substitutions + aln.inDels;
                } else {
                    size_t insertionsBeforeStart = 0, deletionsBeforeStart = 0;
                    for (size_t i = (uint16_t)firstEditLocation + 1; i < cigarStr.size(); i++) {
                        if (cigarStr[i] == 'I') {
                            insertionsBeforeStart++;
                        } else if (cigarStr[i] == 'D') {
                            deletionsBeforeStart++;
                        }
                    }
                    cigarStr = cigarStr.substr(firstEditLocation + 1);
                    aln.partialMatchSize = cigarStr.size();
                    aln.queryRegionStartPos = aln.queryRegionStartPos + firstEditLocation + 1 - deletionsBeforeStart;
                    aln.readRegionStartPos = aln.readRegionStartPos + firstEditLocation + 1 - insertionsBeforeStart;
                    aln.cigar = convertStrToCigar(cigarStr, aln.queryRegionStartPos, aln.queryRegionEndPos);
                    util.parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, aln.editLocations);
                    aln.partialMatchSize = aln.matches + aln.substitutions + aln.inDels;
                    aln.editDistance = aln.substitutions + aln.inDels;
                }
            }
        }
    }

};

#endif