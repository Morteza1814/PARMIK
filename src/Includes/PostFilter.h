#ifndef POSTFILTER_H
#define POSTFILTER_H

#include <iostream>
#include <stack>
#include "Alignment.h"

using namespace std;

class PostFilter {
private:
    uint16_t regionSize;
    uint16_t k_;
    uint16_t contigSize;
    uint16_t minNumExactMatchKmer;
    double identityPercentage;
public: 
    PostFilter(uint16_t R, uint16_t k, uint16_t c, uint16_t m, double i) : regionSize(R), k_(k), contigSize(c), minNumExactMatchKmer(m), identityPercentage(i) {}

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
                if(i > 0 && cigarStr[i-1] != '=')
                    continue;
                return i;
            }
        }
        return -1;
    }

    int getFirstEditLocation (string cigarStr){
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'X' || cigarStr[i] == 'I' || cigarStr[i] == 'D') {
                if((i < cigarStr.length() - 1) && cigarStr[i+1] != '=')
                    continue;
                return i;
            }
        }
        return -1;
    } 

    bool checkminNumExactMatchKmer(const string& cigarStr){
         uint32_t kmerCounts = 0;

        for (size_t i = 0; i <= cigarStr.length() - k_; ++i) {
            string kmer = cigarStr.substr(i, k_);
            if (kmer.find_first_not_of('=') == string::npos) {
                kmerCounts++;
            }
        }
        if(kmerCounts < minNumExactMatchKmer)
            return false;
        return true;
    }

    bool checkAndUpdateBasedOnAlingmentCriteria(Alignment &aln, bool bCheckminNumExactMatchKmer = false)
    {
        Alignment tmpAln = aln;
        bool secondChance = false, foundValidAlignment = false;
        stack<string> leftCigarSegments, rightCigarSegments;
        bool firstOrLastEdit = false, secondChanceFirstOrLastEdit = false;
        Utilities<uint32_t> util;
        // covert cigar to str
        string cigarStr = convertCigarToStr(aln.cigar);
        cigarStr = trimClips(cigarStr);
        if(DEBUG_MODE) cout << "-----------------------------------------------------\n";
        if(cigarStr.length() == 0 || cigarStr.length() < regionSize) {
            if(DEBUG_MODE) cout << "cigar len was smaller than the R at the beginning\n";
            aln.criteriaCode = 0x10;
            return false;
        }
        while (true){
            if(DEBUG_MODE) cout << "cigar string: " << cigarStr << endl;
            if(cigarStr.length() < regionSize) {
                if(foundValidAlignment) {
                    return true;
                }
                if(!secondChance && ((firstOrLastEdit && rightCigarSegments.size() > 0) || (!firstOrLastEdit && leftCigarSegments.size() > 0))) {
                    if(DEBUG_MODE) cout << "second chance" << endl;
                    //give it another chance
                    secondChance = true;
                    secondChanceFirstOrLastEdit = firstOrLastEdit;
                } else {
                    aln.criteriaCode = 0x10;
                    return false;
                }
            }
            if(checkIdentityPercentange(cigarStr)){
                if(bCheckminNumExactMatchKmer){
                    if(!checkminNumExactMatchKmer(cigarStr)){
                        aln.criteriaCode = 0x40;
                        return false;
                    }
                }
                if(DEBUG_MODE) cout << "rightCigarSegments.size() : " << rightCigarSegments.size() << ", leftCigarSegments.size(): " << leftCigarSegments.size() << ", second chance: " << ((secondChance==true)?"T":"F") << ", firstOrLastEdit: " << ((firstOrLastEdit==true)?"T":"F") << endl;
                if(!secondChance && ((firstOrLastEdit && rightCigarSegments.size() > 0) || (!firstOrLastEdit && leftCigarSegments.size() > 0))){
                    secondChance = true;
                    secondChanceFirstOrLastEdit = firstOrLastEdit;
                    foundValidAlignment = true;
                    aln = tmpAln;
                } else if(secondChance && ((secondChanceFirstOrLastEdit && rightCigarSegments.size() > 0) || (!secondChanceFirstOrLastEdit && leftCigarSegments.size() > 0))){
                    if (tmpAln.partialMatchSize > aln.partialMatchSize){
                        aln = tmpAln;
                        foundValidAlignment = true;
                    }
                } else {
                    aln = tmpAln;
                    return true;
                }
            } else {
                if(foundValidAlignment) {
                    return true;
                } 
                // if(!secondChance && ((firstOrLastEdit && rightCigarSegments.size() > 0) || (!firstOrLastEdit && leftCigarSegments.size() > 0))) {
                //     cout << "44444\n";
                //     secondChance = true;
                //     secondChanceFirstOrLastEdit = firstOrLastEdit;
                // } 
                 if (secondChance && ((secondChanceFirstOrLastEdit && rightCigarSegments.size() == 0) || (!secondChanceFirstOrLastEdit && leftCigarSegments.size() == 0))){
                    if(cigarStr.length() < regionSize) {
                        aln.criteriaCode = 0x80;
                        return false;
                    } //else give it the last chances
                }
            }
            if(secondChance) {
                if(DEBUG_MODE) cout << "second chancce" << endl;
                if(DEBUG_MODE) cout << "cigar string: " << cigarStr << endl;
                if (secondChanceFirstOrLastEdit && rightCigarSegments.size() > 0) {
                    string rightCigarSegment = rightCigarSegments.top();
                    rightCigarSegments.pop();
                    if(DEBUG_MODE) cout << "second chance an right and rightCigarSegment: " << rightCigarSegment << endl;
                    cigarStr = cigarStr + rightCigarSegment;
                    if(DEBUG_MODE) cout << "cigar string: " << cigarStr << endl;
                    size_t insertionsInRightCigarSegment = 0, deletionsInRightCigarSegment = 0;
                    for (size_t i = 0; i < rightCigarSegment.size(); i++) {
                        if (rightCigarSegment[i] == 'I') {
                            insertionsInRightCigarSegment++;
                        } else if (rightCigarSegment[i] == 'D') {
                            deletionsInRightCigarSegment++;
                        }
                    }
                    tmpAln.partialMatchSize = cigarStr.size();
                    tmpAln.queryRegionEndPos = tmpAln.queryRegionEndPos + rightCigarSegment.size() - deletionsInRightCigarSegment;
                    tmpAln.readRegionEndPos = tmpAln.readRegionEndPos + rightCigarSegment.size() - insertionsInRightCigarSegment;
                    tmpAln.cigar = convertStrToCigar(cigarStr, tmpAln.queryRegionStartPos, tmpAln.queryRegionEndPos);
                    util.parseCigar(tmpAln.cigar, tmpAln.matches, tmpAln.substitutions, tmpAln.inDels, tmpAln.editLocations);
                    tmpAln.partialMatchSize = tmpAln.matches + tmpAln.substitutions + tmpAln.inDels;
                    tmpAln.editDistance = tmpAln.substitutions + tmpAln.inDels;
                } else if (!secondChanceFirstOrLastEdit && leftCigarSegments.size() > 0) {
                    string leftCigarSegment = leftCigarSegments.top();
                    leftCigarSegments.pop();
                    if(DEBUG_MODE) cout << "second chance an left and leftCigarSegment: " << leftCigarSegment << endl;
                    cigarStr = leftCigarSegment + cigarStr;
                    if(DEBUG_MODE) cout << "cigar string: " << cigarStr << endl;
                    size_t insertionsInLeftCigarSegment = 0, deletionsInLeftCigarSegment = 0;
                    for (size_t i = 0; i < leftCigarSegment.size(); i++) {
                        if (leftCigarSegment[i] == 'I') {
                            insertionsInLeftCigarSegment++;
                        } else if (leftCigarSegment[i] == 'D') {
                            deletionsInLeftCigarSegment++;
                        }
                    }
                    tmpAln.partialMatchSize = cigarStr.size();
                    tmpAln.queryRegionStartPos = tmpAln.queryRegionStartPos - leftCigarSegment.size() + deletionsInLeftCigarSegment;
                    tmpAln.readRegionStartPos = tmpAln.readRegionStartPos - leftCigarSegment.size() + insertionsInLeftCigarSegment;
                    tmpAln.cigar = convertStrToCigar(cigarStr, tmpAln.queryRegionStartPos, tmpAln.queryRegionEndPos);
                    util.parseCigar(tmpAln.cigar, tmpAln.matches, tmpAln.substitutions, tmpAln.inDels, tmpAln.editLocations);
                    tmpAln.partialMatchSize = tmpAln.matches + tmpAln.substitutions + tmpAln.inDels;
                    tmpAln.editDistance = tmpAln.substitutions + tmpAln.inDels;
                } 
                if(checkIdentityPercentange(cigarStr)) {
                    continue;
                }
            }
            int lastEditLocation = getLastEditLocation(cigarStr);
            int firstEditLocation = getFirstEditLocation(cigarStr);
            if(DEBUG_MODE) cout << "lastEditLocation: " << lastEditLocation << ", firstEditLocation: " << firstEditLocation << endl;
            if(DEBUG_MODE) cout << "Before: queryRegionStartPos: " << tmpAln.queryRegionStartPos << ", queryRegionEndPos: " << tmpAln.queryRegionEndPos << ", readRegionStartPos: " << tmpAln.readRegionStartPos << ", readRegionEndPos: " << tmpAln.readRegionEndPos << endl;
            if (lastEditLocation < 0 || firstEditLocation < 0) {
                if(cigarStr.size() > regionSize){
                    return true;
                }
                return false;
            }
            string rightCigarSegment = cigarStr.substr(lastEditLocation);
            string leftCigarSegment = cigarStr.substr(0, firstEditLocation + 1);
            if ((!secondChance && (countMatches(rightCigarSegment) < countMatches(leftCigarSegment))) || (secondChance && !secondChanceFirstOrLastEdit)){
                firstOrLastEdit = false;
                rightCigarSegments.push(rightCigarSegment);
                if(DEBUG_MODE) cout << "rightCigarSegment: " << rightCigarSegment << endl;
                cigarStr = cigarStr.substr(0, lastEditLocation);
                size_t insertionsBeforeEnd = 0, deletionsBeforeEnd = 0;
                for (size_t i = 0; i < (uint16_t)lastEditLocation; i++) {
                    if (cigarStr[i] == 'I') {
                        insertionsBeforeEnd++;
                    } else if (cigarStr[i] == 'D') {
                        deletionsBeforeEnd++;
                    }
                }
                tmpAln.partialMatchSize = cigarStr.size();
                tmpAln.queryRegionEndPos = tmpAln.queryRegionStartPos + lastEditLocation - 1 - deletionsBeforeEnd;
                tmpAln.readRegionEndPos = tmpAln.readRegionStartPos + lastEditLocation - 1 - insertionsBeforeEnd;
                tmpAln.cigar = convertStrToCigar(cigarStr, tmpAln.queryRegionStartPos, tmpAln.queryRegionEndPos);
                util.parseCigar(tmpAln.cigar, tmpAln.matches, tmpAln.substitutions, tmpAln.inDels, tmpAln.editLocations);
                tmpAln.partialMatchSize = tmpAln.matches + tmpAln.substitutions + tmpAln.inDels;
                tmpAln.editDistance = tmpAln.substitutions + tmpAln.inDels;
            } else if ((!secondChance && (countMatches(rightCigarSegment) >= countMatches(leftCigarSegment))) || (secondChance && secondChanceFirstOrLastEdit)){
                firstOrLastEdit = true;
                leftCigarSegments.push(leftCigarSegment);
                if(DEBUG_MODE) cout << "leftCigarSegment: " << leftCigarSegment << endl;
                size_t insertionsBeforeStart = 0, deletionsBeforeStart = 0;
                for (size_t i = 0; i <= (uint16_t) firstEditLocation; i++) {
                    if (cigarStr[i] == 'I') {
                        insertionsBeforeStart++;
                    } else if (cigarStr[i] == 'D') {
                        deletionsBeforeStart++;
                    }
                }
                cigarStr = cigarStr.substr(firstEditLocation + 1);
                tmpAln.partialMatchSize = cigarStr.size();
                tmpAln.queryRegionStartPos = tmpAln.queryRegionStartPos + firstEditLocation + 1 - deletionsBeforeStart;
                tmpAln.readRegionStartPos = tmpAln.readRegionStartPos + firstEditLocation + 1 - insertionsBeforeStart;
                tmpAln.cigar = convertStrToCigar(cigarStr, tmpAln.queryRegionStartPos, tmpAln.queryRegionEndPos);
                util.parseCigar(tmpAln.cigar, tmpAln.matches, tmpAln.substitutions, tmpAln.inDels, tmpAln.editLocations);
                tmpAln.partialMatchSize = tmpAln.matches + tmpAln.substitutions + tmpAln.inDels;
                tmpAln.editDistance = tmpAln.substitutions + tmpAln.inDels;
            }
            if(DEBUG_MODE) cout << "after: queryRegionStartPos: " << tmpAln.queryRegionStartPos << ", queryRegionEndPos: " << tmpAln.queryRegionEndPos << ", readRegionStartPos: " << tmpAln.readRegionStartPos << ", readRegionEndPos: " << tmpAln.readRegionEndPos << endl;
        }
        return false;
    }

    uint16_t countMatches(string cigarStr){
        uint16_t matches = 0;
        for(char c : cigarStr){
            if(c == '=') matches++; 
        }
        return matches;
    }

};

#endif