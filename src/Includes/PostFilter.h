#ifndef POSTFILTER_H
#define POSTFILTER_H

#include <iostream>
#include <stack>
#include "Alignment.h"

#define CHECK_REGION_SIZE 0
#define TURN_POLISHING_OFF 0

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
                // if(i > 0 && cigarStr[i-1] != '=')
                //     continue;
                return i;
            }
        }
        return -1;
    }

    int getFirstEditLocation (string cigarStr){
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'X' || cigarStr[i] == 'I' || cigarStr[i] == 'D') {
                // if((i < cigarStr.length() - 1) && cigarStr[i+1] != '=')
                //     continue;
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
        if(cigarStr.length() < k_ || (cigarStr.length() < regionSize && CHECK_REGION_SIZE)) {
            if(DEBUG_MODE) cout << "cigar len was smaller than the R or K at the beginning\n";
            aln.criteriaCode = 0x10;
            return false;
        }
        bool piCheck = checkIdentityPercentange(cigarStr);
        if (piCheck) gAlignmentFoundWithNoPolish++;
        if(TURN_POLISHING_OFF){
            if(piCheck)
                return true;
        } else {
            while (true){
                if(DEBUG_MODE) cout << "cigar string: " << cigarStr << endl;
                if((cigarStr.length() < regionSize && CHECK_REGION_SIZE) || cigarStr.length() < k_) {
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
                        if(piCheck) cout << "(< k_) query[" << aln.queryID << "] and read[" << aln.readID << "] and cigar " << cigarStr << "did not got a good second chance!!!" << endl;
                        return false;
                    }
                }
                if(checkIdentityPercentange(cigarStr)){
                    if(bCheckminNumExactMatchKmer){
                        if(!checkminNumExactMatchKmer(cigarStr)){
                            aln.criteriaCode = 0x40;
                            if(piCheck) cout << "(minMatchKmer) query[" << aln.queryID << "] and read[" << aln.readID << "] and cigar " << cigarStr << "did not got a good second chance!!!" << endl;
                            return false;
                        }
                    }
                    if(DEBUG_MODE) cout << "rightCigarSegments.size() : " << rightCigarSegments.size() << ", leftCigarSegments.size(): " << leftCigarSegments.size() << ", second chance: " << ((secondChance==true)?"T":"F") << ", firstOrLastEdit: " << ((firstOrLastEdit==true)?"T":"F") << endl;
                    if(!secondChance && ((firstOrLastEdit && rightCigarSegments.size() > 0) || (!firstOrLastEdit && leftCigarSegments.size() > 0))){
                        secondChance = true;
                        secondChanceFirstOrLastEdit = firstOrLastEdit;
                        if(foundValidAlignment){
                            if(tmpAln.partialMatchSize > aln.partialMatchSize)
                                aln = tmpAln;
                        } else 
                            aln = tmpAln;
                        foundValidAlignment = true;
                    } else if(secondChance && ((secondChanceFirstOrLastEdit && rightCigarSegments.size() > 0) || (!secondChanceFirstOrLastEdit && leftCigarSegments.size() > 0))){
                        if(foundValidAlignment){
                            if(tmpAln.partialMatchSize > aln.partialMatchSize)
                                aln = tmpAln;
                        } else 
                            aln = tmpAln;
                        foundValidAlignment = true;
                    } else {
                        if(foundValidAlignment){
                            if(tmpAln.partialMatchSize > aln.partialMatchSize)
                                aln = tmpAln;
                        } else 
                            aln = tmpAln;
                        return true;
                    }
                } else {
                    if(foundValidAlignment) {
                        return true;
                    } 
                    if (secondChance && ((secondChanceFirstOrLastEdit && rightCigarSegments.size() == 0) || (!secondChanceFirstOrLastEdit && leftCigarSegments.size() == 0))){
                        if((cigarStr.length() < regionSize && CHECK_REGION_SIZE) || cigarStr.length() < k_) {
                            aln.criteriaCode = 0x80;
                            if(piCheck) cout << "(< k_) query[" << aln.queryID << "] and read[" << aln.readID << "] and cigar " << cigarStr << "did not got a good second chance!!!" << endl;
                            return false;
                        } //else give it the last chances
                    }
                }
                if(secondChance) {
                    if(DEBUG_MODE) cout << "<<<<<<<<second chance>>>>>>>>" << endl;
                    if(DEBUG_MODE) cout << "second chance : before cigar string: " << cigarStr << endl;
                    if (secondChanceFirstOrLastEdit && rightCigarSegments.size() > 0) {
                        string rightCigarSegment = "";
                        while(rightCigarSegments.size() > 0) {
                            string top = rightCigarSegments.top();
                            if(top.size() == 1 && top[top.size()-1] != '='){
                                rightCigarSegment += top;
                                rightCigarSegments.pop();
                            } else {
                                break;
                            }
                        }
                        if (rightCigarSegments.size() > 0) {
                            rightCigarSegment += rightCigarSegments.top();
                            rightCigarSegments.pop();
                        }
                        if(DEBUG_MODE) cout << "second chance : in right and rightCigarSegment: " << rightCigarSegment << endl;
                        cigarStr = cigarStr + rightCigarSegment;
                        if(DEBUG_MODE) cout << "second chance : cigar string: " << cigarStr << endl;
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
                        string leftCigarSegment = "";
                        while(leftCigarSegments.size() > 0) {
                            string top = leftCigarSegments.top();
                            if(top.size() == 1 && top[0] != '='){
                                leftCigarSegment = top + leftCigarSegment;
                                leftCigarSegments.pop();
                            } else {
                                break;
                            }
                        }
                        if (leftCigarSegments.size() > 0) {
                            leftCigarSegment = leftCigarSegments.top() + leftCigarSegment;
                            leftCigarSegments.pop();
                        }
                        if(DEBUG_MODE) cout << "second chance : in left and leftCigarSegment: " << leftCigarSegment << endl;
                        cigarStr = leftCigarSegment + cigarStr;
                        if(DEBUG_MODE) cout << "second chance : cigar string: " << cigarStr << endl;
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
                        if(DEBUG_MODE) cout << "second chance : after: queryRegionStartPos: " << tmpAln.queryRegionStartPos << ", queryRegionEndPos: " << tmpAln.queryRegionEndPos << ", readRegionStartPos: " << tmpAln.readRegionStartPos << ", readRegionEndPos: " << tmpAln.readRegionEndPos << endl;
                        continue;
                    }
                }
                if(DEBUG_MODE) cout << "<<<<<<<<<<trimming>>>>>>>>>>" << cigarStr << endl;
                int lastEditLocation = getLastEditLocation(cigarStr);
                int firstEditLocation = getFirstEditLocation(cigarStr);
                if(DEBUG_MODE) cout << "trim : lastEditLocation: " << lastEditLocation << ", firstEditLocation: " << firstEditLocation << endl;
                if(DEBUG_MODE) cout << "trim : Before: queryRegionStartPos: " << tmpAln.queryRegionStartPos << ", queryRegionEndPos: " << tmpAln.queryRegionEndPos << ", readRegionStartPos: " << tmpAln.readRegionStartPos << ", readRegionEndPos: " << tmpAln.readRegionEndPos << endl;
                if (lastEditLocation < 0 || firstEditLocation < 0) {
                    if((cigarStr.size() > regionSize && CHECK_REGION_SIZE) || cigarStr.size() > k_) {
                        return true;
                    }
                    if(piCheck) cout << "(EditLocation) query[" << aln.queryID << "] and read[" << aln.readID << "] and cigar " << cigarStr << "did not got a good second chance!!!" << endl;
                    return false;
                }
                string rightCigarSegment = cigarStr.substr(lastEditLocation);
                string leftCigarSegment = cigarStr.substr(0, firstEditLocation + 1);
                if ((!secondChance && (countMatches(rightCigarSegment) < countMatches(leftCigarSegment))) || (secondChance && !secondChanceFirstOrLastEdit)){
                    firstOrLastEdit = false;
                    rightCigarSegments.push(rightCigarSegment);
                    cigarStr = cigarStr.substr(0, lastEditLocation);
                    if(DEBUG_MODE) cout << "trim : cigar string before pop: " << cigarStr << endl;
                    while(cigarStr.length() > 0 && cigarStr[cigarStr.length() - 1] != '=') {
                        string lastChar = "";
                        lastChar += cigarStr[cigarStr.length() - 1];
                        cigarStr = cigarStr.substr(0, cigarStr.size() - 1);
                        rightCigarSegments.push(lastChar);
                        lastEditLocation--;
                    }
                    if(DEBUG_MODE) cout << "trim : cigar string after pop: " << cigarStr << endl;
                    if(DEBUG_MODE) cout << "trim : rightCigarSegment: " << rightCigarSegment << endl;
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
                    cigarStr = cigarStr.substr(firstEditLocation + 1);
                    if(DEBUG_MODE) cout << "trim : cigar string before pop: " << cigarStr << endl;
                    while(cigarStr.length() > 0 && cigarStr[0] != '=') {
                        string firstChar = "";
                        firstChar += cigarStr[0];
                        leftCigarSegment += firstChar;
                        cigarStr = cigarStr.substr(1, cigarStr.size() - 1);
                        leftCigarSegments.push(firstChar);
                        firstEditLocation++;
                    }
                    if(DEBUG_MODE) cout << "trim : cigar string after pop: " << cigarStr << endl;
                    if(DEBUG_MODE) cout << "trim : leftCigarSegment: " << leftCigarSegment << endl;
                    size_t insertionsBeforeStart = 0, deletionsBeforeStart = 0;
                    for (size_t i = 0; i < leftCigarSegment.length(); i++) {
                        if (leftCigarSegment[i] == 'I') {
                            insertionsBeforeStart++;
                        } else if (leftCigarSegment[i] == 'D') {
                            deletionsBeforeStart++;
                        }
                    }
                    tmpAln.partialMatchSize = cigarStr.size();
                    tmpAln.queryRegionStartPos = tmpAln.queryRegionStartPos + firstEditLocation + 1 - deletionsBeforeStart;
                    tmpAln.readRegionStartPos = tmpAln.readRegionStartPos + firstEditLocation + 1 - insertionsBeforeStart;
                    tmpAln.cigar = convertStrToCigar(cigarStr, tmpAln.queryRegionStartPos, tmpAln.queryRegionEndPos);
                    util.parseCigar(tmpAln.cigar, tmpAln.matches, tmpAln.substitutions, tmpAln.inDels, tmpAln.editLocations);
                    tmpAln.partialMatchSize = tmpAln.matches + tmpAln.substitutions + tmpAln.inDels;
                    tmpAln.editDistance = tmpAln.substitutions + tmpAln.inDels;
                }
                if(DEBUG_MODE) cout << "trim : after: queryRegionStartPos: " << tmpAln.queryRegionStartPos << ", queryRegionEndPos: " << tmpAln.queryRegionEndPos << ", readRegionStartPos: " << tmpAln.readRegionStartPos << ", readRegionEndPos: " << tmpAln.readRegionEndPos << endl;
            }
        }
        if(piCheck) cout << "(end) query[" << aln.queryID << "] and read[" << aln.readID << "] and cigar " << cigarStr << "did not got a good second chance!!!" << endl;
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