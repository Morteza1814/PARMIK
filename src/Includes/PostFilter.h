#ifndef POSTFILTER_H
#define POSTFILTER_H

#include <iostream>

using namespace std;

class PostFilter {
private:
    size_t regionSize;
    size_t allowedEditDistance;
    size_t contigSize;
    size_t minExactMatchLength;
public: 
    PostFilter(size_t R, size_t a, size_t c, size_t m) : regionSize(R), allowedEditDistance(a), contigSize(c), minExactMatchLength(m) {}

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

    // bool hasMinConsecutiveMatches(uint32_t queryS, const string& cigarStr) {
    //     uint32_t consecutiveMatchCount = 0;

    //     // check if it is a cigarStr
    //     if (str2 == "cigarStr") {
    //         return hasMinConsecutiveMatchesWithCigarStr(cigarStr);
    //     }
    
    //     // Ensure that both strings have the same length for exact matching
    //     if (str1.length() != str2.length()) {
    //         cerr << "Error: Sequences have different lengths." << endl;
    //         return false;
    //     }

    //     // Iterate through each character in the strings and count consecutive exact matches
    //     uint32_t startPos = 0;
    //     for (size_t i = 0; i < str1.length(); i++) {
    //         if (str1[i] == str2[i]) {
    //             consecutiveMatchCount++;
    //             if (consecutiveMatchCount >= minExactMatchLength) {
    //                 if ((startPos + queryS + minExactMatchLength <= regionSize) || ((startPos + queryS >= contigSize - regionSize) && (startPos + queryS + minExactMatchLength <= contigSize))) {
    //                     cout << "11startPos: " << startPos << " queryS: " << queryS << " consecutiveMatchCount: " << consecutiveMatchCount << endl;
    //                     return true;  // Found enough consecutive matches inth front or back region
    //                 }else{
    //                     startPos++;
    //                     consecutiveMatchCount--;
    //                     cout << "22startPos: " << startPos << " queryS: " << queryS << " consecutiveMatchCount: " << consecutiveMatchCount << endl;
    //                 }
    //             }
    //         } else {
    //             consecutiveMatchCount = 0;  // Reset count if consecutive match is broken
    //             startPos = i+1;
    //         }
    //     }

    //     // Return false if the number of consecutive exact matches is less than minConsecutiveMatch
    //     return false;
    // }

    bool checkAlingmentCriteria(uint32_t editDistance, uint32_t alnLen, uint32_t queryS, string cigarStr, uint32_t mismatches, uint32_t inDels, uint32_t& criteriaCode)
    {
        /*criteriaCode
        0000 -> accepted
        0x01 -> front region dismissed because region's starting position higher than 0 (each missed bp is 1 edit distance)
        0x02 -> back region dismissed because region's starting position lower than contigSize - regionSize (each missed bp is 1 edit distance)
        0x04 -> No min exact match region in the partial match region
        0x08 -> Edit distance does not match the criteria
        0x10 -> Alignment len does not match the criteria
        0x20 -> region starts from 100 + allowedEditDistance
        */
        //check the alignment len and edit distance
        // cout << "queryS: " << queryS << " editDistance: " << editDistance << " alnLen: " << alnLen << endl;
        if((editDistance > allowedEditDistance))
        {
            criteriaCode |= 0x08;
            return false;
        }
        if(((uint32_t)alnLen < (uint32_t)regionSize))
        {
            criteriaCode |= 0x10;
            return false;
        }

        //check whether the region starts at most from 100 + 2
        if(queryS > (contigSize - regionSize + allowedEditDistance)){ //first bp of back region starts from 100
            criteriaCode |= 0x20;
            return false;
        }
        //check whether the regions has at least minExactMatchLength consecutive matches
        if(!hasMinConsecutiveMatches(cigarStr)){
            // cout << "!hasMinConsecutiveMatches" << "and queryS: " << queryS << "and queryAligned: " << queryAligned << endl;
            criteriaCode |= 0x04;
            return false;
        }

        bool frontRegionDismissed = false, backRegionDismissed = false;

        if(allowedEditDistance >= queryS){
            auto editsAllwedInFrontRegion = allowedEditDistance - queryS;
            if(mismatches + inDels > editsAllwedInFrontRegion)
                frontRegionDismissed = true;
        } else {
            frontRegionDismissed = true;
        }
        if (frontRegionDismissed)
            criteriaCode |= 0x01;

        if(queryS + alnLen >= contigSize - allowedEditDistance){ //query start position + alignment len should be larger than conig size - allowed edit
            auto editsAllwedInBackRegion = queryS + alnLen - (contigSize - allowedEditDistance);
            if(mismatches + inDels > editsAllwedInBackRegion)
                backRegionDismissed = true;
        }else {
            backRegionDismissed = true;
        }
        if (backRegionDismissed)
            criteriaCode |= 0x02;

        if(frontRegionDismissed && backRegionDismissed)
            return false;
        return true;
    }

    //  void testCheckBlastEditPositions(uint32_t contigSize, uint32_t regionSize, uint32_t editDistance, uint32_t minExactMatchLength,
    // uint32_t queryS, string qAln, string readAln, uint32_t alnLen, uint32_t blastMismatches, uint32_t blastInDel){
    //     BlastReader::Blast blastAlignment;
    //     Config cfg;
    //     contigSize = contigSize;
    //     regionSize = regionSize;
    //     editDistance = editDistance;
    //     minExactMatchLength = minExactMatchLength;
    //     blastAlignment.queryAligned = qAln;
    //     blastAlignment.readAligned = readAln;
    //     mismatches = blastMismatches;
    //     inDels = blastInDel;
    //     blastAlignment.queryS = queryS;
    //     blastAlignment.AlignmentLength = alnLen;
    //     if(checkBlastEditPositions(blastAlignment, cfg))
    //         cout << "blast alignment fits to the criteria" << endl;
    //     else
    //         cout << "blast alignment does not fit to the criteria" << endl;
    // }

};

#endif