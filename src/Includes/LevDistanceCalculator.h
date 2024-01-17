#ifndef LEVDISTANCECALCULATOR_H
#define LEVDISTANCECALCULATOR_H

#include "Configs.h"
#include "SamReader.h"
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))

typedef struct LevenshteinDistanceAlignment
{
    int readID = -1;            //read ID of the alignment
    int queryID = -1;            //query ID of the alignment
    string read;                //the region in read evaluated for the partial match
    string query;               //the region in query evaluated for the partial match
    string alignedRead;         //the region in read aligned for the partial match
    string alignedQuery;        //the region in query aligned for the partial match
    string editDistanceTypes;   //edit distance types of the region of the partial match
    unsigned int editDistance = 0;           //number of edit distances detected in the region
    int partialMatchSize = 0;       //partial match region size
    vector<unsigned int> editPositions;  //the positions of the edit distances in the partial match region
    string cigar;               //CIGAR string of the alignment
    uint32_t numberOfSub = 0;       //number of substitutions
    uint32_t numberOfInDel = 0;       //number of InDel
    uint32_t numberOfMatches = 0;       //number of matched bp
    uint16_t readRegionStartPos = 0;
    uint16_t readRegionEndPos = 0;
    uint16_t queryRegionStartPos = 0;
    uint16_t queryRegionEndPos = 0;
    uint16_t flag = 0;                  // determines the strand for now
}LevAlign;

class LevDistanceCalculator {
public:
    vector<unsigned int> getEditDistancePositions(string editDistanceType)
    {
        vector<unsigned int> editPositions;
        // if (cfg.isVerboseLog) cout << "edtype = " << editDistanceType <<endl;
        for (size_t i = 0; i < editDistanceType.size(); i++)
        {
            if (editDistanceType[i] != '-')
            {
                editPositions.push_back(i);
            }
        }
        // for (size_t i = 0; i < editPositions.size(); i++)
        // {
        // 	if (cfg.isVerboseLog) cout << "edit [" << i << "] pos = " << editPositions[i] << endl;
        // }cout << endl;
        return editPositions;
    }

    unsigned int getEditDistance(string editDistanceType)
    {
        uint32_t ed = 0;
        // if (cfg.isVerboseLog) cout << "edtype = " << editDistanceType <<endl;
        for (size_t i = 0; i < editDistanceType.size(); i++)
        {
            if (editDistanceType[i] != '-')
            {
                ed++;
            }
        }
        return ed;
    }
    void printMatrix(const vector<vector<double>>& matrix) {
        for (const auto& row : matrix) {
            for (const auto& element : row) {
                cout << element << "\t";
            }
            cout << endl;
        }
}

    unsigned int edit_distance(LevAlign *la, double subPenalty, double inDelPenalty)
    {
        const int len1 = la->read.size(), len2 = la->query.size();
        vector<vector<double>> d(len1 + 1, vector<double>(len2 + 1));
        // cout << "indel penalty: " << inDelPenalty << endl;
        d[0][0] = 0;
        for (double i = 1; i <= len1; ++i)
            d[i][0] = i * inDelPenalty;
        for (double i = 1; i <= len2; ++i)
            d[0][i] = i * inDelPenalty;

        // if (len1 != len2)
        // {
        //     cout << "\nQ and R MEM size is not equal!" << endl;
        // }
        // make the alignment 'N' insensitive
        string query = la->query, read = la->read;
        for (int i = 0; i < len1; i++)
        {
            if (query[i] != read[i] && query[i] == 'N')
            {
                read[i] = 'N';
            }
            else if (query[i] != read[i] && read[i] == 'N')
            {
                query[i] = 'N';
            }
        }
        // calculate the edit distance
        for (int i = 1; i <= len1; ++i)
            for (int j = 1; j <= len2; ++j)
        {
                d[i][j] = MIN3(d[i - 1][j] + inDelPenalty, d[i][j - 1] + inDelPenalty, d[i - 1][j - 1] + (read[i - 1] == query[j - 1] ? 0 : subPenalty));
        }
        // printMatrix(d);
        int i = len1, j = len2;
        // if (cfg.isVerboseLog) cout << "\nedit distance : " << d[len1][len2] << endl;
        string alignedRead = la->read, alignedQuery = la->query;
        string editDistanceTypes;
        // back tracking to find the alignment and ed types
        while (i > 0 && j > 0)
        {
            if (d[i - 1][j - 1] <= d[i][j - 1] && d[i - 1][j - 1] <= d[i - 1][j])
            {
                if (d[i - 1][j - 1] != d[i][j])
                    editDistanceTypes.insert(0, 1, 's');
                else
                    editDistanceTypes.insert(0, 1, '-');
                i--;
                j--;
            }
            else if (d[i - 1][j] <= d[i - 1][j - 1] && d[i - 1][j] <= d[i][j - 1])
            {
                editDistanceTypes.insert(0, 1, 'd');
                alignedQuery.insert(j, 1, '-');
                i--;
            }
            else
            {
                editDistanceTypes.insert(0, 1, 'i');
                alignedRead.insert(i, 1, '+');
                j--;
            }
        }
        // indel at the beginning of the strings if there are more bp in one of the sequences
        alignedRead.insert(0, j, '+');
        editDistanceTypes.insert(0, j, 'i');
        alignedQuery.insert(0, i, '-');
        editDistanceTypes.insert(0, i, 'd');
        la->alignedRead = alignedRead;
        la->alignedQuery = alignedQuery;
        la->editDistanceTypes = editDistanceTypes;
        la->editDistance = getEditDistance(editDistanceTypes);
        //remove the indDels at the end
        while(la->editDistanceTypes[la->alignedQuery.size()-1] == 'i' || la->editDistanceTypes[la->alignedQuery.size()-1] == 'd')
        {
            la->editDistanceTypes = la->editDistanceTypes.substr(0, la->editDistanceTypes.size()-1);
            la->alignedQuery = la->alignedQuery.substr(0, la->alignedQuery.size()-1);
            la->alignedRead = la->alignedRead.substr(0, la->alignedRead.size()-1);
            la->editDistance--;
        }

        LevDistanceCalculator ldc;
        SamReader s("");
        uint32_t numberOfSub = 0, numberOfInDel = 0;
        for (size_t j = 0; j < la->editDistanceTypes.size(); j++)
        {
            if (la->editDistanceTypes[j] == 's')
            {
                numberOfSub++;
            } else if (la->editDistanceTypes[j] == 'i' || la->editDistanceTypes[j] == 'd')
            {
                numberOfInDel++;
            }
        }
        la->numberOfSub = numberOfSub;
        la->numberOfInDel = numberOfInDel;
        la->cigar = ldc.alignmentToCIGAR(la->editDistanceTypes);
        la->numberOfMatches = s.countMatches(la->cigar);
        return la->editDistance;
    }

    int getInDelCount (string editDistanceType)
    {
        int indelCount = 0;
        for (size_t i = 0; i < editDistanceType.size(); i++)
        {
            if (editDistanceType[i] == 'i' || editDistanceType[i] == 'd')
            {
                indelCount++;
            }
        }
        return indelCount;
    }
    
    string alignmentToCIGAR(const string& alignmentResult) {
        string cigar;
        string alignment;
        for (size_t i = 0; i < alignmentResult.size(); ++i) {
            if (alignmentResult[i] == 's') 
                alignment += '-';
            else
                alignment += alignmentResult[i];
        }
        char currentOp = alignment[0];
        int opCount = 1;
        char op = 'M';
        for (size_t i = 1; i < alignment.size(); ++i) {
            if (alignment[i] == currentOp) {
                opCount++;
            } else {
                if (currentOp == '-')
                {
                    op = 'M';
                } else if (currentOp == 'i')
                {
                    op = 'I';
                } else if (currentOp == 'd')
                {
                    op = 'D';
                }
                cigar += to_string(opCount) + op;
                currentOp = alignment[i];
                opCount = 1;
            }
        }
        if (currentOp == '-')
        {
            op = 'M';
        } else if (currentOp == 'i')
        {
            op = 'I';
        } else if (currentOp == 'd')
        {
            op = 'D';
        }
        cigar += to_string(opCount) + op; // Add the last set of operations
        return cigar;
    }
};

#endif