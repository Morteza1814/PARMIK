#ifndef SEEDMATCHEXTENDER_H
#define SEEDMATCHEXTENDER_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "NucleotideEncoder.h"
#include "IndexContainer.h"
#include "LevDistanceCalculator.h"
#include "Utils.h"
#include "SamReader.h"
#include "tsl/robin_map.h"

using namespace tsl;
using namespace std;

template <typename contigIndT, typename kmerT>
class SeedMatchExtender
{
private:
    NucleotideEncoder<kmerT> encoder_;
    size_t k_;
    char frontRegion = 'F';
    char backRegion = 'B';
    size_t regionSize;
    bool isVerboseLog_ = false;
public:
    SeedMatchExtender(size_t k, size_t R, bool isVerboseLog) : encoder_(k), k_(k), regionSize(R), isVerboseLog_(isVerboseLog) {}

    string readOneContig(const string &filename, contigIndT n)
    {
        Utilities<contigIndT> utl;
        ifstream file(filename);
        if (!file.is_open())
        {
            cerr << "Error opening file: " << filename << endl;
            return "";
        }

        string line;
        for (contigIndT i = 0; i <= n; i++)
        {
            if (!getline(file, line))
            {
                cerr << "Error reading line " << n + 1 << " from file: " << filename << endl;
                return "";
            }
            if (i == n)
            {
                contigIndT id = utl.extractContigId(line);
                if (id != n)
                {
                    cerr << "Error contig id does not match " << to_string(id) << " with n: " << to_string(n) << endl;
                    return "";
                }
            }
            if (!getline(file, line))
            {
                cerr << "Error reading line " << n + 1 << " from file: " << filename << endl;
                return "";
            }
        }

        file.close();
        return line;
    }

    map<contigIndT, string> readContigs(const string &filename, set<contigIndT> contigIdSet)
    {
        ifstream file(filename);
        Utilities<contigIndT> utl;
        map<contigIndT, string> contigs;
        if (!file.is_open())
        {
            cerr << "Error opening file: " << filename << endl;
            return contigs;
        }
        unsigned int setInd = 0, i = 0;
        string line;
        while (getline(file, line))
        {
            contigIndT setElement;
            if (setInd >= 0 && setInd < contigIdSet.size())
            {
                auto it = contigIdSet.begin();
                advance(it, setInd);
                setElement = *it;
            }
            if ((setElement * 2) == i)
            {
                contigIndT id = utl.extractContigId(line);
                if (!getline(file, line))
                {
                    cerr << "Error reading contig: " << setElement << " from file: " << filename << endl;
                    return contigs;
                }
                i++;
                contigs[id] = line;
                setInd++;
                if (setInd == contigIdSet.size() - 1)
                {
                    break;
                }
            }
            i++;
        }
        file.close();
        return contigs;
    }

    vector<kmerT> extractReadKmers(const string &sequence)
    {
        vector<kmerT> kmers;
        if (sequence.size() < k_)
        {
            return kmers;
        }

        for (size_t i = 0; i <= sequence.size() - k_; i += 1)
        {
            string kmerString = sequence.substr(i, k_);
            kmerT encodedKmer = encoder_.encode(kmerString);
            kmers.push_back(encodedKmer);
        }

        return kmers;
    }

    void extractKmersFromRegions(const string &sequence, vector<kmerT> &regionKmers, char region)
    {
        if (sequence.size() < k_)
        {
            cout << "sequence : " << sequence << endl;
            throw invalid_argument("less than k sequence size");
        }
        if (isVerboseLog_) cout << "<<<<<<<<<<<<<<<<<<<<region>>>>>>>>>>>>>>>>>>>> : " << region << endl;
        if (region == frontRegion)
        {
            for (size_t i = 0; i <= regionSize - k_; i++)
            {
                string kmerString = sequence.substr(i, k_);
                kmerT encodedKmer = encoder_.encode(kmerString);
                regionKmers.push_back(encodedKmer);
                // cout << "kmer = " << kmerString << endl;
            }
        }
        else if (region == backRegion)
        {
            for (size_t i = (sequence.size() - regionSize); i <= sequence.size() - k_; i++)
            {
                string kmerString = sequence.substr(i, k_);
                // cout << "kmer = " << kmerString << endl;
                kmerT encodedKmer = encoder_.encode(kmerString);
                regionKmers.push_back(encodedKmer);
            }
        }
        else
        {
            cerr << "wrong region char!!" << endl;
        }
    }

    vector<pair<unsigned int, unsigned int>> findCommonKmers(const vector<kmerT> &qKmers, const vector<kmerT> &rKmers, char region, unsigned int contigSize)
    {
        vector<pair<unsigned int, unsigned int>> commonPositions;
        for (unsigned int i = 0; i < qKmers.size(); i++)
        {
            for (unsigned int j = 0; j < rKmers.size(); j++)
            {
                if (qKmers[i] == rKmers[j])
                {
                    // cout << "i = " << i << ", j = " << j << endl;
                    unsigned int regionOffset = (region == backRegion) ? (contigSize - regionSize) : 0;
                    // cout << "region : " << region << ", regionOffset = " << regionOffset << endl;
                    commonPositions.push_back(make_pair(i + regionOffset, j));
                }
            }
        }
        return commonPositions;
    }

    LevAlign extendPreMem(int readMemStartInd, int queryMemStartInd, string read, string query, uint32_t allowedEditDistance)
    {
        LevDistanceCalculator ldc;
        LevAlign preLa;
        unsigned int prefixRegionSize = 0;
        if (readMemStartInd == 0 || queryMemStartInd == 0) // no Prefix region check is needed
        {
            prefixRegionSize = 0;
        }
        else if (readMemStartInd > queryMemStartInd) // read Prefix region is larger
        {
            prefixRegionSize = queryMemStartInd;
        }
        else // query Prefix region is larger
        {
            prefixRegionSize = readMemStartInd;
        }
        if (prefixRegionSize > 0)
        {
            preLa.readRegionStartPos = readMemStartInd - prefixRegionSize;
            preLa.queryRegionStartPos = queryMemStartInd - prefixRegionSize;
            preLa.read = read.substr(preLa.readRegionStartPos, prefixRegionSize);
            preLa.query = query.substr(preLa.queryRegionStartPos, prefixRegionSize);
            unsigned int preEd = ldc.edit_distance(&preLa);
            vector<unsigned int> editPos= ldc.getEditDistancePositions(preLa.editDistanceTypes);
            if (preEd <= allowedEditDistance) // accept the whole prefix
            {
                preLa.partialMatchSize = preLa.editDistanceTypes.size();
                for (size_t i = 0; i < editPos.size(); i++)
                {
                    preLa.editPositions.push_back(editPos[i]);
                }
            }
            else
            {
                unsigned int regionStartPos = editPos[editPos.size() - 1 - allowedEditDistance] + 1;
                for (size_t i = editPos.size()- allowedEditDistance; i < editPos.size(); i++)
                {
                    preLa.editPositions.push_back(editPos[i] - regionStartPos);
                }
                preLa.partialMatchSize = preLa.editDistanceTypes.substr(regionStartPos).size();
                // cout << "pre regionStartPos: " << regionStartPos << endl;
                preLa.editDistanceTypes = preLa.editDistanceTypes.substr(regionStartPos);
                preLa.partialMatchSize = preLa.editDistanceTypes.size();
                preLa.alignedQuery = preLa.alignedQuery.substr(regionStartPos);
                preLa.alignedRead = preLa.alignedRead.substr(regionStartPos);
                preLa.editDistance = allowedEditDistance;
                // cout << "pre partialMatchSize: " << preLa.partialMatchSize << endl;
            }
        }
        return preLa;
    }

    LevAlign extendSufMem(unsigned int readMemEndInd, unsigned int queryMemEndInd, string read, string query, uint32_t allowedEditDistance, uint32_t contigSize)
    {
        LevDistanceCalculator ldc;
        LevAlign sufLa;
        unsigned int suffixRegionSize = 0;
        if (readMemEndInd == (contigSize - 1) || queryMemEndInd == (contigSize - 1)) // no Suffix region check is needed
        {
            suffixRegionSize = 0;
        }
        else if (readMemEndInd > queryMemEndInd) // query Suffix region is larger
        {
            suffixRegionSize = contigSize - (readMemEndInd + 1);
        }
        else // query Suffix region is larger
        {
            suffixRegionSize = contigSize - (queryMemEndInd + 1);
        }
        if (suffixRegionSize > 0)
        {
            sufLa.readRegionEndPos = readMemEndInd + suffixRegionSize;
            sufLa.queryRegionEndPos = queryMemEndInd + suffixRegionSize;
            sufLa.read = read.substr(readMemEndInd + 1, suffixRegionSize);
            sufLa.query = query.substr(queryMemEndInd + 1, suffixRegionSize);
            unsigned int sufEd = ldc.edit_distance(&sufLa);
            vector<unsigned int> editPos = ldc.getEditDistancePositions(sufLa.editDistanceTypes);
            // cout << "suf query: " << sufLa.query << endl;
            // cout << "suf ed: " << sufEd << endl;
            // cout << "sufLa.query : " <<  sufLa.query << endl;
            // cout <<  "sufLa.alignedQuery : " << sufLa.alignedQuery << endl;
            // cout << "sufLa.read : " <<  sufLa.read << endl;
            // cout <<  "sufLa.alignedRead : " << sufLa.alignedRead << endl;
            if (sufEd <= allowedEditDistance) // accept the whole prefix
            {
                sufLa.partialMatchSize = sufLa.editDistanceTypes.size();
                for (size_t i = 0; i < editPos.size(); i++)
                {
                    sufLa.editPositions.push_back(editPos[i]);
                }
            }
            else
            {
                unsigned int regionEndPos = editPos[allowedEditDistance] - 1;
                // cout << "suf regionEndPos: " << regionEndPos << endl;
                for (size_t i = 0; i < allowedEditDistance; i++)
                {
                    sufLa.editPositions.push_back(editPos[i]);
                }
                sufLa.editDistanceTypes = sufLa.editDistanceTypes.substr(0, regionEndPos + 1);
                sufLa.partialMatchSize = sufLa.editDistanceTypes.size();
                sufLa.alignedQuery = sufLa.alignedQuery.substr(0, regionEndPos + 1);
                sufLa.alignedRead = sufLa.alignedRead.substr(0, regionEndPos + 1);
                sufLa.editDistance = allowedEditDistance;
                // cout << "suf partialMatchSize: " << sufLa.partialMatchSize << endl;
            }
        }
        return sufLa;
    }

    string extendMEM(unsigned int readSeedStartPos, unsigned int querySeedStartPos, string read, string query, uint32_t &readMemStartInd, uint32_t &queryMemStartInd, uint32_t &readMemEndInd, uint32_t &queryMemEndInd, unsigned int contigSize)
    {
        readMemStartInd = readSeedStartPos;
        queryMemStartInd = querySeedStartPos;
        readMemEndInd = (readSeedStartPos + k_ - 1);
        queryMemEndInd = (querySeedStartPos + k_ - 1);
        // extend mem to left
        for (int i = readSeedStartPos - 1, j = querySeedStartPos - 1; i >= 0 && j >= 0; i--, j--)
        {
            if (read[i] == query[j])
            {
                // maxExactMatchLen++;
                readMemStartInd = i;
                queryMemStartInd = j;
            }
            else
            {
                break;
            }
        }
        // extend mem to right
        for (unsigned int i = (readSeedStartPos + k_), j = (querySeedStartPos + k_); i < contigSize && j < contigSize; i++, j++)
        {
            if (read[i] == query[j])
            {
                readMemEndInd = i;
                queryMemEndInd = j;
            }
            else
            {
                break;
            }
        }
        return read.substr(readMemStartInd, (readMemEndInd - readMemStartInd + 1));
    }

    LevAlign findMaximalPartialMatch(uint32_t allowedEditDistance, uint32_t regionSize, LevAlign preLa, LevAlign sufLa, LevAlign memLa)
    {
        uint32_t maxPM = 0;
        LevAlign midLa;
        // cout << "regionSize : " << regionSize << endl;
        // whole region is accepted
        if (preLa.editPositions.size() + sufLa.editPositions.size() <= allowedEditDistance)
        {
            // cout << "edit is nice : " << preLa.editPositions.size() << " , " << sufLa.editPositions.size() << endl;
            // cout << sufLa.editDistanceTypes << endl;
            // cout << sufLa.alignedQuery << endl;
            midLa.editDistance = preLa.editPositions.size() + sufLa.editPositions.size();
            midLa.editPositions.insert(midLa.editPositions.end(), preLa.editPositions.begin(), preLa.editPositions.end());
            midLa.editPositions.insert(midLa.editPositions.end(), sufLa.editPositions.begin(), sufLa.editPositions.end());
            midLa.partialMatchSize = regionSize;
            midLa.alignedQuery = preLa.alignedQuery.substr(0, preLa.alignedQuery.size() - memLa.alignedQuery.size()) + memLa.alignedQuery + sufLa.alignedQuery.substr(memLa.alignedQuery.size());
            midLa.alignedRead = preLa.alignedRead.substr(0, preLa.alignedRead.size() - memLa.alignedRead.size()) + memLa.alignedRead + sufLa.alignedRead.substr(memLa.alignedRead.size());
            midLa.editDistanceTypes = preLa.editDistanceTypes.substr(0, preLa.editDistanceTypes.size() - memLa.editDistanceTypes.size()) + memLa.editDistanceTypes + sufLa.editDistanceTypes.substr(memLa.editDistanceTypes.size());
            return midLa;
        }
        else
        {
            int start = 0, i = 0;
            uint32_t ed = 0, j = 0, end = 0;
            i = preLa.editPositions.size() - 1;
            // cout << "i 0: " << i << endl;
            bool startFromBeginning = false, finishAtEnd = false;
            start = preLa.editPositions[i]; // begining of mem
            while (i >= 0 && ed <= allowedEditDistance)
            {
                // cout << "i 1.5: " << i << ", ed : " << ed << endl;
                i--;
                ed++;
            }
            // cout << "i 1: " << i << endl;
            if (ed >= allowedEditDistance)
            {
                ed = allowedEditDistance;
                if (i < 0)
                {
                    // start = (preLa.editPositions[0] == 0) ? 1 : 0;
                    start = 0;
                    i = 0;
                    // if (preLa.editPositions[0] > 0)
                    // {
                        startFromBeginning = true;
                    // }
                }
                else
                {
                    start = preLa.editPositions[i] + 1;
                }
                end = sufLa.editPositions[0]; // end is at the end of MEM
                // cout << "end -1 => " << end << endl;
            }
            else
            {
                if (i < 0)
                {
                    i = 0;
                }
                // start = (preLa.editPositions[0] == 0) ? 1 : 0;
                start = 0;
                // cout << "startak : " << start << ", ed : " << ed << endl;
                // if (preLa.editPositions[0] > 0)
                // {
                    startFromBeginning = true;
                // }
            }
            // continue in suffix region            
            if (ed < allowedEditDistance)
            {
                while (j < sufLa.editPositions.size() && ed < allowedEditDistance)
                {
                    j++;
                    ed++;
                }
                if (ed >= allowedEditDistance)
                {
                    ed = allowedEditDistance;
                    if (j >= sufLa.editPositions.size())
                    {
                        // end = (sufLa.editPositions[sufLa.editPositions.size() - 1] == (sufLa.partialMatchSize - 1)) ? sufLa.partialMatchSize - 1 : sufLa.partialMatchSize;
                        end = sufLa.partialMatchSize - memLa.partialMatchSize;
                        // cout << "end 0 => " << end << endl;
                        j = sufLa.editPositions.size() - 1;
                        // if (sufLa.editPositions[sufLa.editPositions.size() - 1] < (uint32_t)(sufLa.partialMatchSize - 1))
                        // {
                            finishAtEnd = true;
                        // }
                    }
                    else
                    {
                        end = sufLa.editPositions[j];
                        // cout << "end 1 => " << end << endl;
                    }
                }
                else
                {
                    // end = (sufLa.editPositions[sufLa.editPositions.size() - 1] == (sufLa.partialMatchSize - 1)) ? sufLa.partialMatchSize - 1 : sufLa.partialMatchSize;
                    end = sufLa.partialMatchSize - memLa.partialMatchSize;
                    // cout << "end 2 => " << end << endl;
                    j = sufLa.editPositions.size() - 1;
                    // if (sufLa.editPositions[sufLa.editPositions.size() - 1] < (uint32_t)(sufLa.partialMatchSize - 1))
                    // {
                        finishAtEnd = true;
                    // }
                }
            }
            // cout << "preLa.partialMatchSiz : " << preLa.partialMatchSize << endl;
            maxPM = (preLa.partialMatchSize - start) + (end);
            // cout << "end : " << end << ", start = " << start << ", ed = " << ed << endl;
            // cout << "i = " << i << ", j = " << j << endl;
            // cout << "maxPM = " << maxPM << endl;
            midLa.editDistance = ed;
            for (size_t m = i; m < preLa.editPositions.size(); m++)
            {
                midLa.editPositions.push_back(preLa.editPositions[m]);
            }
            for (size_t m = i; m < sufLa.editPositions.size(); m++)
            {
                midLa.editPositions.push_back(sufLa.editPositions[m]);
            }
            midLa.partialMatchSize = maxPM;
            midLa.alignedQuery = preLa.alignedQuery.substr(start, preLa.alignedQuery.size() - memLa.alignedQuery.size() - start) + memLa.alignedQuery + sufLa.alignedQuery.substr(memLa.alignedQuery.size(), end);
            midLa.alignedRead = preLa.alignedRead.substr(start, preLa.alignedRead.size() - memLa.alignedRead.size() - start) + memLa.alignedRead + sufLa.alignedRead.substr(memLa.alignedRead.size(), end);
            midLa.editDistanceTypes = preLa.editDistanceTypes.substr(start, preLa.editDistanceTypes.size() - memLa.editDistanceTypes.size() - start) + memLa.editDistanceTypes + sufLa.editDistanceTypes.substr(memLa.editDistanceTypes.size(), end);
            if (!startFromBeginning) i++;
            if (!finishAtEnd) j++;
            while ((uint32_t) i < preLa.editPositions.size() && j <= sufLa.editPositions.size())
            {
                // cout << "i : " << i << ", j : " << j << endl;
                start = preLa.editPositions[i] + 1;
                bool timeToBreak = false;
                if (j >= sufLa.editPositions.size())
                {
                    // end = (sufLa.editPositions[sufLa.editPositions.size() - 1] == (sufLa.partialMatchSize - 1)) ? sufLa.partialMatchSize - 1 : sufLa.partialMatchSize;
                    end = sufLa.partialMatchSize - memLa.partialMatchSize;
                    j = sufLa.editPositions.size() - 1;
                    timeToBreak = true;
                }
                else
                {
                    end = sufLa.editPositions[j];
                }
                uint32_t maxPM_t = (preLa.partialMatchSize - start) + (end);
                // cout << "maxPM_t = " << maxPM_t << endl;
                // cout << "preLa.partialMatchSize  = " << preLa.partialMatchSize << ", start " << start << ", end = " << end << ",  memLa.partialMatchSize = " <<  memLa.partialMatchSize << endl;
                if (maxPM_t > maxPM)
                {
                    // cout << "maxPM_t >  maxPM" << endl;
                    maxPM = maxPM_t;
                    midLa.editDistance = ed;
                    for (size_t m = i; m < preLa.editPositions.size(); m++)
                    {
                        midLa.editPositions.push_back(preLa.editPositions[m]);
                    }
                    for (size_t m = 0 ; m <= j; m++)
                    {
                        midLa.editPositions.push_back(sufLa.editPositions[m]);
                    }
                    midLa.partialMatchSize = maxPM;
                    midLa.alignedQuery = preLa.alignedQuery.substr(start, preLa.alignedQuery.size() - memLa.alignedQuery.size() - start) + memLa.alignedQuery + sufLa.alignedQuery.substr(memLa.alignedQuery.size(), end);
                    midLa.alignedRead = preLa.alignedRead.substr(start, preLa.alignedRead.size() - memLa.alignedRead.size() - start) + memLa.alignedRead + sufLa.alignedRead.substr(memLa.alignedRead.size(), end);
                    midLa.editDistanceTypes = preLa.editDistanceTypes.substr(start, preLa.editDistanceTypes.size() - memLa.editDistanceTypes.size() - start) + memLa.editDistanceTypes + sufLa.editDistanceTypes.substr(memLa.editDistanceTypes.size(), end);
                    // cout << "end = " << end << ", start = " << start << ", ed = " << ed << endl;
                }
                if (timeToBreak)
                {
                    break;
                }
                
                i++;
                j++;
            }
        }
        return midLa;
    }

    LevAlign comparePartialMatchRes(vector<LevAlign>& pms)
    {
        LevAlign bestAln;
        LevDistanceCalculator ldc;
        SamReader s("");
        for (size_t i = 0; i < pms.size(); i++)
        {
            int numberOfSub = 0, numberOfInDel = 0;
            for (size_t j = 0; j < pms[i].editDistanceTypes.size(); j++)
            {
                if (pms[i].editDistanceTypes[j] == 's')
                {
                    numberOfSub++;
                } else if (pms[i].editDistanceTypes[j] == 'i' || pms[i].editDistanceTypes[j] == 'd')
                {
                    numberOfInDel++;
                }
            }
            pms[i].numberOfSub = numberOfSub;
            pms[i].numberOfInDel = numberOfInDel;
            pms[i].cigar = ldc.alignmentToCIGAR(pms[i].editDistanceTypes);
            pms[i].numberOfMatches = s.countMatches(pms[i].cigar);
            if (pms[i].numberOfMatches > bestAln.numberOfMatches)
            {
                bestAln = pms[i];
            }
            else if (pms[i].numberOfMatches == bestAln.numberOfMatches)
            {
                if (pms[i].numberOfInDel < bestAln.numberOfInDel && pms[i].numberOfSub <= bestAln.numberOfSub + 1) // this is a preference on sub over InDel
                {
                    bestAln = pms[i];
                } else if (pms[i].numberOfInDel == bestAln.numberOfInDel && pms[i].numberOfSub < bestAln.numberOfSub)
                {
                    bestAln = pms[i];
                }
            }
        }
        return bestAln;
    }

    string concatenateNTimes(char ch, int n)
    {
        string result;
        result.reserve(n); // Preallocate memory to avoid reallocation

        for (int i = 0; i < n; ++i)
        {
            result += ch; // Concatenate the character to the string
        }

        return result;
    }

    uint32_t hammingDistanceWithMarkers(LevAlign& hamLa) {
        // Check if the sequences have the same length
        if (hamLa.read.length() != hamLa.query.length()) {
            cerr << "Error: Sequences have different lengths." << endl;
            return 100; 
        }

        // Iterate through the sequences and compare corresponding characters
        for (size_t i = 0; i < hamLa.read.length(); ++i) {
            if (hamLa.read[i] != hamLa.query[i]) {
                hamLa.editDistanceTypes += 's'; // Indicates substitution
                hamLa.editDistance++;
            } else {
                hamLa.editDistanceTypes += '-'; // Indicates exact match
            }
        }
 
        return hamLa.editDistance;
    }

    LevAlign extendSeed(string query, string read, vector<kmerT> regionKmers, uint32_t allowedEditDistance, uint32_t contigSize, char region)
    {
        vector<kmerT> readKmers = extractReadKmers(read);
        LevAlign maxLa;
        LevDistanceCalculator ldc;
        vector<pair<unsigned int, unsigned int>> commonPositions = findCommonKmers(regionKmers, readKmers, region, contigSize);
        // cout << "number of common positions: " << commonPositions.size() << endl;
        if (!commonPositions.empty())
        {
            map<uint32_t, uint32_t> readMemStartEnd, queryMemStartEnd;
            uint32_t skippedCommonPositions = 0;
            for (const auto &pos : commonPositions)
            {
                // cout << "Sequence 1 position: " << pos.first << ", Sequence 2 position: " << pos.second << endl;
                // if (regionKmers[pos.first] != readKmers[pos.second])
                // {
                //     cerr << "wrong k-mer compare occured!!!" << endl;
                //     continue;
                // }
                uint32_t readMemStartInd = 0, queryMemStartInd = 0, readMemEndInd = 0, queryMemEndInd = 0;
                // cout << "----MEM extension-----" << endl;
                string memStr = extendMEM(pos.second, pos.first, read, query, readMemStartInd, queryMemStartInd, readMemEndInd, queryMemEndInd, contigSize);
                //check for skipping rechecking the same MEM
                auto it_r = readMemStartEnd.find(readMemStartInd);
                auto it_q = queryMemStartEnd.find(queryMemStartInd);
                if (it_r != readMemStartEnd.end() && it_q != queryMemStartEnd.end())
                {
                    if (it_r->second == readMemEndInd && it_q->second == queryMemEndInd)
                    {
                        skippedCommonPositions++;
                        continue;
                    } 

                } else
                {
                    readMemStartEnd[readMemStartInd] = readMemEndInd;
                    queryMemStartEnd[queryMemStartInd] = queryMemEndInd;
                }
                
                // cout << "query seed start: " << pos.first << ", read seed start: " << pos.second << endl;
                // cout << "query MEM start: " << readMemStartInd << ", end: " << readMemEndInd << endl;
                // cout << "query MEM start: " << queryMemStartInd << ", end: " << queryMemEndInd << endl;
                // cout << "MEM size: " << memStr.size() << ", MEM = " << memStr << endl;
                LevAlign memLa;
                memLa.partialMatchSize = memStr.size();
                memLa.editDistanceTypes = concatenateNTimes('-', memStr.size());
                memLa.alignedQuery = memStr;
                memLa.alignedRead = memStr;
                memLa.read = read;
                memLa.query = query;
                memLa.cigar = ldc.alignmentToCIGAR(memLa.editDistanceTypes);
                memLa.readRegionStartPos = readMemStartInd;
                memLa.readRegionEndPos = readMemEndInd;
                memLa.queryRegionStartPos = queryMemStartInd;
                memLa.queryRegionEndPos = queryMemEndInd;
                // cout << "MEM edit types: " << memLa.editDistanceTypes << endl;
                if (memStr.size() == contigSize)
                {
                    return memLa;
                }
                /*	Finding edits in prefix MEM region*/
                // cout << "=====preMEM extension=====" << endl;
                LevAlign preLa = extendPreMem(readMemStartInd, queryMemStartInd, read, query, allowedEditDistance);
                preLa.partialMatchSize += memStr.size();
                preLa.alignedQuery = preLa.alignedQuery + memStr;
                preLa.alignedRead = preLa.alignedRead + memStr;
                preLa.editDistanceTypes = preLa.editDistanceTypes + memLa.editDistanceTypes;
                // cout << "preLa edit: " << preLa.editDistanceTypes << endl;
                /*	Finding edits in sudfix MEM region*/
                // cout << "=====sufMEM extension=====" << endl;
                LevAlign sufLa = extendSufMem(readMemEndInd, queryMemEndInd, read, query, allowedEditDistance, contigSize);
                sufLa.partialMatchSize += memStr.size();
                sufLa.alignedQuery = memStr + sufLa.alignedQuery;
                sufLa.alignedRead = memStr + sufLa.alignedRead;
                sufLa.editDistanceTypes = memLa.editDistanceTypes + sufLa.editDistanceTypes;
                // cout << "sufLa edit: " << sufLa.editDistanceTypes << endl;
                /*	Finding edits in middle MEM region*/
                if ((readMemStartInd == 0 || queryMemStartInd == 0) && (readMemEndInd == (contigSize - 1) || queryMemEndInd == (contigSize - 1)))
                {
                    // no prefix and suffix
                    // only MEM
                    //TODO: the hamming distance can be inplemented here too
                    vector<LevAlign> pms;
                    pms.push_back(memLa);
                    pms.push_back(maxLa);
                    maxLa = comparePartialMatchRes(pms);
                    // cout << "<<<<only MEM>>>>" << endl;
                }
                else if ((readMemStartInd > 0 && queryMemStartInd > 0) > 0 && (readMemEndInd == (contigSize - 1) || queryMemEndInd == (contigSize - 1)))
                {
                    // just prefix
                    //check the hamming distance
                    LevAlign hamLa;
                    preLa.readRegionEndPos = memLa.readRegionEndPos;
                    preLa.queryRegionEndPos = memLa.queryRegionEndPos;
                    hamLa.read = read.substr(preLa.readRegionStartPos, (preLa.readRegionEndPos - preLa.readRegionStartPos + 1));
                    hamLa.query = query.substr(preLa.queryRegionStartPos, (preLa.queryRegionEndPos - preLa.queryRegionStartPos + 1));
                    uint32_t ed = hammingDistanceWithMarkers(hamLa);
                    vector<LevAlign> pms;
                    pms.push_back(preLa);
                    if(ed <= allowedEditDistance)
                        pms.push_back(hamLa);
                    pms.push_back(maxLa);
                    maxLa = comparePartialMatchRes(pms);
                    // cout << "<<<<only prefix>>>>" << endl;
                }
                else if ((readMemStartInd == 0 || queryMemStartInd == 0) && (readMemEndInd < (contigSize - 1) && queryMemEndInd < (contigSize - 1)))
                {
                    // just suffix
                    //check the hamming distance
                    LevAlign hamLa;
                    sufLa.readRegionStartPos = memLa.readRegionStartPos;
                    sufLa.queryRegionStartPos = memLa.queryRegionStartPos;
                    hamLa.read = read.substr(sufLa.readRegionStartPos, (sufLa.readRegionEndPos - sufLa.readRegionStartPos + 1));
                    hamLa.query = query.substr(sufLa.queryRegionStartPos, (sufLa.queryRegionEndPos - sufLa.queryRegionStartPos + 1));
                    uint32_t ed = hammingDistanceWithMarkers(hamLa);
                    vector<LevAlign> pms;
                    pms.push_back(sufLa);
                    if(ed <= allowedEditDistance)
                        pms.push_back(hamLa);
                    pms.push_back(maxLa);
                    maxLa = comparePartialMatchRes(pms);
                    // cout << "<<<<only suffix>>>>" << endl;
                }
                else if ((readMemStartInd > 0 && queryMemStartInd > 0) && (readMemEndInd < (contigSize - 1) && queryMemEndInd < (contigSize - 1)))
                {
                    // both prefix and suffix present
                    // cout << "----middle calculation-----" << endl;
                    int preMemRegionSizeDiff = ((readMemStartInd > queryMemStartInd) ? (readMemStartInd - queryMemStartInd) : (queryMemStartInd - readMemStartInd));
                    // int memRegionStart = ((readMemStartInd>queryMemStartInd)?readMemStartInd:queryMemStartInd) - 1;
                    LevAlign midLa = findMaximalPartialMatch(allowedEditDistance, (contigSize - preMemRegionSizeDiff), preLa, sufLa, memLa);\
                    // cout << "<<<<only midle>>>>" << endl;
                    //check the hamming distance
                    LevAlign hamLa;
                    if(readMemStartInd > queryMemStartInd)
                    {
                        hamLa.read = read.substr(preMemRegionSizeDiff, (contigSize - preMemRegionSizeDiff));
                        hamLa.query = query.substr(0, (contigSize - preMemRegionSizeDiff));
                    } else
                    {
                        hamLa.read = read.substr(0, (contigSize - preMemRegionSizeDiff));
                        hamLa.query = query.substr(preMemRegionSizeDiff, (contigSize - preMemRegionSizeDiff));
                    }
                    uint32_t ed = hammingDistanceWithMarkers(hamLa);
                    vector<LevAlign> pms;
                    pms.push_back(midLa);
                    if(ed <= allowedEditDistance)
                        pms.push_back(hamLa);  
                    pms.push_back(maxLa);
                    maxLa = comparePartialMatchRes(pms);
                }
            }
            // cout << "skippedCommonPositions : " << skippedCommonPositions << endl;
        }
        else
        {
            // cout << "No common k-mers found." << endl;
        }
        maxLa.read = read;
        maxLa.query = query;
        maxLa.cigar = ldc.alignmentToCIGAR(maxLa.editDistanceTypes);
        return maxLa;
    }
    
    map<contigIndT, string> readContigsFromMap(tsl::robin_map<uint32_t, string>& reads,  set<contigIndT>& contigIdSet)
    {
        contigIndT setElement;
        unsigned int setInd = 0;
        map<contigIndT, string> readContigs;
        while (true)
        {        
            if (setInd >= 0 && setInd < contigIdSet.size())
            {
                auto it = contigIdSet.begin();
                advance(it, setInd);
                setElement = *it;
                setInd++;
                readContigs[setElement] = reads[setElement];
            }
            else
            {
                break;
            }
        }
        return readContigs;
    }

    bool checkAlingmentCriteria(LevAlign l, uint32_t maxAllowedEdit)
    {
        return ((l.editDistance <= maxAllowedEdit) && (((uint32_t)l.partialMatchSize >= (uint32_t)regionSize) || ((uint32_t)(l.numberOfMatches + l.numberOfInDel + l.numberOfSub) >= (uint32_t)regionSize)));
    }

    void findPartiaMatches(tsl::robin_map <uint32_t, string>& reads, map<uint32_t,string>& queries, IndexContainer<contigIndT, contigIndT>& frontMinThCheapSeedReads, IndexContainer<contigIndT, contigIndT>& backMinThCheapSeedReads, contigIndT queryCount, uint32_t allowedEditDistance, uint32_t contigSize, map<contigIndT, LevAlign> &pmres, IndexContainer<contigIndT, contigIndT>& alignments)
    {      
        Utilities<double> utildouble;   
        set<double> seedAndExtend_times;
        uint32_t pmNotFoundMatchQueries = 0;
        for (size_t i = 0; i < queryCount; i++)
        {
            auto itq = queries.find(i);
            if (itq == queries.end())
                continue;
            string query = queries[i];
            if (query.find('n') != string::npos || query.find('N') != string::npos)
                continue;
            // long long readFromFile_ExeTime = 0, seedExtension_ExeTime = 0;
            auto se_start = chrono::high_resolution_clock::now();
            LevAlign bestAlignmentForQuery;
            //check the pmres list (list of best alignments) for the reverse strand to extract the best alignment of the forward strand
            auto it = pmres.find(i);
            if (it != pmres.end())
            {
                bestAlignmentForQuery = it->second;
            }
            
            // check front region of query
            // read the candidate reads of cheap k-mer filter of front
            auto frontReadSet = frontMinThCheapSeedReads.get(i);
            // auto start = chrono::high_resolution_clock::now();
            // map<contigIndT, string> frontCandidateReads = readContigs(readAddress, frontReadSet);
            map<contigIndT, string> frontCandidateReads = readContigsFromMap(reads, frontReadSet);
            // auto end = chrono::high_resolution_clock::now();
            // auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            // readFromFile_ExeTime += duration.count();
            vector<kmerT> frontKmers;
            extractKmersFromRegions(query, frontKmers, frontRegion);
            LevAlign frontLa;
            // start = chrono::high_resolution_clock::now();
            for (auto it = frontCandidateReads.begin(); it != frontCandidateReads.end(); it++)
            {
                frontLa = extendSeed(query, it->second, frontKmers, allowedEditDistance, contigSize, frontRegion);
                frontLa.readID = it->first;
                vector<LevAlign> pms;
                pms.push_back(frontLa);
                pms.push_back(bestAlignmentForQuery);
                bestAlignmentForQuery = comparePartialMatchRes(pms);
                if(checkAlingmentCriteria(frontLa, allowedEditDistance))
                    alignments.put(i, it->first);
            }
            // end = chrono::high_resolution_clock::now();
            // duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            // seedExtension_ExeTime += duration.count();
            // check back region of query
            // read the candidate reads of cheap k-mer filter
            auto backReadSet = backMinThCheapSeedReads.get(i);
            // start = chrono::high_resolution_clock::now();
            // map<contigIndT, string> backCandidateReads = readContigs(readAddress, backReadSet);
            map<contigIndT, string> backCandidateReads = readContigsFromMap(reads, backReadSet);
            // end = chrono::high_resolution_clock::now();
            // duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            // readFromFile_ExeTime += duration.count();
            vector<kmerT> backKmers;
            extractKmersFromRegions(query, backKmers, backRegion);
            LevAlign backLa;
            // start = chrono::high_resolution_clock::now();
            for (auto it = backCandidateReads.begin(); it != backCandidateReads.end(); it++)
            {
                backLa = extendSeed(query, it->second, backKmers, allowedEditDistance, contigSize, backRegion);
                backLa.readID = it->first;
                vector<LevAlign> pms;
                pms.push_back(backLa);
                pms.push_back(bestAlignmentForQuery);
                bestAlignmentForQuery = comparePartialMatchRes(pms);
                if(checkAlingmentCriteria(backLa, allowedEditDistance))
                    alignments.put(i, it->first);
            }
            // end = chrono::high_resolution_clock::now();
            // duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            // seedExtension_ExeTime += duration.count();
            auto se_end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::nanoseconds>(se_end - se_start);
            long long seedAndExtendForQuery_microsecond = static_cast<double>(duration.count())  / 1'000.0;
            seedAndExtend_times.insert((double)(seedAndExtendForQuery_microsecond));
            // long long seedExtensionForQuery_second = static_cast<double>(seedExtension_ExeTime)  / 1'000'000'000.0;
            // long long readFromFileForQuery_second = static_cast<double>(readFromFile_ExeTime)  / 1'000'000'000.0;
            //report best alignment
            if (isVerboseLog_) cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            if (isVerboseLog_) cout << "Q : " << bestAlignmentForQuery.query << endl;
            if (isVerboseLog_) cout << "R : " << bestAlignmentForQuery.read << endl;
            if (isVerboseLog_) cout << "q : " << bestAlignmentForQuery.alignedQuery << endl;
            if (isVerboseLog_) cout << "r : " << bestAlignmentForQuery.alignedRead << endl;
            if (isVerboseLog_) cout << "E : " << bestAlignmentForQuery.editDistanceTypes << endl;
            if (isVerboseLog_) cout << "partial match size : " << bestAlignmentForQuery.partialMatchSize << ", ed : " << bestAlignmentForQuery.editDistance << ", CIGAR : " << bestAlignmentForQuery.cigar << endl;
            if (isVerboseLog_) cout << "Q: " << i << ", R: " << bestAlignmentForQuery.readID << endl;
            if (isVerboseLog_) cout << "seed find and extension for this q took : " << seedAndExtendForQuery_microsecond << " microseconds" << endl;
            // cout << "seed extension for this q took : " << seedExtensionForQuery_second << " seconds" << endl;
            // cout << "read from file for this q took : " << readFromFileForQuery_second << " seconds" << endl;
            if(bestAlignmentForQuery.readID >= 0 && checkAlingmentCriteria(bestAlignmentForQuery, allowedEditDistance))
                pmres[i] = bestAlignmentForQuery;
            else
            {
                pmNotFoundMatchQueries++;
                // cout << "not found match for this query!" << endl;
            }
        }
        tuple<double, double, double> seedAndExtend_timesTuple = utildouble.calculateStatistics(seedAndExtend_times);
        printf("seed find and extension time (micro second) for queries up to q(%d) => [average: %f, median: %f, mean: %f]\n", queryCount, get<0>(seedAndExtend_timesTuple), get<1>(seedAndExtend_timesTuple), get<2>(seedAndExtend_timesTuple));
        cout << "# of queries that pm did not find match : " << pmNotFoundMatchQueries << endl;
    }
};

#endif