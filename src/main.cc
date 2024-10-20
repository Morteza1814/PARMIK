#include "Includes/IndexContainer.h"
#include "Includes/InvertedIndexBuilder.h"
#include "Includes/KmersFrequencyCounter.h"
#include "Includes/BaseLinePartialMatcher.h"
#include "Includes/CheapKmerPartialMatcher.h"
#include "Includes/args.hxx"
#include "Includes/Configs.h"
#include <iostream>
#include "Includes/Utils.h"
#include "Includes/IndexContainer.h"
#include "Includes/Container.h"
#include "Includes/SeedMatchExtender.h"
#include "Includes/SamReader.h"
#include "Includes/CompareWithBWA.h"
#include "Includes/CompareWithBlast.h"
#include "Includes/CompareWithGT.h"
#include "Includes/CheckKmersFrequency.h"
#include "Includes/Aligner.h"
#include "Includes/Alignment.h"
#include "Includes/SSW_BaseLine.h"
#include "Includes/CompareWithBaseLine.h"
#include "Includes/ExpensiveKmersFNEvaluator.h"
#include "Includes/EvaluateSecondChance.h"
#include <filesystem>
#include <cmath>

#define PARMIK_MODE_INDEX           0
#define PARMIK_MODE_ALIGN           1
#define PARMIK_MODE_COMPARE         2
#define PARMIK_MODE_BASELINE        3
#define PARMIK_MODE_CMP_BASELINE    4

namespace fs = std::filesystem;

int argParse(int argc, char** argv, Config &cfg){
	args::ArgumentParser parser("=========================Arguments===========================", "======================================================");
    args::ValueFlag<int> parmikModeArg(parser, "", "PARMIK mode",                                   {'a', "mode"});
    args::ValueFlag<string> otherToolAddressArg(parser, "", "Other Tool Alignment File Address",    {'b', "toolFileAddress"});
	args::ValueFlag<int> contigSizeArg(parser, "", "Contig Size",                                   {'c', "contigSize"});
    args::ValueFlag<int> percentageIdentityArg(parser, "", "Percentage Identity",                   {'d', "percentageIdentity"});
	args::ValueFlag<int> editDistanceArg(parser, "", "Max Edit Distance (i/d/s)",                   {'e', "editDistance"});
	args::ValueFlag<string> offlineIndexAddressArg(parser, "", "Inexpensive K-mer Index Address",   {'f', "ikiAddress"});
    args::HelpFlag help(parser, "help", "Help",                                                     {'h', "help"});
	args::ValueFlag<int> readsCountArg(parser, "", "Number of Metageomic Reads",                    {'i', "readCount"});
	args::ValueFlag<int> queryCountArg(parser, "", "Number of Queries",                             {'j', "queryCount"});
	args::ValueFlag<int> kmerLengthArg(parser, "", "K-mer Length",                                  {'k', "kmerLen"});
    args::ValueFlag<string> otherToolArg(parser, "", "The Other Tool  (bwa, blast, etc)",           {'l', "otherTool"});
    args::ValueFlag<int> minExactMatchLenArg(parser, "", "Minimum Exact Match Length",              {'m', "minExactMatchLen"});
    args::ValueFlag<string> kmerRangesFileAddressArg(parser, "", "K-mer Ranges File Address",       {'n', "kmerRangesFileAddress"});
	args::ValueFlag<string> outputDirArg(parser, "", "OutputDir",                                   {'o', "outputDir"});
    args::ValueFlag<string> penaltyFileAddressArg(parser, "", "Penalty File Address",               {'p', "penaltyFileAddress"});
	args::ValueFlag<string> queryFileAddressArg(parser, "", "Query File Address",                   {'q', "query"});
	args::ValueFlag<string> readDatabaseAddressArg(parser, "", "Metageomic Read Data Base Address", {'r', "read"});
	args::ValueFlag<int> regionSizeArg(parser, "", "Region Size",                                   {'s', "regionSize"});
    args::ValueFlag<int> cheapKmerThresholdArg(parser, "", "Cheap Kmer Threshold",                  {'t', "cheapKmerThreshold"});
    args::Flag isSecondChanceOff(parser, "", "Turn Second Chance Off",                              {'u', "isSecondChanceOff"});
	args::Flag isVerboseLogArg(parser, "", "Verbose Logging",                                       {'v', "verboseLog"});
    args::ValueFlag<int> numThreadsArg(parser, "", "Number of Threads",                             {'w', "numThreads"});
	args::Flag isIndexOfflineArg(parser, "", "Is the read index offline",                           {'x', "isIndexOffline"});
    args::ValueFlag<string> baselineBaseAddressArg(parser, "", "BaseLine file base address",        {'z', "baselineBaseAddress"});
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e)
    {
        cout << e.what();
        return 0;
    }
    catch (const args::Help&)
    {
        cout << parser;
        return 0;
    }
    catch (const args::ParseError& e)
    {
        cerr << e.what() << endl;
        cerr << parser;
        return 0;
    }
    if (parmikModeArg) {cfg.parmikMode = args::get(parmikModeArg); } else {cout << "no parmik mode is determined!"<< endl; return 0;}
	if (readDatabaseAddressArg) {cfg.readDatabaseAddress = args::get(readDatabaseAddressArg); } else {cout << "no readDatabaseAddress!"<< endl; return 0;}
	if (queryFileAddressArg) {cfg.queryFileAddress = args::get(queryFileAddressArg);} else {cout << "no queryFileAddress!"<< endl; return 0;}
    if (penaltyFileAddressArg) {cfg.penaltyFileAddress = args::get(penaltyFileAddressArg);} else {cfg.penaltyFileAddress = ""; cout << "no penaltyFileAddress!"<< endl;}
	if (outputDirArg) {cfg.outputDir = args::get(outputDirArg);} else {cout << "no outputDirArg!"<< endl; if(cfg.parmikMode != PARMIK_MODE_INDEX) return 0;}
    if (kmerRangesFileAddressArg) {cfg.kmerRangesFileAddress = args::get(kmerRangesFileAddressArg);} else {cout << "no kmerRangesFileAddressArg!"<< endl; cfg.kmerRangesFileAddress = "";}
    if (percentageIdentityArg) {cfg.percentageIdentity = args::get(percentageIdentityArg);cfg.percentageIdentity = (double)(cfg.percentageIdentity/100);} else {cout << "no percentageIdentity (default = 0.9)!"<< endl;}
	if (readsCountArg) {cfg.readsCount = args::get(readsCountArg); } else {cfg.readsCount = NUMBER_OF_READS;}
	if (queryCountArg) {cfg.queryCount = args::get(queryCountArg); } else {cfg.queryCount = NUMBER_OF_QUERIES;}
	if (kmerLengthArg) {cfg.kmerLength = args::get(kmerLengthArg); } else {cfg.kmerLength = KMER_SZ;}
	if (regionSizeArg) {cfg.regionSize = args::get(regionSizeArg); } else {cfg.regionSize = REGION_SZ;}
	if (contigSizeArg) {cfg.contigSize = args::get(contigSizeArg); } else {cfg.contigSize = CONTIG_SZ;}
    if (cheapKmerThresholdArg) {cfg.cheapKmerThreshold = args::get(cheapKmerThresholdArg); } else {cout << "no cheapKmerThreshold!"<< endl; return 0;}
    if (numThreadsArg) {cfg.numThreads = args::get(numThreadsArg); } else {cfg.numThreads = 1;};
    if (isIndexOfflineArg) {cfg.isIndexOffline = true; } else {cfg.isIndexOffline = false;}
    if (isVerboseLogArg) {cfg.isVerboseLog = true; } else {cfg.isVerboseLog = false;}
    if (isSecondChanceOff) {cfg.isSecondChanceOff = true; } else {cfg.isSecondChanceOff = false;}
    if (offlineIndexAddressArg) {cfg.offlineIndexAddress = args::get(offlineIndexAddressArg); } else {cout << "no offlineIndexAddress!"<< endl; return 0;}
    if (otherToolAddressArg) {cfg.otherToolOutputFileAddress = args::get(otherToolAddressArg); } else {cout << "no otherToolOutputFileAddress!"<< endl; if(cfg.parmikMode == PARMIK_MODE_CMP_BASELINE || cfg.parmikMode == PARMIK_MODE_COMPARE) return 0;}
    if (baselineBaseAddressArg) {cfg.baselineBaseAddress = args::get(baselineBaseAddressArg); } else {cout << "no baselineBaseAddress!"<< endl; if(cfg.parmikMode == PARMIK_MODE_CMP_BASELINE){return 0;}}
    if (otherToolArg) {cfg.otherTool = args::get(otherToolArg);} else {cout << "no otherToolArg!"<< endl;if(cfg.parmikMode == PARMIK_MODE_COMPARE) return 0;}
	if (editDistanceArg) {cfg.editDistance = args::get(editDistanceArg); } else {cfg.editDistance = NUMBER_OF_ALLOWED_EDIT_DISTANCES;}
    if (minExactMatchLenArg) {cfg.minExactMatchLen = args::get(minExactMatchLenArg); } else {cfg.minExactMatchLen = 0;}
	cfg.readFileName = cfg.readDatabaseAddress.substr(cfg.readDatabaseAddress.find_last_of("/\\") + 1);
	cfg.queryFileName = cfg.queryFileAddress.substr(cfg.queryFileAddress.find_last_of("/\\") + 1);
    assert(cfg.regionSize > cfg.kmerLength && "region size should be larger than k-mer length");
	return 1;
}

set<uint32_t> convertMapToSet(const map<uint32_t, uint32_t>& inputMap) {
    set<uint32_t> resultSet;

    for (const auto& pair : inputMap) {
        resultSet.insert(pair.second);
    }

    return resultSet;
}

map<uint32_t, uint32_t> checkRates(IndexContainer<uint32_t, uint32_t> baselineSeeds, IndexContainer<uint32_t, uint32_t> experimentSeeds, map<uint32_t, uint32_t> &seedsFalseRate)
{
    map<uint32_t, uint32_t> seedsTrueRate;
    for(auto it = baselineSeeds.begin(); it != baselineSeeds.end(); it++)
    {
        bool found = false;
        if (seedsFalseRate.find(it->first) == seedsFalseRate.end())
        {
            seedsFalseRate[it->first] = 0;
        }
        if (seedsTrueRate.find(it->first) == seedsTrueRate.end())
        {
            seedsTrueRate[it->first] = 0;
        }
        auto range = experimentSeeds.getRange(it->first);
        for (auto itt = range.first; itt != range.second; itt++) 
        {
            if (it->second == itt->second)
            {
                found = true;
                break;
            } 
        }
        if (!found)
        {
            // cout << "false item : first = " << it->first << ", second = " << it->second << endl;
            seedsFalseRate[it->first]++;
        } else 
        {
            seedsTrueRate[it->first]++;
        }
    }
    return seedsTrueRate;
}

uint32_t sumMapValues(const map<uint32_t, uint32_t>& inputMap) {
    uint32_t sum = 0; // Initialize sum to the default value of Value

    for (const auto& pair : inputMap) {
        sum += pair.second;
    }

    return sum;
}

void convertSamToLev(SamReader::Sam parmikSamAlignment, LevAlign& parmikAlignment)
{
    SamReader sam("");
    parmikAlignment.readID = parmikSamAlignment.readId;
    parmikAlignment.queryID = parmikSamAlignment.queryId;
    parmikAlignment.editDistance = parmikSamAlignment.editDistance;
    parmikAlignment.cigar = parmikSamAlignment.cigar;
    parmikAlignment.numberOfSub = parmikSamAlignment.editDistance;
    parmikAlignment.numberOfMatches = sam.countMatches(parmikSamAlignment.cigar);
    parmikAlignment.numberOfInDel = sam.countInsertions(parmikSamAlignment.cigar) + sam.countDeletions(parmikSamAlignment.cigar);
    parmikAlignment.flag = parmikSamAlignment.flag;                
}

vector<Penalty> readPenalties(const string& readPenaltiesFileAddress)
{
    vector<Penalty> penalties;
    ifstream readPenaltiesFile(readPenaltiesFileAddress);
    if (!readPenaltiesFile.is_open())
    {
        cout << "Error opening file " << readPenaltiesFileAddress << endl;
        Penalty p;
        penalties.push_back(p);
        cout << "--------------Penalties--------------" << endl;
        cout << "m: " << p.matchPenalty << ", s: " << p.mismatchPenalty << ", go: " << p.gapOpenPenalty << ", ge: " << p.gapExtendPenalty << endl;
        cout << "------------------------------------" << endl;
        return penalties;
    }
    string line;
    cout << "--------------Penalties--------------" << endl;
    while (getline(readPenaltiesFile, line))
    {
        Penalty p;
        stringstream ss(line);
        ss >> p.matchPenalty;
        ss >> p.mismatchPenalty;
        ss >> p.gapOpenPenalty;
        ss >> p.gapExtendPenalty;
        cout << "m: " << p.matchPenalty << ", s: " << p.mismatchPenalty << ", go: " << p.gapOpenPenalty << ", ge: " << p.gapExtendPenalty << endl;
        penalties.push_back(p);
    }
    readPenaltiesFile.close();
    if (penalties.size() == 0)
    {
        Penalty p;
        cout << "m: " << p.matchPenalty << ", s: " << p.mismatchPenalty << ", go: " << p.gapOpenPenalty << ", ge: " << p.gapExtendPenalty << endl;
        penalties.push_back(p);
    }
    cout << "------------------------------------" << endl;
    return penalties;
}

string getPenaltiesSubstr(vector<Penalty> penalties)
{
    string substr = "";
    for (auto p : penalties)
    {
        substr += "_" + to_string(p.matchPenalty) + "" + to_string(p.mismatchPenalty) + "" + to_string(p.gapOpenPenalty) + "" + to_string(p.gapExtendPenalty) + "";
    }
    return substr;
}

string getTextAfterLastSlash(const string& filePath) {
    size_t pos = filePath.find_last_of('/');
    if (pos != string::npos) {
        return filePath.substr(pos + 1);
    }
    return filePath; // If no '/' found, return the original string
}

void createDir(string path){
    try {
        if (fs::create_directory(path)) {
            std::cout << "Directory created successfully: " << path << std::endl;
        } else {
            std::cout << "Directory already exists or could not be created: " << path << std::endl;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "General error: " << e.what() << std::endl;
    }
}

int run(int argc, char *argv[]) {
    Config cfg;
	ofstream out;
    Utilities<uint32_t> util;
	if (!argParse(argc, argv, cfg))
    	return 0;
	cout << "**********************CONFIGURATIONS*****************************" << endl;
    cout << left << setw(30) << "PARMIK mode: " << cfg.parmikMode << endl;
	cout << left << setw(30) << "readDatabaseAddress: " << cfg.readDatabaseAddress << endl;
	cout << left << setw(30) << "readFileName: " << cfg.readFileName << endl;
	cout << left << setw(30) << "queryFileAddress: " << cfg.queryFileAddress << endl;
	cout << left << setw(30) << "queryFileName: " << cfg.queryFileName << endl;
    cout << left << setw(30) << "otherToolOutputFileAddress: " << cfg.otherToolOutputFileAddress << endl;
	cout << left << setw(30) << "outputDir: " << cfg.outputDir << endl;
    cout << left << setw(30) << "penaltyFileAddress: " << cfg.penaltyFileAddress << endl;
    cout << left << setw(30) << "kmerRangesFileAddress: " << cfg.kmerRangesFileAddress << endl;
	cout << left << setw(30) << "readsCount: " << cfg.readsCount << endl; 
	cout << left << setw(30) << "queryCount: " << cfg.queryCount << endl;
	cout << left << setw(30) << "kmerLength: " << cfg.kmerLength << endl;
	cout << left << setw(30) << "regionSize: " << cfg.regionSize << endl;
    cout << left << setw(30) << "cheapKmerThreshold: " << cfg.cheapKmerThreshold << endl;
    cout << left << setw(30) << "percentageIdentity: " << cfg.percentageIdentity << endl;
    cout << left << setw(30) << "minExactMatchLen: " << cfg.minExactMatchLen << endl;
    cout << left << setw(30) << "numThreads: " << cfg.numThreads << endl;
    cout << left << setw(30) << "isSecondChanceOff: " << ((cfg.isSecondChanceOff == 1) ? "T":"F") << endl;
    uint32_t minNumExactMatchKmer = 0;
    if (cfg.parmikMode != PARMIK_MODE_INDEX && cfg.parmikMode != PARMIK_MODE_BASELINE) {
        if(cfg.minExactMatchLen > 0) {
            minNumExactMatchKmer = cfg.minExactMatchLen - (cfg.kmerLength - 1);
        } else {
            cfg.editDistance = (uint32_t)(round(cfg.regionSize * (1-cfg.percentageIdentity)));
            uint16_t exactMatchSize = cfg.regionSize - cfg.editDistance;
            uint32_t regionKmers = (uint32_t)(floor(exactMatchSize/cfg.kmerLength));
            cout << "regionKmers: " << regionKmers << ", exactMatchSize: " << exactMatchSize << ", editDistance: " << cfg.editDistance << endl;
            assert(regionKmers > cfg.editDistance && "k-mers in a region must be greater than E for cheap k-mer matching");
            minNumExactMatchKmer = (uint32_t)(regionKmers + floor(exactMatchSize%cfg.kmerLength));
        }
        assert(minNumExactMatchKmer > 0 && "minNumExactMatchKmer must be greater than 0 for cheap k-mer matching");
    }
    cout << left << setw(30) << "minNumExactMatchKmer: " << minNumExactMatchKmer << endl;
    cout << left << setw(30) << "isIndexOffline: " << cfg.isIndexOffline << endl;
    cout << left << setw(30) << "offlineIndexAddress: " << cfg.offlineIndexAddress << endl;
    cout << left << setw(30) << "baselineBaseAddress: " << cfg.baselineBaseAddress << endl;
    cout << left << setw(30) << "max editDistance in min region: " << cfg.editDistance << endl;
    cout << left << setw(30) << "otherTool: " << cfg.otherTool << endl;
	cout << "**********************CONFIGURATIONS*****************************" << endl;

    try { 
        //construct the baseline
        // IndexContainer<uint32_t, uint32_t> baselineQueriesFrontSeeds; 
        // IndexContainer<uint32_t, uint32_t> baselineQueriesBackSeeds;
        // BaseLinePartialMatcher<uint32_t, uint64_t, uint32_t> blpm(cfg.minExactMatchLen, cfg.regionSize);
        // blpm.constructBaseline(cfg.readDatabaseAddress, cfg.readsCount, cfg.queryFileAddress, cfg.queryCount,cfg.minExactMatchLen, cfg.cheapKmerThreshold, baselineQueriesFrontSeeds, baselineQueriesBackSeeds, cfg.isIndexOffline, cfg.offlineIndexAddress);
        // run the experiment
        Container<uint32_t, uint32_t> minThCheapSeedReads, revMinThCheapSeedReads;
        // map<uint32_t, LevAlign> pmr;
        tsl::robin_map <uint32_t, string> reads, queries;
        uint32_t queryCount = 0;
        vector<Penalty> penalties = readPenalties(cfg.penaltyFileAddress);
        // experiment dir name
        string expDirName = "IKT" + to_string(cfg.cheapKmerThreshold) + "_K" + to_string(cfg.kmerLength) + "_PI" + to_string((uint32_t)floor(cfg.percentageIdentity*100)) + "_M" + to_string(minNumExactMatchKmer) + "_T" + to_string(cfg.numThreads) + "_SC" + ((cfg.isSecondChanceOff == 1) ? "0":"1") + "_P" + getPenaltiesSubstr(penalties);
        string expDir = cfg.outputDir + "/" + expDirName;
        cout << "experiment dir: " << expDir << endl;
        uint32_t queryBaseIndex = 0;
        if (cfg.parmikMode != PARMIK_MODE_INDEX)
        {
            if (cfg.parmikMode != PARMIK_MODE_COMPARE) createDir(expDir);
            queryCount = util.readContigsFromFile(cfg.queryFileAddress, cfg.queryCount, queries, queryBaseIndex);
            cout << "queryCount : " << queryCount << endl;
            // read the penalties
        }
        uint32_t readBaseIndex = 0;
        uint32_t readCount = util.readContigsFromFile(cfg.readDatabaseAddress, cfg.readsCount, reads, readBaseIndex);
        cout << "readCount : " << readCount << endl;
        // unordered_map<uint32_t, unordered_set<LevAlign>> alignments;
        string offlineCheapIndexAddress = cfg.offlineIndexAddress + "ck_T" + to_string(cfg.cheapKmerThreshold) + "_K"+ to_string(cfg.kmerLength) + "_r" + to_string(cfg.readsCount);
        string offlineExpensiveIndexAddress = cfg.offlineIndexAddress + "ek_T" + to_string(cfg.cheapKmerThreshold) + "_K"+ to_string(cfg.kmerLength) + "_r" + to_string(cfg.readsCount);
        cout << "offlineCheapIndexAddress: " << offlineCheapIndexAddress << endl;
        cout << "offlineExpensiveIndexAddress: " << offlineExpensiveIndexAddress << endl;
        // string parmikAlignmentsAddress = cfg.outputDir + "/aln/pmAln_" + "R" + to_string(cfg.regionSize) + "_PI" + to_string((uint32_t)floor(cfg.percentageIdentity*100)) + "_L" + to_string(cfg.minExactMatchLen) + "_M" + to_string(minNumExactMatchKmer) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
        string parmikAlignmentsAddress = expDir + "/" + "parmikAln.txt";
        string parmikExpensiveKmerFNsAddress = expDir + "/" + "parmikExpensiveKmerFNs.txt";
        cout << "pmrkAlignmentsAddress: " << parmikAlignmentsAddress << endl;
        //baseline
        if (cfg.parmikMode == PARMIK_MODE_BASELINE) {
            string queryFileName = getTextAfterLastSlash(cfg.queryFileAddress);
            if(queryFileName == cfg.queryFileAddress){
                cerr << "Error: query file name is not valid" << endl;
                return 1;
            }
            string baselineAlignmentsAddress = cfg.outputDir + "/BL_Aln" + "_RS" + to_string(cfg.kmerLength) + "_PI" + to_string((uint32_t)floor(cfg.percentageIdentity*100)) + "_P" + getPenaltiesSubstr(penalties) + "_Q" + queryFileName + ".txt";
            SSW_BaseLine aligner(cfg.kmerLength, cfg.percentageIdentity); 
            aligner.findPartiaMatches(reads, queries, queryCount, true, baselineAlignmentsAddress, penalties, queryBaseIndex);
            //do it again for the reverse strand
            tsl::robin_map <uint32_t, string> revQueries = util.reverseComplementMapValues(queries);
            aligner.findPartiaMatches(reads, revQueries, queryCount, false, baselineAlignmentsAddress, penalties, queryBaseIndex);
        } else if (cfg.parmikMode == PARMIK_MODE_INDEX || cfg.parmikMode == PARMIK_MODE_ALIGN) {
            if (cfg.kmerLength <= 16)
            {
                Container<uint32_t, uint32_t> cheapKmers, expensiveKmers;
                if(!cfg.isIndexOffline)
                {
                    //create the partial matching inverted read index
                    InvertedIndexBuilder<uint32_t, uint32_t> builder(cfg.kmerLength, cfg.kmerLength - 1);
                    // Build the inverted index
                    IndexContainer<uint32_t, uint32_t> invertedIndex = builder.build(cfg.readDatabaseAddress, cfg.readsCount, cfg.isIndexOffline);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // invertedIndex.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    //collect cheap kmers
                    KmersFrequencyCounter<uint32_t, uint32_t> kfc(cfg.cheapKmerThreshold);
                    kfc.collectCheapKmers(cheapKmers, expensiveKmers, invertedIndex, cfg.kmerRangesFileAddress);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // cheapKmers.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    invertedIndex.clear();
                    clock_t ser_start_time = clock();
                    cheapKmers.serialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Serialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                    if(cfg.cheapKmerThreshold > 0) expensiveKmers.serialize(offlineExpensiveIndexAddress);
                } else 
                {
                    clock_t ser_start_time = clock();
                    cheapKmers.deserialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Deserialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                    if(cfg.cheapKmerThreshold > 0) expensiveKmers.deserialize(offlineExpensiveIndexAddress);
                }
                cout << "-------End of loading k-mer indices-------" << endl;
                if (cfg.parmikMode == PARMIK_MODE_ALIGN)
                {
                    //do partial matching based on cheap k-mers
                    CheapKmerPartialMatcher<uint32_t, uint32_t, uint32_t> ckpm(cfg.kmerLength, cfg.contigSize, minNumExactMatchKmer, cfg.isVerboseLog);
                    ckpm.cheapSeedFilter(cheapKmers, expensiveKmers, queries, minThCheapSeedReads, parmikExpensiveKmerFNsAddress);
                    //get the reverse complement of queries
                    tsl::robin_map <uint32_t, string> revQueries = util.reverseComplementMapValues(queries);
                    ckpm.cheapSeedFilter(cheapKmers, expensiveKmers, revQueries, revMinThCheapSeedReads, parmikExpensiveKmerFNsAddress);
                    // ckpm.printArrays();
                    // SeedMatchExtender<uint32_t, uint64_t> pm(cfg.minExactMatchLen, cfg.regionSize, cfg.isVerboseLog, cfg.editDistance, cfg.contigSize, cfg.inDelPenalty, cfg.subPenalty);
                    Aligner <uint32_t> aligner(cfg.regionSize, cfg.editDistance, cfg.contigSize, cfg.kmerLength, minNumExactMatchKmer, cfg.percentageIdentity, cfg.isSecondChanceOff, cfg.numThreads);
                    aligner.findPartiaMatches(reads, queries, minThCheapSeedReads, queryCount, true, parmikAlignmentsAddress, penalties);
                    //do it again for the reverse strand
                    aligner.findPartiaMatches(reads, revQueries, revMinThCheapSeedReads, queryCount, false, parmikAlignmentsAddress, penalties);
                }
            } else
            {
                Container<uint64_t, uint32_t> cheapKmers, expensiveKmers;
                if(!cfg.isIndexOffline)
                {
                    //create the partial matching inverted read index
                    InvertedIndexBuilder<uint64_t, uint32_t> builder(cfg.kmerLength, cfg.kmerLength - 1);
                    // Build the inverted index
                    IndexContainer<uint64_t, uint32_t> invertedIndex = builder.build(cfg.readDatabaseAddress, cfg.readsCount, cfg.isIndexOffline);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // invertedIndex.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    //collect cheap kmers
                    KmersFrequencyCounter<uint64_t, uint32_t> kfc(cfg.cheapKmerThreshold);
                    kfc.collectCheapKmers(cheapKmers, expensiveKmers, invertedIndex, cfg.kmerRangesFileAddress);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // cheapKmers.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    invertedIndex.clear();
                    clock_t ser_start_time = clock();
                    cheapKmers.serialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Serialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                    if(cfg.cheapKmerThreshold > 0) expensiveKmers.serialize(offlineExpensiveIndexAddress);
                } else 
                {
                    clock_t ser_start_time = clock();
                    cheapKmers.deserialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Deserialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                    if(cfg.cheapKmerThreshold > 0) expensiveKmers.deserialize(offlineExpensiveIndexAddress);
                }
                cout << "-------End of loading k-mer indices-------" << endl;
                if (cfg.parmikMode == PARMIK_MODE_ALIGN)
                {
                    // do partial matching based on cheap k-mers
                    CheapKmerPartialMatcher<uint32_t, uint64_t, uint32_t> ckpm(cfg.kmerLength, cfg.contigSize, minNumExactMatchKmer, cfg.isVerboseLog);
                    ckpm.cheapSeedFilter(cheapKmers, expensiveKmers, queries, minThCheapSeedReads, parmikExpensiveKmerFNsAddress);
                    //get the reverse complement of queries
                    tsl::robin_map <uint32_t, string> revQueries = util.reverseComplementMapValues(queries);
                    ckpm.cheapSeedFilter(cheapKmers, expensiveKmers, revQueries, revMinThCheapSeedReads, parmikExpensiveKmerFNsAddress);
                    // ckpm.printArrays();
                    // SeedMatchExtender<uint32_t, uint64_t> pm(cfg.minExactMatchLen, cfg.regionSize, cfg.isVerboseLog, cfg.editDistance, cfg.contigSize, cfg.inDelPenalty, cfg.subPenalty);
                    Aligner <uint32_t> aligner(cfg.regionSize, cfg.editDistance, cfg.contigSize, cfg.kmerLength, minNumExactMatchKmer, cfg.percentageIdentity, cfg.isSecondChanceOff, cfg.numThreads);
                    aligner.findPartiaMatches(reads, queries, minThCheapSeedReads, queryCount, true, parmikAlignmentsAddress, penalties);
                    //do it again for the reverse strand
                    aligner.findPartiaMatches(reads, revQueries, revMinThCheapSeedReads, queryCount, false, parmikAlignmentsAddress, penalties);
                }
            }
            if(cfg.parmikMode == PARMIK_MODE_ALIGN){
                //print global variables
                // cout << "----------global variables----------" << endl;
                // cout << "gAlignmentFoundWithNoPolish: " << gAlignmentFoundWithNoPolish << endl;
                // cout << "gAlignmentFoundWithPolish: " << gAlignmentFoundWithPolish << endl;
                // cout << "gAlignmentDumped: " << gAlignmentDumped << endl;
                // cout << "gQueriesHaveAtLeastOneAlignment: " << gQueriesHaveAtLeastOneAlignment << endl;
                // cout << "gAlignmentFoundWithPolishLargerThanBest: " << gAlignmentFoundWithPolishLargerThanBest << endl;
            }
        } else if (cfg.parmikMode == PARMIK_MODE_CMP_BASELINE){
            string queryFileName = getTextAfterLastSlash(cfg.queryFileAddress);
            if(queryFileName == cfg.queryFileAddress){
                cerr << "Error: query file name is not valid" << endl;
                return 1;
            }
            string comparisonResultsFileAddress = expDir + "/" + "cmp_Baseline_" + cfg.otherTool + ".txt";
            string alnReportAddressBase = expDir + "/" + "cmp_Baseline_" + cfg.otherTool + "_";
            string baselineBaseAddress = cfg.baselineBaseAddress + "/BL_Aln" + "_RS" + "11" + "_PI85" + "_P" + getPenaltiesSubstr(penalties) + "_Q" + queryFileName;
            CompareWithBaseLine blCmp(cfg.percentageIdentity);
            if(cfg.otherTool == "parmik" || cfg.otherTool == "PARMIK") {
                cfg.otherToolOutputFileAddress = parmikAlignmentsAddress;
                cout << "updated otherToolOutputFileAddress alignments file: " << parmikAlignmentsAddress << endl;
            }
            blCmp.compareWithBaseLine(cfg, reads, queries, comparisonResultsFileAddress, queryCount, alnReportAddressBase, baselineBaseAddress);
        } else if (cfg.parmikMode == PARMIK_MODE_COMPARE){
            vector<pair<uint32_t, uint32_t>> alnPmVsOtherAlnSizesMap;
            string comparisonResultsFileAddress = expDir + "/" + "cmp_parmik_" + cfg.otherTool + ".txt";
            if(cfg.otherTool == "BWA" || cfg.otherTool == "bwa")
            {
                ComparatorWithBWA cwb(cfg.percentageIdentity);
                cwb.comparePmWithBwa(cfg, reads, queries, queryCount, comparisonResultsFileAddress);
            } else if(cfg.otherTool == "BLAST" || cfg.otherTool == "blast")
            {
                CompareWithBlast cwb(cfg.percentageIdentity);
                cwb.comparePmWithBlast(cfg, reads, queries, queryCount, comparisonResultsFileAddress);
            }
        }

    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << endl;
        return 1;
    }

    return 0;
}

void checkParmikFNalignments(int argc, char *argv[]){
    // checkKmerFreq(string missedMatchesFileAddress, string kmerFrequLocFileAddress, string queryFileAddress, string readFileAddress, uint32_t queryCount, uint32_t readCount);
    CheckKmerFrequency<uint64_t> ckf(15, 50, 30);
    ckf.checkKmerFreq(argv[1], argv[2], argv[3], argv[4], stod(argv[5]), stod(argv[6]));
}

void testAligner(int argc, char *argv[]){
    Alignment aln;
    aln.query = argv[1];
    aln.read = argv[2];
    // PostFilter pf(50, 2, 150, 30, 0.9);
    Aligner <uint32_t> aligner(stod(argv[3]), 2, 150, stod(argv[4]), stod(argv[5]), stod(argv[6]), false, 1);
    vector<Penalty> penalties = readPenalties(argv[7]);
    aln = aligner.alignDifferentPenaltyScores(argv[1], argv[2], 1, 1, 1, penalties);
    cout << "aln.partialMatchSize: " << aln.partialMatchSize << endl;
    cout << "aln.cigar: " << aln.cigar << endl;
    cout << "aln.queryRegionStartPos: " << aln.queryRegionStartPos << endl;
    cout << "aln.queryRegionEndPos: " << aln.queryRegionEndPos << endl;
    cout << "aln.readRegionStartPos: " << aln.readRegionStartPos << endl;
    cout << "aln.readRegionEndPos: " << aln.readRegionEndPos << endl;
    cout << "aln.criteriaCode: " << aln.criteriaCode << endl;
    cout << "aln.editDistance: " << aln.editDistance << endl;
    cout << "aln.substitutions: " << aln.substitutions << endl;
    cout << "aln.inDels: " << aln.inDels << endl;
}

// Wrapper function
void expensiveKmerFNEval(int argc, char* argv[]) {
    if (argc != 9) {
        cerr << "Usage: " << argv[0] << " parmikAlnFileAddress queryCount readsCount expKmersFNFileAddress offlineExpensiveIndexAddress readDatabaseAddress k outputDir" << endl;
        exit(EXIT_FAILURE);
    }

    string parmikAlnFileAddress = argv[1];
    uint32_t queryCount = static_cast<uint32_t>(stoul(argv[2])); // Convert string to unsigned long and then to uint32_t
    uint32_t readsCount = static_cast<uint32_t>(stoul(argv[3])); // Convert string to unsigned long and then to uint32_t
    string expKmersFNFileAddress = argv[4];
    string offlineExpensiveIndexAddress = argv[5];
    string readDatabaseAddress = argv[6];
    uint16_t k = static_cast<uint16_t>(stoul(argv[7])); // Convert string to unsigned long and then to uint16_t
    string outputDir = argv[8];

    ExpensiveKmerFNEvalutor eke;
    // Call the original function with the converted arguments
    eke.expensiveKmerFNEvalutor(parmikAlnFileAddress, queryCount, readsCount, expKmersFNFileAddress, offlineExpensiveIndexAddress, readDatabaseAddress, k, outputDir);
}

void evaluateSecondChance(int argc, char *argv[]){
    if (argc != 6) {
        cerr << "Usage: " << argv[0] << "parmikAlnFileAddressWithSC parmikAlnFileAddressWithoutSC queryCount outputDir queryFileAddress" << endl;
        exit(EXIT_FAILURE);
    }
    string parmikAlnFileAddressWithSC = argv[1];
    string parmikAlnFileAddressWithoutSC = argv[2];
    uint32_t queryCount = static_cast<uint32_t>(stoul(argv[3]));
    string outputDir = argv[4];
    string queryFileAddress = argv[5];
    EvaluateSecondChance esc;
    esc.evaluateSecondChance(parmikAlnFileAddressWithSC, parmikAlnFileAddressWithoutSC, queryCount, outputDir, queryFileAddress);
}

int main(int argc, char *argv[])
{
    // testCheckBlastEditPositionsWrapper(argc, argv);
    // checkParmikFNalignments(argc, argv);
    // expensiveKmerFNEval(argc, argv);
    // evaluateSecondChance(argc, argv);
    if(DEBUG_MODE) testAligner(argc, argv);
    if(EXE_MODE) run(argc, argv);
    return 0;
}