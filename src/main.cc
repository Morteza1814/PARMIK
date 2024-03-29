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
#include <cmath>

#define PARMIK_MODE_INDEX   0
#define PARMIK_MODE_ALIGN   1
#define PARMIK_MODE_COMPARE 2

int argParse(int argc, char** argv, Config &cfg){
	args::ArgumentParser parser("=========================Arguments===========================", "======================================================");
    args::ValueFlag<int> parmikModeArg(parser, "", "PARMIK mode",                               {'a', "mode"});
    args::ValueFlag<string> otherToolAddressArg(parser, "", "Other tool output",                {'b', "tool"});
	args::ValueFlag<int> contigSizeArg(parser, "", "Contig Size",                               {'c', "contigSize"});
    args::ValueFlag<int> identityPercentageArg(parser, "", "Identity Percentage",               {'d', "identityPercentage"});
	args::ValueFlag<int> editDistanceArg(parser, "", "Max Edit Distance (i/d/s)",               {'e', "editDistance"});
	args::ValueFlag<string> offlineIndexAddressArg(parser, "", "Offline Index Address",         {'f', "offlineIndex"});
    args::HelpFlag help(parser, "help", "Help",                                                 {'h', "help"});
	args::ValueFlag<int> readsCountArg(parser, "", "Number of Reads",                           {'i', "readCount"});
	args::ValueFlag<int> queryCountArg(parser, "", "Number of Queries",                         {'j', "queryCount"});
	args::ValueFlag<int> kmerLengthArg(parser, "", "Kmer Length",                               {'k', "kmerLen"});
    args::ValueFlag<string> otherToolArg(parser, "", "The Other Tool  (bwa, blast, etc)",       {'l', "otherTool"});
    args::ValueFlag<int> minExactMatchLenArg(parser, "", "Minimum Exact Match Length",          {'m', "minExactMatchLen"});
    args::ValueFlag<string> kmerRangesFileAddressArg(parser, "", "k-mer Ranges File Address",   {'n', "kmerRangesFileAddress"});
	args::ValueFlag<string> outputDirArg(parser, "", "OutputDir",                               {'o', "outputDir"});
    args::ValueFlag<string> penaltyFileAddressArg(parser, "", "Penalty File Address",           {'p', "penaltyFileAddress"});
	args::ValueFlag<string> queryFileAddressArg(parser, "", "Query File Address",               {'q', "query"});
	args::ValueFlag<string> readDatabaseAddressArg(parser, "", "Read Data Base Address",        {'r', "read"});
	args::ValueFlag<int> regionSizeArg(parser, "", "Region Size",                               {'s', "regionSize"});
    args::ValueFlag<int> cheapKmerThresholdArg(parser, "", "Cheap Kmer Threshold",              {'t', "cheapKmerThreshold"});
	args::Flag isVerboseLogArg(parser, "", "Verbose Logging",                                   {'v', "verboseLog"});
	args::Flag isIndexOfflineArg(parser, "", "Is the read index offline",                       {'x', "isIndexOffline"});
    // args::ValueFlag<int> overlapSizeArg(parser, "", "Overlap Size",                         {'z', "overlapSize"});
	// args::Flag isFreqAndMemReportArg(parser, "", "Report Seed Frequencies and Avg MEM sizes", {'z', "memfreq"});
	// args::Flag isReverseStrandArg(parser, "", "Check the Reverse Strand of Queries", {'y', "reverseStrand"});
	// args::Flag isOnlyBestAlignmentArg(parser, "", "Only Output the Largest Alignment for Each Query", {'b', "bestAlign"});
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
	if (readDatabaseAddressArg) {cfg.readDatabaseAddress = args::get(readDatabaseAddressArg); } else {cout << "no readDatabaseAddress!"<< endl; return 0;}
	if (queryFileAddressArg) {cfg.queryFileAddress = args::get(queryFileAddressArg);} else {cout << "no queryFileAddress!"<< endl; return 0;}
    if (penaltyFileAddressArg) {cfg.penaltyFileAddress = args::get(penaltyFileAddressArg);} else {cfg.penaltyFileAddress = ""; cout << "no penaltyFileAddress!"<< endl;}
	if (outputDirArg) {cfg.outputDir = args::get(outputDirArg);} else {cout << "no outputDirArg!"<< endl; cfg.noOutputFileDump = true;}
    if (kmerRangesFileAddressArg) {cfg.kmerRangesFileAddress = args::get(kmerRangesFileAddressArg);} else {cout << "no kmerRangesFileAddressArg!"<< endl; cfg.kmerRangesFileAddress = "";}
    if (identityPercentageArg) {cfg.identityPercentage = args::get(identityPercentageArg);cfg.identityPercentage = (double)(cfg.identityPercentage/100);} else {cout << "no identityPercentage (default = 0.9)!"<< endl;}
	if (readsCountArg) {cfg.readsCount = args::get(readsCountArg); } else {cfg.readsCount = NUMBER_OF_READS;}
	if (queryCountArg) {cfg.queryCount = args::get(queryCountArg); } else {cfg.queryCount = NUMBER_OF_QUERIES;}
	if (kmerLengthArg) {cfg.kmerLength = args::get(kmerLengthArg); } else {cfg.kmerLength = KMER_SZ;}
	if (regionSizeArg) {cfg.regionSize = args::get(regionSizeArg); } else {cfg.regionSize = REGION_SZ;}
	if (contigSizeArg) {cfg.contigSize = args::get(contigSizeArg); } else {cfg.contigSize = CONTIG_SZ;}
    if (cheapKmerThresholdArg) {cfg.cheapKmerThreshold = args::get(cheapKmerThresholdArg); } else {cout << "no cheapKmerThreshold!"<< endl; return 0;}
    if (isIndexOfflineArg) {cfg.isIndexOffline = true; } else {cfg.isIndexOffline = false;}
    if (isVerboseLogArg) {cfg.isVerboseLog = true; } else {cfg.isVerboseLog = false;}
    if (offlineIndexAddressArg) {cfg.offlineIndexAddress = args::get(offlineIndexAddressArg); } else {cout << "no offlineIndexAddress!"<< endl; return 0;}
    if (otherToolAddressArg) {cfg.otherToolOutputFileAddress = args::get(otherToolAddressArg); } else {cout << "no otherToolOutputFileAddress!"<< endl; return 0;}
    if (parmikModeArg) {cfg.parmikMode = args::get(parmikModeArg); } else {cout << "no parmik mode is determined!"<< endl; return 0;}
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

void convertSamToAln(SamReader::Sam parmikSamAlignment, Alignment& parmikAlignment)
{
    SamReader sam("");
    parmikAlignment.readID = parmikSamAlignment.readId;
    parmikAlignment.queryID = parmikSamAlignment.queryId;
    parmikAlignment.cigar = parmikSamAlignment.cigar;
    parmikAlignment.substitutions = sam.countSubstitutions(parmikSamAlignment.cigar);
    parmikAlignment.matches = sam.countMatches(parmikSamAlignment.cigar);
    parmikAlignment.inDels = sam.countInsertions(parmikSamAlignment.cigar) + sam.countDeletions(parmikSamAlignment.cigar);
    parmikAlignment.editDistance  = parmikAlignment.substitutions + parmikAlignment.inDels;
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
    cout << left << setw(30) << "identityPercentage: " << cfg.identityPercentage << endl;
    cout << left << setw(30) << "minExactMatchLen: " << cfg.minExactMatchLen << endl;
    uint32_t minNumExactMatchKmer = 0;
    if(cfg.minExactMatchLen > 0) {
        minNumExactMatchKmer = cfg.minExactMatchLen - (cfg.kmerLength - 1);
    } else {
        cfg.editDistance = (uint32_t)(round(cfg.regionSize * (1-cfg.identityPercentage)));
        uint16_t exactMatchSize = cfg.regionSize - cfg.editDistance;
        minNumExactMatchKmer = (uint32_t)(floor((exactMatchSize/cfg.kmerLength) + (exactMatchSize%cfg.kmerLength)));
        cfg.editDistance = minNumExactMatchKmer - 1;
    }
    assert(cfg.minExactMatchLen == 0 && "min exact match len should be 0 for this version");
    assert(minNumExactMatchKmer > 0 && "minNumExactMatchKmer must be greater than 0 for cheap k-mer matching");
    cout << left << setw(30) << "minNumExactMatchKmer: " << minNumExactMatchKmer << endl;
    cout << left << setw(30) << "isIndexOffline: " << cfg.isIndexOffline << endl;
    cout << left << setw(30) << "offlineIndexAddress: " << cfg.offlineIndexAddress << endl;
    cout << left << setw(30) << "max editDistance in min region: " << cfg.editDistance << endl;
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
        if (cfg.parmikMode != PARMIK_MODE_INDEX)
        {
            queryCount = util.readContigsFromFile(cfg.queryFileAddress, cfg.queryCount, queries);
            cout << "queryCount : " << queryCount << endl;
            // read the penalties
        }
        uint32_t readCount = util.readContigsFromFile(cfg.readDatabaseAddress, cfg.readsCount, reads);
        cout << "readCount : " << readCount << endl;
        // unordered_map<uint32_t, unordered_set<LevAlign>> alignments;
        string offlineCheapIndexAddress = cfg.offlineIndexAddress + "ck_T" + to_string(cfg.cheapKmerThreshold) + "_K"+ to_string(cfg.kmerLength) + "_r" + to_string(cfg.readsCount);
        string parmikAlignmentsAddress = cfg.outputDir + "/aln/pmAln_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(minNumExactMatchKmer) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
        if (cfg.parmikMode != PARMIK_MODE_COMPARE)
        {
            if (cfg.kmerLength <= 16)
            {
                Container<uint32_t, uint32_t> cheapKmers;
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
                    kfc.collectCheapKmers(cheapKmers, invertedIndex, cfg.kmerRangesFileAddress);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // cheapKmers.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    invertedIndex.clear();
                    clock_t ser_start_time = clock();
                    cheapKmers.serialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Serialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                } else 
                {
                    clock_t ser_start_time = clock();
                    cheapKmers.deserialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Deserialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                }
                if (cfg.parmikMode == PARMIK_MODE_ALIGN)
                {
                    //do partial matching based on cheap k-mers
                    CheapKmerPartialMatcher<uint32_t, uint32_t, uint32_t> ckpm(cfg.kmerLength, cfg.contigSize, minNumExactMatchKmer, cfg.isVerboseLog);
                    ckpm.cheapSeedFilter(cheapKmers, queries, minThCheapSeedReads);
                    //get the reverse complement of queries
                    tsl::robin_map <uint32_t, string> revQueries = util.reverseComplementMapValues(queries);
                    ckpm.cheapSeedFilter(cheapKmers, revQueries, revMinThCheapSeedReads);
                    // ckpm.printArrays();
                    // SeedMatchExtender<uint32_t, uint64_t> pm(cfg.minExactMatchLen, cfg.regionSize, cfg.isVerboseLog, cfg.editDistance, cfg.contigSize, cfg.inDelPenalty, cfg.subPenalty);
                    Aligner <uint32_t> aligner(cfg.regionSize, cfg.editDistance, cfg.contigSize, cfg.kmerLength, minNumExactMatchKmer, cfg.identityPercentage);
                    // pm.findPartiaMatches(reads, queries, minThCheapSeedReads, backMinThCheapSeedReads, queryCount, pmr, true, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, queries, minThCheapSeedReads, queryCount, true, parmikAlignmentsAddress, penalties);
                    //do it again for the reverse strand
                    // pm.findPartiaMatches(reads, revQueries, revMinThCheapSeedReads, revBackMinThCheapSeedReads, queryCount, pmr, false, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, revQueries, revMinThCheapSeedReads, queryCount, false, parmikAlignmentsAddress, penalties);
                }
            } else
            {
                Container<uint64_t, uint32_t> cheapKmers;
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
                    kfc.collectCheapKmers(cheapKmers, invertedIndex, cfg.kmerRangesFileAddress);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // cheapKmers.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Cheap K-mer Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    invertedIndex.clear();
                    clock_t ser_start_time = clock();
                    cheapKmers.serialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Serialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                } else 
                {
                    clock_t ser_start_time = clock();
                    cheapKmers.deserialize(offlineCheapIndexAddress);
                    cout << left << setw(30) << fixed << setprecision(2) << "Cheap Kmers Deserialization Time: " << (double)(clock() - ser_start_time)/CLOCKS_PER_SEC << " seconds" << endl;
                }
                if (cfg.parmikMode == PARMIK_MODE_ALIGN)
                {
                    // do partial matching based on cheap k-mers
                    CheapKmerPartialMatcher<uint32_t, uint64_t, uint32_t> ckpm(cfg.kmerLength, cfg.contigSize, minNumExactMatchKmer, cfg.isVerboseLog);
                    ckpm.cheapSeedFilter(cheapKmers, queries, minThCheapSeedReads);
                    //get the reverse complement of queries
                    tsl::robin_map <uint32_t, string> revQueries = util.reverseComplementMapValues(queries);
                    ckpm.cheapSeedFilter(cheapKmers, revQueries, revMinThCheapSeedReads);
                    // ckpm.printArrays();
                    // SeedMatchExtender<uint32_t, uint64_t> pm(cfg.minExactMatchLen, cfg.regionSize, cfg.isVerboseLog, cfg.editDistance, cfg.contigSize, cfg.inDelPenalty, cfg.subPenalty);
                    Aligner <uint32_t> aligner(cfg.regionSize, cfg.editDistance, cfg.contigSize, cfg.kmerLength, minNumExactMatchKmer, cfg.identityPercentage);
                    // pm.findPartiaMatches(reads, queries, minThCheapSeedReads, backMinThCheapSeedReads, queryCount, pmr, true, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, queries, minThCheapSeedReads, queryCount, true, parmikAlignmentsAddress, penalties);
                    //do it again for the reverse strand
                    // pm.findPartiaMatches(reads, revQueries, revMinThCheapSeedReads, revBackMinThCheapSeedReads, queryCount, pmr, false, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, revQueries, revMinThCheapSeedReads, queryCount, false, parmikAlignmentsAddress, penalties);
                }
            }
        
            //report the histogram of the alignments
            // if(cfg.parmikMode == PARMIK_MODE_ALIGN && !cfg.noOutputFileDump)
            // {
            //     ofstream matches(cfg.outputDir + "/aln/PM_best_matches.txt");
            //     //report PM final best match results
            //     matches << "Q , R" << endl;
            //     for (const auto& pair : pmr) {
            //         matches << pair.first << " , " << pair.second.readID << endl;
            //     }
            //     matches.close();
            // }
        } else
        {
            if(cfg.noOutputFileDump)
            {
                cout << "no noOutputFileDump!!" << endl;
                return 0;
            }
            //load the parmik alignments from the sam formatted file
            SamReader parmikSam(parmikAlignmentsAddress);
            vector<SamReader::Sam> parmikSamAlignments = parmikSam.parseFile(queryCount);
            // IndexContainer<uint32_t, LevAlign> parmikMultiAlignments;
            IndexContainer<uint32_t, Alignment> parmikMultiAlignments;
            for (const SamReader::Sam& aln : parmikSamAlignments) 
            {
                // LevAlign l;
                Alignment l;
                // convertSamToLev(aln, l);
                convertSamToAln(aln, l);
                parmikMultiAlignments.put(aln.queryId, l);
            }
            // cout<<"finished \n";
            //check the alignment results with another aligner
            string comparisonResultsFileAddress = cfg.outputDir + "/cmp/" + cfg.otherTool + "/cmp_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(minNumExactMatchKmer) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
            string alnPerQueryFileAddress = cfg.outputDir + "/cmp/" + cfg.otherTool + "/AlnPerQ_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(minNumExactMatchKmer) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
            string parmikFnReadsFileAddress = cfg.outputDir + "/cmp/" + cfg.otherTool + "/ParmikFnReads_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(minNumExactMatchKmer) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
            string bestAlnCmpFileAddress = cfg.outputDir + "/cmp/" + cfg.otherTool + "/BestAlnCmp_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(minNumExactMatchKmer) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
            vector<std::pair<uint32_t, uint32_t>> alnPmVsOtherAlnSizesMap;
            // (BWA)
            if(cfg.otherTool == "BWA" || cfg.otherTool == "bwa")
            {
                // ComparatorWithBWA cwb;
                // cwb.comparePmWithBWA(cfg, reads, queries, comparisonResultsFileAddress, parmikMultiAlignments, alnPmVsOtherAlnSizesMap, queryCount, alnPerQueryFileAddress);
            } else if(cfg.otherTool == "BLAST" || cfg.otherTool == "blast")
            {
                CompareWithBlast cwb;
                cwb.comparePmWithBlast(cfg, reads, queries, comparisonResultsFileAddress, parmikMultiAlignments, alnPmVsOtherAlnSizesMap, queryCount, minNumExactMatchKmer, alnPerQueryFileAddress, parmikFnReadsFileAddress, bestAlnCmpFileAddress);
            } else if(cfg.otherTool == "GT" || cfg.otherTool == "gt")
            {
                SamReader gtSam(cfg.otherToolOutputFileAddress);
                vector<SamReader::Sam> gtSamAlignments = gtSam.parseFile(queryCount);
                IndexContainer<uint32_t, LevAlign> gtMultiAlignments;
                for (const SamReader::Sam& aln : gtSamAlignments) 
                {
                    LevAlign l;
                    convertSamToLev(aln, l);
                    gtMultiAlignments.put(aln.queryId, l);
                }
                // CompareWithGroundTruth cwgt;
                // cwgt.comparePmWithGroundTruth(gtMultiAlignments, cfg, reads, queries, comparisonResultsFileAddress, parmikMultiAlignments, alnPmVsOtherAlnSizesMap, queryCount, alnPerQueryFileAddress);
            
            }
         
            // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<BWA vs PM Alignments Sizes started!>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
            // ofstream alnPmVsBwaAlnSizes(cfg.outputDir + "alnPmVsBwaAlnSizes.txt");
            // for (const auto& pair : alnPmVsOtherAlnSizesMap) 
            // {
            //     alnPmVsBwaAlnSizes << pair.first << " " << pair.second << endl;
            // }
            // alnPmVsBwaAlnSizes.close();
            // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<BWA vs PM Alignments Sizes finished!>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        }

    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << endl;
        return 1;
    }

    return 0;
}

// void testCheckBlastEditPositionsWrapper(int argc, char *argv[]){
//     CompareWithBlast cwb;
//     // testCheckBlastEditPositions(uint32_t contigSize, uint32_t regionSize, uint32_t editDistance, uint32_t minExactMatchLen, uint32_t queryS, string qAln, sting readAln, uint32_t blastMismatches, uint32_t blastInDel);
//     cwb.testCheckBlastEditPositions(stod(argv[1]), stod(argv[2]), stod(argv[3]), stod(argv[4]), stod(argv[5]), argv[6], argv[7], stod(argv[8]), stod(argv[9]), stod(argv[10]));
// }

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
    Aligner <uint32_t> aligner(stod(argv[3]), 2, 150, stod(argv[4]), stod(argv[5]), stod(argv[6]));
    // aligner.align(aln, stod(argv[3]), stod(argv[4]), stod(argv[5]), stod(argv[6]));
    // // bool criteriaCheck = aligner.checkAlingmentCriteria(aln);
    // bool criteriaCheck = pf.checkAndUpdateBasedOnAlingmentCriteria(aln);
    // cout << ((criteriaCheck == true) ? "aln meets criteria" : "aln not meet criteria") << endl;
    // vector<Penalty> penalties = readPenalties("/u/rgq5aw/GIT/PARMIK/experiments/parmik/CK_SSW_FLTR/PenaltySets/2284_1111_2444_2288_2848");
    vector<Penalty> penalties = readPenalties(argv[7]);
    aln = aligner.alignDifferentPenaltyScores(argv[1], argv[2], 1, 1, 1, penalties);
    cout << "aln.partialMatchSize: " << aln.partialMatchSize << endl;
    cout << "aln.cigar: " << aln.cigar << endl;
    cout << "aln.queryRegionStartPos: " << aln.queryRegionStartPos << endl;
    cout << "aln.queryRegionEndPos: " << aln.queryRegionEndPos << endl;
    cout << "aln.readRegionStartPos: " << aln.readRegionStartPos << endl;
    cout << "aln.readRegionEndPos: " << aln.readRegionEndPos << endl;
}

int main(int argc, char *argv[])
{
    // testCheckBlastEditPositionsWrapper(argc, argv);
    // checkParmikFNalignments(argc, argv);
    if(DEBUG_MODE) testAligner(argc, argv);
    if(EXE_MODE) run(argc, argv);
    return 0;
}