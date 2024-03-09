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

#define PARMIK_MODE_INDEX   0
#define PARMIK_MODE_ALIGN   1
#define PARMIK_MODE_COMPARE 2

int argParse(int argc, char** argv, Config &cfg){
	args::ArgumentParser parser("=========================Arguments===========================", "======================================================");
    args::ValueFlag<int> parmikModeArg(parser, "", "PARMIK mode",                              {'a', "mode"});
    args::ValueFlag<string> bwaSamAddressArg(parser, "", "Other tool output",              {'b', "tool"});
	args::ValueFlag<int> contigSizeArg(parser, "", "Contig Size",                           {'c', "contigSize"});
	args::ValueFlag<int> editDistanceArg(parser, "", "Edit Distance (i/d/s)",               {'e', "editDistance"});
	args::ValueFlag<string> offlineIndexAddressArg(parser, "", "Offline Index Address",     {'f', "offlineIndex"});
    args::HelpFlag help(parser, "help", "Help",                                             {'h', "help"});
	args::ValueFlag<int> readsCountArg(parser, "", "Number of Reads",                       {'i', "readCount"});
	args::ValueFlag<int> queryCountArg(parser, "", "Number of Queries",                     {'j', "queryCount"});
	args::ValueFlag<int> kmerLengthArg(parser, "", "Kmer Length",                           {'k', "kmerLen"});
    args::ValueFlag<string> otherToolArg(parser, "", "The Other Tool  (bwa, blast, etc)",   {'l', "otherTool"});
	args::ValueFlag<int> minExactMatchLenArg(parser, "", "Minimum Exact Match Len",         {'m', "minExactMatchLen"});
	args::ValueFlag<string> outputDirArg(parser, "", "OutputDir",                           {'o', "outputDir"});
    args::ValueFlag<string> penaltyFileAddressArg(parser, "", "Penalty File Address",       {'p', "penaltyFileAddress"});
	args::ValueFlag<string> queryFileAddressArg(parser, "", "Query File Address",           {'q', "query"});
	args::ValueFlag<string> readDatabaseAddressArg(parser, "", "Read Data Base Address",    {'r', "read"});
	args::ValueFlag<int> regionSizeArg(parser, "", "Region Size",                           {'s', "regionSize"});
    args::ValueFlag<int> cheapKmerThresholdArg(parser, "", "Cheap Kmer Threshold",          {'t', "cheapKmerThreshold"});
	args::Flag isVerboseLogArg(parser, "", "Verbose Logging",                               {'v', "verboseLog"});
	args::Flag isIndexOfflineArg(parser, "", "Is the read inndex offline",                  {'x', "isIndexOffline"});
    args::ValueFlag<int> overlapSizeArg(parser, "", "Overlap Size",                         {'z', "overlapSize"});
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
	if (readsCountArg) {cfg.readsCount = args::get(readsCountArg); } else {cfg.readsCount = NUMBER_OF_READS;}
	if (queryCountArg) {cfg.queryCount = args::get(queryCountArg); } else {cfg.queryCount = NUMBER_OF_QUERIES;}
	if (kmerLengthArg) {cfg.kmerLength = args::get(kmerLengthArg); } else {cfg.kmerLength = KMER_SZ;}
	if (minExactMatchLenArg) {cfg.minExactMatchLen = args::get(minExactMatchLenArg); } else {cfg.minExactMatchLen = MIN_EXACT_MATCH_LEN;}
	if (regionSizeArg) {cfg.regionSize = args::get(regionSizeArg); } else {cfg.regionSize = REGION_SZ;}
	if (contigSizeArg) {cfg.contigSize = args::get(contigSizeArg); } else {cfg.contigSize = CONTIG_SZ;}
    if (cheapKmerThresholdArg) {cfg.cheapKmerThreshold = args::get(cheapKmerThresholdArg); } else {cout << "no cheapKmerThreshold!"<< endl; return 0;}
    if (overlapSizeArg) {cfg.overlapSize = args::get(overlapSizeArg); } else {cout << "no overlapSize!"<< endl; return 0;}
    if (isIndexOfflineArg) {cfg.isIndexOffline = true; } else {cfg.isIndexOffline = false;}
    if (isVerboseLogArg) {cfg.isVerboseLog = true; } else {cfg.isVerboseLog = false;}
    if (offlineIndexAddressArg) {cfg.offlineIndexAddress = args::get(offlineIndexAddressArg); } else {cout << "no offlineIndexAddress!"<< endl; return 0;}
    if (bwaSamAddressArg) {cfg.otherToolOutputFileAddress = args::get(bwaSamAddressArg); } else {cout << "no otherToolOutputFileAddress!"<< endl; return 0;}
    if (parmikModeArg) {cfg.parmikMode = args::get(parmikModeArg); } else {cout << "no parmik mode is determined!"<< endl; return 0;}
    if (otherToolArg) {cfg.otherTool = args::get(otherToolArg);} else {cout << "no otherToolArg!"<< endl;if(cfg.parmikMode == PARMIK_MODE_COMPARE) return 0;}
	if (editDistanceArg) {cfg.editDistance = args::get(editDistanceArg); } else {cfg.editDistance = NUMBER_OF_ALLOWED_EDIT_DISTANCES;}
	cfg.readFileName = cfg.readDatabaseAddress.substr(cfg.readDatabaseAddress.find_last_of("/\\") + 1);
	cfg.queryFileName = cfg.queryFileAddress.substr(cfg.queryFileAddress.find_last_of("/\\") + 1);
	// cfg.numberOfKmers = (cfg.contigSize - cfg.kmerLength + 1);
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

void testOnePair(int argc, char *argv[])
{
    string read =  argv[1];
    string query = argv[2];
    double idPen = stod(argv[3]);
    double subPen = stod(argv[4]);
    SeedMatchExtender<uint32_t, uint64_t> pm(30, 50, true, 2, 150, idPen, subPen);
    cout << "read  : " << read << endl;
    cout << "query : " << query << endl;
    vector<uint64_t> frontKmers;
    char frontRegion = 'F';
    pm.extractKmersFromRegions(query, frontKmers, frontRegion);
    LevAlign fla = pm.extendSeed(query, read, frontKmers, frontRegion);
    cout << "front partialMatchSize : " << fla.partialMatchSize << ", editDistance : " << fla.editDistance << endl;
    cout << "or R : " << fla.read << endl;
    cout << "or Q : " << fla.query << endl;
    cout << "al R : " << fla.alignedRead << endl;
    cout << "al Q : " << fla.alignedQuery << endl;
    cout << "ed T : " << fla.editDistanceTypes << endl;

    vector<uint64_t> backKmers;
    char backRegion = 'B';
    pm.extractKmersFromRegions(query, backKmers, backRegion);
    LevAlign bla = pm.extendSeed(query, read, backKmers, backRegion);
    cout << "back partialMatchSize : " << bla.partialMatchSize << ", editDistance : " << bla.editDistance << endl;
    cout << "or R : " << bla.read << endl;
    cout << "or Q : " << bla.query << endl;
    cout << "al R : " << bla.alignedRead << endl;
    cout << "al Q : " << bla.alignedQuery << endl;
    cout << "ed T : " << bla.editDistanceTypes << endl;
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

vector<Penalty> readPenalties(string& readPenaltiesFileAddress)
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
	cout << left << setw(30) << "readsCount: " << cfg.readsCount << endl; 
	cout << left << setw(30) << "queryCount: " << cfg.queryCount << endl;
	cout << left << setw(30) << "kmerLength: " << cfg.kmerLength << endl;
	cout << left << setw(30) << "regionSize: " << cfg.regionSize << endl;
	cout << left << setw(30) << "minExactMatchLen: " << cfg.minExactMatchLen << endl;
    cout << left << setw(30) << "overlapsize: " << cfg.overlapSize << endl;
    cout << left << setw(30) << "cheapKmerThreshold: " << cfg.cheapKmerThreshold << endl;
    uint32_t minNumExactMatchKmer = ((cfg.minExactMatchLen-cfg.overlapSize)/(cfg.kmerLength - cfg.overlapSize));
    cout << left << setw(30) << "minNumExactMatchKmer: " << minNumExactMatchKmer << endl;
    cout << left << setw(30) << "isIndexOffline: " << cfg.isIndexOffline << endl;
    cout << left << setw(30) << "offlineIndexAddress: " << cfg.offlineIndexAddress << endl;
    cout << left << setw(30) << "isVeboseLog: " << cfg.isVerboseLog << endl;
	cout << "**********************CONFIGURATIONS*****************************" << endl;

    try { 
        //construct the baseline
        // IndexContainer<uint32_t, uint32_t> baselineQueriesFrontSeeds; 
        // IndexContainer<uint32_t, uint32_t> baselineQueriesBackSeeds;
        // BaseLinePartialMatcher<uint32_t, uint64_t, uint32_t> blpm(cfg.minExactMatchLen, cfg.regionSize);
        // blpm.constructBaseline(cfg.readDatabaseAddress, cfg.readsCount, cfg.queryFileAddress, cfg.queryCount,cfg.minExactMatchLen, cfg.cheapKmerThreshold, baselineQueriesFrontSeeds, baselineQueriesBackSeeds, cfg.isIndexOffline, cfg.offlineIndexAddress);
        // run the experiment
        IndexContainer<uint32_t, uint32_t> frontMinThCheapSeedReads, backMinThCheapSeedReads, revFrontMinThCheapSeedReads, revBackMinThCheapSeedReads;
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
        string parmikAlignmentsAddress = cfg.outputDir + "/aln/pmAln_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(cfg.minExactMatchLen) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
        if (cfg.parmikMode != PARMIK_MODE_COMPARE)
        {
            if (cfg.kmerLength <= 16)
            {
                Container<uint32_t, uint32_t> cheapKmers;
                if(!cfg.isIndexOffline)
                {
                    //create the partial matching inverted read index
                    InvertedIndexBuilder<uint32_t, uint32_t> builder(cfg.kmerLength, cfg.overlapSize);
                    // Build the inverted index
                    IndexContainer<uint32_t, uint32_t> invertedIndex = builder.build(cfg.readDatabaseAddress, cfg.readsCount, cfg.isIndexOffline);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // invertedIndex.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    //collect cheap kmers
                    KmersFrequencyCounter<uint32_t, uint32_t> kfc(cfg.cheapKmerThreshold);
                    kfc.collectCheapKmers(cheapKmers, invertedIndex);
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
                    CheapKmerPartialMatcher<uint32_t, uint32_t, uint32_t> ckpm50(cfg.kmerLength, cfg.regionSize, minNumExactMatchKmer, cfg.isVerboseLog);
                    ckpm50.cheapSeedFilter(cheapKmers, queries, frontMinThCheapSeedReads, backMinThCheapSeedReads);
                    //get the reverse complement of queries
                    tsl::robin_map <uint32_t, string> revQueries = util.reverseComplementMapValues(queries);
                    ckpm50.cheapSeedFilter(cheapKmers, revQueries, revFrontMinThCheapSeedReads, revBackMinThCheapSeedReads);
                    // ckpm50.printArrays();
                    // SeedMatchExtender<uint32_t, uint64_t> pm(cfg.minExactMatchLen, cfg.regionSize, cfg.isVerboseLog, cfg.editDistance, cfg.contigSize, cfg.inDelPenalty, cfg.subPenalty);
                    Aligner <uint32_t> aligner(cfg.regionSize, cfg.editDistance, cfg.contigSize, cfg.minExactMatchLen);
                    // pm.findPartiaMatches(reads, queries, frontMinThCheapSeedReads, backMinThCheapSeedReads, queryCount, pmr, true, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, queries, frontMinThCheapSeedReads, backMinThCheapSeedReads, queryCount, true, parmikAlignmentsAddress, penalties);
                    //do it again for the reverse strand
                    // pm.findPartiaMatches(reads, revQueries, revFrontMinThCheapSeedReads, revBackMinThCheapSeedReads, queryCount, pmr, false, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, revQueries, revFrontMinThCheapSeedReads, revBackMinThCheapSeedReads, queryCount, false, parmikAlignmentsAddress, penalties);
                }
            } else
            {
                Container<uint64_t, uint32_t> cheapKmers;
                if(!cfg.isIndexOffline)
                {
                    //create the partial matching inverted read index
                    InvertedIndexBuilder<uint64_t, uint32_t> builder(cfg.kmerLength, cfg.overlapSize);
                    // Build the inverted index
                    IndexContainer<uint64_t, uint32_t> invertedIndex = builder.build(cfg.readDatabaseAddress, cfg.readsCount, cfg.isIndexOffline);
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    // invertedIndex.getSize();
                    // cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Inverted Index Size>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                    //collect cheap kmers
                    KmersFrequencyCounter<uint64_t, uint32_t> kfc(cfg.cheapKmerThreshold);
                    kfc.collectCheapKmers(cheapKmers, invertedIndex);
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
                    CheapKmerPartialMatcher<uint32_t, uint64_t, uint32_t> ckpm50(cfg.kmerLength, cfg.regionSize, minNumExactMatchKmer, cfg.isVerboseLog);
                    ckpm50.cheapSeedFilter(cheapKmers, queries, frontMinThCheapSeedReads, backMinThCheapSeedReads);
                    //get the reverse complement of queries
                    tsl::robin_map <uint32_t, string> revQueries = util.reverseComplementMapValues(queries);
                    ckpm50.cheapSeedFilter(cheapKmers, revQueries, revFrontMinThCheapSeedReads, revBackMinThCheapSeedReads);
                    // ckpm50.printArrays();
                    // SeedMatchExtender<uint32_t, uint64_t> pm(cfg.minExactMatchLen, cfg.regionSize, cfg.isVerboseLog, cfg.editDistance, cfg.contigSize, cfg.inDelPenalty, cfg.subPenalty);
                    Aligner <uint32_t> aligner(cfg.regionSize, cfg.editDistance, cfg.contigSize, cfg.minExactMatchLen);
                    // pm.findPartiaMatches(reads, queries, frontMinThCheapSeedReads, backMinThCheapSeedReads, queryCount, pmr, true, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, queries, frontMinThCheapSeedReads, backMinThCheapSeedReads, queryCount, true, parmikAlignmentsAddress, penalties);
                    //do it again for the reverse strand
                    // pm.findPartiaMatches(reads, revQueries, revFrontMinThCheapSeedReads, revBackMinThCheapSeedReads, queryCount, pmr, false, parmikAlignmentsAddress);
                    aligner.findPartiaMatches(reads, revQueries, revFrontMinThCheapSeedReads, revBackMinThCheapSeedReads, queryCount, false, parmikAlignmentsAddress, penalties);
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
            string comparisonResultsFileAddress = cfg.outputDir + "/cmp/" + cfg.otherTool + "/cmp_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(cfg.minExactMatchLen) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
            string alnPerQueryFileAddress = cfg.outputDir + "/cmp/" + cfg.otherTool + "/AlnPerQ_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(cfg.minExactMatchLen) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
            string parmikFnReadsFileAddress = cfg.outputDir + "/cmp/" + cfg.otherTool + "/ParmikFnReads_" + "R" + to_string(cfg.regionSize) + "_M" + to_string(cfg.minExactMatchLen) + "_E" + to_string(cfg.editDistance) + "_K" + to_string(cfg.kmerLength) + "_T" + to_string(cfg.cheapKmerThreshold) + "_P" + getPenaltiesSubstr(penalties) + ".txt";
            vector<std::pair<uint32_t, uint32_t>> alnPmVsOtherAlnSizesMap;
            // (BWA)
            if(cfg.otherTool == "BWA" || cfg.otherTool == "bwa")
            {
                // ComparatorWithBWA cwb;
                // cwb.comparePmWithBWA(cfg, reads, queries, comparisonResultsFileAddress, parmikMultiAlignments, alnPmVsOtherAlnSizesMap, queryCount, alnPerQueryFileAddress);
            } else if(cfg.otherTool == "BLAST" || cfg.otherTool == "blast")
            {
                CompareWithBlast cwb;
                cwb.comparePmWithBlast(cfg, reads, queries, comparisonResultsFileAddress, parmikMultiAlignments, alnPmVsOtherAlnSizesMap, queryCount, alnPerQueryFileAddress, parmikFnReadsFileAddress);
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

        // if(cfg.isVerboseLog)
        // {
        //     //check FN, FP,...
        //     cout << "<<<<<<<<<<<<<<<<<<<<<<<<Partial Matching Seeding FN results>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        //     map<uint32_t, uint32_t> frontSeedsFNRates;
        //     map<uint32_t, uint32_t> backSeedsFNRates;
        //     map<uint32_t, uint32_t> frontSeedsTPRates = checkRates(baselineQueriesFrontSeeds, frontMinThCheapSeedReads, frontSeedsFNRates);
        //     map<uint32_t, uint32_t> backSeedsTPRates = checkRates(baselineQueriesBackSeeds, backMinThCheapSeedReads, backSeedsFNRates);
        //     tuple<uint32_t, uint32_t, uint32_t> frontSeedsFNRateTuple = util.calculateStatistics(convertMapToSet(frontSeedsFNRates));
        //     tuple<uint32_t, uint32_t, uint32_t> backSeedsFNRateTuple = util.calculateStatistics(convertMapToSet(backSeedsFNRates));
        //     uint32_t frontSeedsFNRate = sumMapValues(frontSeedsFNRates);
        //     uint32_t backSeedsFNRate= sumMapValues(backSeedsFNRates);
        //     printf("Number of Seeds FN Rate  => Front [total: %d, average: %d, median: %d, mean: %d] - Back [total: %d, average: %d, median: %d, mean: %d]\n", frontSeedsFNRate, get<0>(frontSeedsFNRateTuple), get<1>(frontSeedsFNRateTuple), get<2>(frontSeedsFNRateTuple), backSeedsFNRate, get<0>(backSeedsFNRateTuple), get<1>(backSeedsFNRateTuple), get<2>(backSeedsFNRateTuple));
        //     uint32_t frontSeedsTPRate = sumMapValues(frontSeedsTPRates);
        //     uint32_t backSeedsTPRate= sumMapValues(backSeedsTPRates);
        //     tuple<uint32_t, uint32_t, uint32_t> frontSeedsTPRatesTuple = util.calculateStatistics(convertMapToSet(frontSeedsTPRates));
        //     tuple<uint32_t, uint32_t, uint32_t> backSeedsTPRatesTuple = util.calculateStatistics(convertMapToSet(backSeedsTPRates));
        //     printf("Number of Seeds TP Rate  => Front [total: %d, average: %d, median: %d, mean: %d] - Back [total: %d, average: %d, median: %d, mean: %d]\n", frontSeedsTPRate, get<0>(frontSeedsTPRatesTuple), get<1>(frontSeedsTPRatesTuple), get<2>(frontSeedsTPRatesTuple), backSeedsTPRate, get<0>(backSeedsTPRatesTuple), get<1>(backSeedsTPRatesTuple), get<2>(backSeedsTPRatesTuple));
        //     cout << "<<<<<<<<<<<<<<<<<<<<<<<<Partial Matching Seeding FP results>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        //     map<uint32_t, uint32_t> frontSeedsFPRates;
        //     map<uint32_t, uint32_t> backSeedsFPRates;
        //     checkRates(frontMinThCheapSeedReads, baselineQueriesFrontSeeds, frontSeedsFPRates);
        //     checkRates(backMinThCheapSeedReads, baselineQueriesBackSeeds, backSeedsFPRates);
        //     uint32_t frontSeedsFPRate = sumMapValues(frontSeedsFPRates);
        //     uint32_t backSeedsFPRate= sumMapValues(backSeedsFPRates);
        //     tuple<uint32_t, uint32_t, uint32_t> frontSeedsFPRateTuple = util.calculateStatistics(convertMapToSet(frontSeedsFPRates));
        //     tuple<uint32_t, uint32_t, uint32_t> backSeedsFPRateTuple = util.calculateStatistics(convertMapToSet(backSeedsFPRates));
        //     printf("Number of Seeds FP Rate  => Front [total: %d, average: %d, median: %d, mean: %d] - Back [total: %d, average: %d, median: %d, mean: %d]\n", frontSeedsFPRate , get<0>(frontSeedsFPRateTuple), get<1>(frontSeedsFPRateTuple), get<2>(frontSeedsFPRateTuple), backSeedsFPRate, get<0>(backSeedsFPRateTuple), get<1>(backSeedsFPRateTuple), get<2>(backSeedsFPRateTuple));      
        // }
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
    PostFilter pf(50, 2, 150, 30);
    Aligner <uint32_t> aligner(50, 2, 150, 30);
    aligner.align(aln, stod(argv[3]), stod(argv[4]), stod(argv[5]), stod(argv[6]));
    // bool criteriaCheck = aligner.checkAlingmentCriteria(aln);
    bool criteriaCheck = pf.checkAlingmentCriteria(aln.editDistance, aln.partialMatchSize, aln.queryRegionStartPos, aligner.convertCigarToStr(aln.cigar), "cigarStr", aln.substitutions, aln.inDels, aln.criteriaCode);
    cout << ((criteriaCheck == true) ? "aln meets criteria" : "aln not meet criteria") << endl;
    cout << "criteriaCode: " << aln.criteriaCode << endl;
    if (criteriaCheck || (!criteriaCheck && (aln.criteriaCode <= 3)))
        cout << "dump it!!!" << endl;

}

// void testHasMinConsecutiveMatches(int argc, char *argv[]){
//     CompareWithBlast cwb;
//     Config cfg;
//     cfg.minExactMatchLen = stod(argv[3]);
//     cfg.regionSize = stod(argv[4]);
//     cfg.contigSize = stod(argv[5]);
//     bool hasMinConsecutiveMatches = cwb.hasMinConsecutiveMatches(stod(argv[6]), argv[1], argv[2], cfg);
//     cout << ((hasMinConsecutiveMatches == true) ? "hasMinConsecutiveMatches" : "no hasMinConsecutiveMatches") << endl;
// }

int main(int argc, char *argv[])
{
    // testCheckBlastEditPositionsWrapper(argc, argv);
    // testOnePair(argc, argv);
    // checkParmikFNalignments(argc, argv);
    // testAligner(argc, argv);
    // testHasMinConsecutiveMatches(argc, argv);
    run(argc, argv);
    return 0;
}