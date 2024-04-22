#ifndef SSW_BASELINE_H
#define SSW_BASELINE_H

class SSW_BaseLine {
private:
    uint16_t regionSize; // min region size to be considered for the alignement
    double identityPercentage;
public:
    SSW_BaseLine(uint16_t R, double i) : regionSize(R), identityPercentage(i) {}

    bool checkIdentityPercentange(string cigarStr){
        uint32_t matches = getMatchesCount(cigarStr);
        double identity =  (double) matches / (double) cigarStr.size();
        if(DEBUG_MODE) cout << "matches: " << matches << ", cigarStr.size : " << cigarStr.size() << ", identity: " << identity << endl;
        if(identity < identityPercentage)
            return false;
        return true;
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

    void smithWatermanAligner(Alignment &aln, uint16_t matchPen, uint16_t subPen, uint16_t gapoPen, uint16_t gapextPen)
    {
        Utilities<uint32_t> util;
        int32_t maskLen = strlen(aln.query.c_str())/2;
        maskLen = maskLen < 15 ? 15 : maskLen;

        // Declares a default Aligner
        StripedSmithWaterman::Aligner aligner(matchPen, subPen, gapoPen, gapextPen);
        // Declares a default filter
        StripedSmithWaterman::Filter filter;
        // Declares an alignment that stores the result
        StripedSmithWaterman::Alignment alignment;
        // Aligns the query to the ref
        aligner.Align(aln.query.c_str(), aln.read.c_str(), aln.read.size(), filter, &alignment, maskLen); 

        aln.readRegionStartPos = alignment.ref_begin;
        aln.readRegionEndPos = alignment.ref_end;
        aln.queryRegionStartPos = alignment.query_begin;
        aln.queryRegionEndPos = alignment.query_end;
        aln.cigar = alignment.cigar_string;
        aln.editDistance = alignment.mismatches;
        aln.score = alignment.sw_score;
        // cout << "aln.query_begin: " << alignment.query_begin << ", aln.query_end: " << alignment.query_end << ", aln.ref_begin: " << alignment.ref_begin << ", aln.ref_end: " << alignment.ref_end << ", cigar_string: " << alignment.cigar_string << endl;
        util.parseCigar(aln.cigar, aln.matches, aln.substitutions, aln.inDels, aln.editLocations);
        aln.partialMatchSize = aln.matches + aln.substitutions + aln.inDels;
        // for(auto it = aln.editLocations.begin(); it!= aln.editLocations.end(); it++){
        //     cout << *it << " ";
        // }
    }

    Alignment alignDifferentPenaltyScores(string query, string read, uint32_t queryID, uint32_t readID, bool isForwardStran, vector<Penalty> &penalties)
    {
        Alignment bestAlignment;
        if (penalties.size() == 0)
        {
            Alignment aln;
            aln.read = read;
            aln.query = query;
            aln.readID = readID;
            aln.queryID = queryID;
            if (!isForwardStran) aln.flag = 16;
            smithWatermanAligner(aln, 1, 1, 1, 1);
            string cigarStr = convertCigarToStr(aln.cigar);
            cigarStr = trimClips(cigarStr);
            if(cigarStr.length() >= regionSize && checkIdentityPercentange(cigarStr)) {
                if(DEBUG_MODE) cout << "cigar len was smaller than the R or K at the beginning\n";
                return aln;
            } else {
                //return an empty alignment
                return bestAlignment;
            }
        }
        vector<string> repeatedCigars;
        for (auto penalty : penalties)
        {
            if(DEBUG_MODE) 
            {
                cout << "<<<<<<<<<<<<Match: " << penalty.matchPenalty << " Mismatch: " << penalty.mismatchPenalty << " GapOpen: " << penalty.gapOpenPenalty << " GapExtend: " << penalty.gapExtendPenalty << ">>>>>>>>>>>>>" << endl;
            }
            Alignment aln;
            aln.readID = readID;
            aln.queryID = queryID;
            aln.read = read;
            aln.query = query;
            if (!isForwardStran) aln.flag = 16;
            align(aln, penalty.matchPenalty, penalty.mismatchPenalty, penalty.gapOpenPenalty, penalty.gapExtendPenalty);
            if(find(repeatedCigars.begin(), repeatedCigars.end(), aln.cigar) != repeatedCigars.end())
                continue;
            else
                repeatedCigars.push_back(aln.cigar);
            //check the alignment based on the criteria
            string cigarStr = convertCigarToStr(aln.cigar);
            cigarStr = trimClips(cigarStr);
            if(cigarStr.length() >= regionSize && checkIdentityPercentange(cigarStr)) {
                if(DEBUG_MODE) cout << "cigar len was smaller than the R or K at the beginning\n";
                if (aln.partialMatchSize > bestAlignment.partialMatchSize || (aln.partialMatchSize == bestAlignment.partialMatchSize && aln.editDistance < bestAlignment.editDistance))
                    bestAlignment = aln;
            }
        }
        return bestAlignment;
    }

    void findPartiaMatches(tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, uint32_t queryCount, bool isForwardStrand, string parmikAlignments, vector<Penalty> penalties)
    {
        ofstream pAln(parmikAlignments, ios::app);
        set<uint32_t> matchesPerQuery;
        cout << "Starting alignment for all queries [" << (isForwardStrand ? ("fwd"):("rev")) << "]..." << endl;
        for (size_t i = 0; i < queryCount; i++)
        {
            if(i % 10000 == 0)
                cout << i << " queries processed..." << endl;
            auto itq = queries.find(i);
            if (itq == queries.end())
                continue;
            string query = itq->second;
            if (query.find('n') != string::npos || query.find('N') != string::npos)
                continue;
            // read the candidate reads of cheap k-mer filter of front
            tsl::robin_map <uint32_t, Alignment> alignments;
            auto start = chrono::high_resolution_clock::now();
            for (size_t j = 0; j < reads.size(); j++)
            {
                auto itr = reads.find(j);
                if (itr == reads.end())
                    continue;
                string read = itr->second;
                Alignment aln = alignDifferentPenaltyScores(query, read, i, j, isForwardStrand, penalties);
                if (aln.partialMatchSize > 0){
                    alignments.insert(make_pair(it->first, aln));
                }
            }
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
            cout << "query [" << i << "]: # of matches found " << alignments.size() << " and it took: " << static_cast<double>(duration.count()) / 1'000'000.0 << "ms" << endl;
            //dump the alignments
            uint32_t matchesAccepted = 0;
            for (auto it = alignments.begin(); it!= alignments.end(); it++)
            {
                dumpSam(pAln, it->second);
                matchesAccepted++;
            }
            matchesPerQuery.insert(alignments.size());
            // cout << "queryID: " << i << ", " << (isForwardStrand ? ("fwd"):("rev")) <<", total matches: " << alignments.size() << i << ", matches accepted: " << matchesAccepted << ", matches passed with starting pos over ED: " << matchesPassedWithStartingPosOverED << endl;
        }
        Utilities<uint32_t> util; 
        tuple<uint32_t, uint32_t, uint32_t> matchesPerQueryTuple = util.calculateStatistics(matchesPerQuery);
        printf("PARMIK's matches per query up to q (%d) => [average: %d, median: %d, sum: %d]\n", queryCount, get<0>(matchesPerQueryTuple), get<1>(matchesPerQueryTuple), get<2>(matchesPerQueryTuple));
    }


};

#endif