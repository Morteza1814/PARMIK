CC = g++
CFLAGS = -std=c++20 -O3 -Wall -fopenmp 
# CFLAGS = -Wall -std=c++11 -Iinclude
all: parmik

parmik: main.o ssw_cpp.o ssw.o
	$(CC) $(CFLAGS) -o $@  $^ -lm -lz

main.o:  Includes/IndexContainer.h Includes/InvertedIndexBuilder.h Includes/KmersFrequencyCounter.h Includes/CheapKmerPartialMatcher.h Includes/args.hxx Includes/Configs.h Includes/BaseLinePartialMatcher.h Includes/Utils.h Includes/Container.h Includes/nlohmann/json.hpp Includes/IndexFile.h Includes/SeedMatchExtender.h Includes/LevDistanceCalculator.h Includes/CompareWithBWA.h Includes/BlastReader.h Includes/CompareWithBlast.h Includes/CheckKmersFrequency.h Includes/sw/ssw_cpp.h Includes/PostFilter.h Includes/Alignment.h Includes/SSW_BaseLine.h
	$(CC) $(CFLAGS) -o $@ main.cc -c -lm -lz

ssw_cpp.o: Includes/sw/ssw_cpp.cpp Includes/sw/ssw_cpp.h Includes/sw/ssw.h
	$(CXX) -c -o $@ $< $(CFLAGS) 

ssw.o: Includes/sw/ssw.c Includes/sw/ssw.h
	$(CC) -c -o $@ $< $(CFLAGS)

clean: 
	$(RM) *.o parmik

.PHONY: clean
