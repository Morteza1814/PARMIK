CC = g++
CFLAGS = -std=c++20 -O3 -Wall 
# CFLAGS = -Wall -std=c++11 -Iinclude
all: parmik

parmik: main.o
	$(CC) $(CFLAGS) -o $@  $^

main.o:  Includes/IndexContainer.h Includes/InvertedIndexBuilder.h Includes/KmersFrequencyCounter.h Includes/CheapKmerPartialMatcher.h Includes/args.hxx Includes/Configs.h Includes/BaseLinePartialMatcher.h Includes/Utils.h Includes/Container.h Includes/nlohmann/json.hpp Includes/IndexFile.h Includes/SeedMatchExtender.h Includes/LevDistanceCalculator.h Includes/CompareWithBWA.h Includes/BlastReader.h Includes/CompareWithBlast.h Includes/CheckKmersFrequency.h
	$(CC) $(CFLAGS) -o $@ main.cc -c

clean: 
	$(RM) *.o parmik

.PHONY: clean
