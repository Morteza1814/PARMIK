CC = g++
CFLAGS = -std=c++20 -O3 -Wall -fopenmp 
OBJDIR = build
SRCDIR = src

all: parmik

parmik: $(OBJDIR)/main.o $(OBJDIR)/ssw_cpp.o $(OBJDIR)/ssw.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -lz

$(OBJDIR)/main.o: $(SRCDIR)/Includes/IndexContainer.h $(SRCDIR)/Includes/InvertedIndexBuilder.h $(SRCDIR)/Includes/KmersFrequencyCounter.h $(SRCDIR)/Includes/CheapKmerPartialMatcher.h $(SRCDIR)/Includes/args.hxx $(SRCDIR)/Includes/Configs.h $(SRCDIR)/Includes/BaseLinePartialMatcher.h $(SRCDIR)/Includes/Utils.h $(SRCDIR)/Includes/Container.h $(SRCDIR)/Includes/nlohmann/json.hpp $(SRCDIR)/Includes/IndexFile.h $(SRCDIR)/Includes/SeedMatchExtender.h $(SRCDIR)/Includes/LevDistanceCalculator.h $(SRCDIR)/Includes/CompareWithBWA.h $(SRCDIR)/Includes/BlastReader.h $(SRCDIR)/Includes/CompareWithBlast.h $(SRCDIR)/Includes/CheckKmersFrequency.h $(SRCDIR)/Includes/sw/ssw_cpp.h $(SRCDIR)/Includes/PostFilter.h $(SRCDIR)/Includes/Alignment.h $(SRCDIR)/Includes/SSW_BaseLine.h $(SRCDIR)/Includes/CompareWithBaseLine.h $(SRCDIR)/Includes/ExpensiveKmersFNEvaluator.h $(SRCDIR)/Includes/EvaluateSecondChance.h
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -o $@ $(SRCDIR)/main.cc -c -lm -lz

$(OBJDIR)/ssw_cpp.o: $(SRCDIR)/Includes/sw/ssw_cpp.cpp $(SRCDIR)/Includes/sw/ssw_cpp.h $(SRCDIR)/Includes/sw/ssw.h
	@mkdir -p $(OBJDIR)
	$(CXX) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/ssw.o: $(SRCDIR)/Includes/sw/ssw.c $(SRCDIR)/Includes/sw/ssw.h
	@mkdir -p $(OBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	$(RM) $(OBJDIR)/*.o parmik
	@rmdir $(OBJDIR)

.PHONY: clean
