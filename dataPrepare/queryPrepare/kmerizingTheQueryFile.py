import sys

def stringContainSetOfCharacters (str, set):
        return 1 in [c in str for c in set]

def stringContainNCharacter(input_str):
        # Check if 'N' or 'n' is in the input string
        if 'N' in input_str or 'n' in input_str:
                return True
        return False
    
args = []
repeatedQueries=0
uniqueQueries = []
queryOffset=0
queryCount=0 

args = sys.argv
conting_len = int(sys.argv[1])
genomeFileAddress = sys.argv[2]
outputKmersFileAddress = sys.argv[3]
genomeName = sys.argv[4]

print("conting_len = ", conting_len)
print("genomeFileAddress = ", genomeFileAddress)
print("outputKmersFileAddress = ", outputKmersFileAddress)
print("genomeName = ", genomeName)

f = open(genomeFileAddress, "r")
fo = open(outputKmersFileAddress, "a")
line = f.readline()
if not line:
        print("End of file!")
        exit
inputFileLen = len(line)
print("Length of the line: ", len(line))
queriesWithN = 0
while(queryOffset <= inputFileLen-conting_len):
        query = line[queryOffset:queryOffset+conting_len]
        if len(query) - query.count('\n') != conting_len:
                print("invalid query size: ", len(query))
                exit()
        if stringContainNCharacter(query):
               queriesWithN+=1
        if query in uniqueQueries:
                repeatedQueries+=1
        else:
                uniqueQueries.append(query)
        if stringContainSetOfCharacters(query, "ACGTNnacgt"):                        
                sequence = ">" + genomeName + "." + str(queryOffset) + "\n" +  query + "\n"
                fo.write(sequence)
                queryCount+=1
        queryOffset+=1

print("uniqueQueries len = ", len(uniqueQueries))
print("Number of Queries = ", queryCount)
print("repeatedQueries = ", repeatedQueries)
print("queriesWithN = ", queriesWithN)
f.close()
fo.close()
