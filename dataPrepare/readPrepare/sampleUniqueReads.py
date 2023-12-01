import sys

def stringContainNCharacter(input_str):
        # Check if 'N' or 'n' is in the input string
        if 'N' in input_str or 'n' in input_str:
                return True
        return False

CONTIG_LEN = 150
lines = []
args = []
readCount = 0
queryCount=0
repeatedReads=0
linesRead=0
uniqueReads = {}

args = sys.argv
metagenomicFileName = sys.argv[1]
expectedQueryCount = int(sys.argv[2])
skipN = int(sys.argv[3])
headerLineSubstr = sys.argv[4]
readDatabaseAddress = sys.argv[5]
sampleFileAddress = sys.argv[6]

print("MetagenomicFileName = ", metagenomicFileName)
print("MaxexpectedQueryCount = ", expectedQueryCount)
print("skip Ns? = ", skipN)
print("headerLineSubstr = ", headerLineSubstr)
print("readDatabaseAddress = ", readDatabaseAddress)
print("outputUniqueReadsFileAddress = ", sampleFileAddress)

f = open(readDatabaseAddress, "r")
fo = open(sampleFileAddress, "a")

linesContainingN = 0
while(True):
        line = f.readline()
        linesRead+=1
        if not line:
                print("End of file!")
                break
        if headerLineSubstr not in line and len(line) > 2:
                if (len(line) - line.count('\n')) != CONTIG_LEN:
                        print("invalid read size!!!")
                        break
                readCount += 1
                if skipN == 1 and stringContainNCharacter(line):
                        linesContainingN+=1
                        continue
                sequence = ">" + metagenomicFileName + "." + str(queryCount) + "\n" +  line
                if line not in uniqueReads:
                        uniqueReads[line] = queryCount    
                        fo.write(sequence)
                        queryCount+=1
                else:
                        repeatedReads+=1
        if queryCount >= expectedQueryCount:
                break

print("linesRead is : ", linesRead)
print("linesContainingN is : ", linesContainingN)
print("readCount is : ", readCount)
print("queryCount is : ", queryCount)
print("repeatedReads is : ", repeatedReads)

f.close()
fo.close()
