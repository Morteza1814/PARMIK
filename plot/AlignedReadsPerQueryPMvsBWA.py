import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sys

# Read the file and parse data
query_numbers = []
pm_reads = []
bwa_reads = []
fileName = sys.argv[1]
outFileName = sys.argv[2]

def calculate_average(arr):
    total_sum = sum(arr)
    length = len(arr)
    if length == 0:
        return 0  # Handle the case where the array is empty to avoid division by zero
    else:
        return total_sum / length
    

with open(fileName, 'r') as file:
    for line in file:
        data = line.strip().split()
        query_numbers.append(int(data[0]))
        pm_reads.append(int(data[1]))
        bwa_reads.append(int(data[2]))

# Create the figure
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot PM reads as bars
ax1.bar(query_numbers, pm_reads, color='dodgerblue', label='PM Reads', alpha=0.7)
# ax1.plot(query_numbers, pm_reads, color='blue', marker=',', linewidth=0.5, label='PM Reads')
ax1.set_xlabel('Query Number', fontsize=20)
ax1.set_ylabel('PARMIK Reads', color='dodgerblue', fontsize=20)
ax1.tick_params(axis='y', labelcolor='dodgerblue', labelsize=15)
ax1.tick_params(axis='x', labelsize=15)

ax2 = ax1.twinx()
ax2.plot(query_numbers, bwa_reads, color='red', marker=',', linewidth=0.5, label='BWA Reads')
ax2.set_ylabel('BWA Reads (0/1)', color='red', fontsize=20)
ax2.set_yticks([0, 1])
ax2.set_yticklabels(['0', '1'])

# Adjust the y-axis ticks for BWA reads
ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
ax2.tick_params(axis='y', labelcolor='red', labelsize=15)

# Add legends
# ax1.legend(loc='upper left')
# ax2.legend(loc='upper right')

# plt.title('Number of Reads per Query for PM (blue) and BWA (red)')

plt.savefig(outFileName)

total_sum = sum(pm_reads)
print("The sum of number of reads PM has found:", total_sum)
print("The average of number of reads per query PM has found:", calculate_average(pm_reads))
# plt.show()
