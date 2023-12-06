import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import sys
import numpy as np

# Read the file and parse data
query_numbers = []
pm_reads = []
bwa_reads = []
blast_reads = []
gt_reads = []
bwa_fileName = sys.argv[1]
blast_fileName = sys.argv[2]
gt_fileName = sys.argv[3]
outFileName = sys.argv[4]

def calculate_average(arr):
    total_sum = sum(arr)
    length = len(arr)
    if length == 0:
        return 0  # Handle the case where the array is empty to avoid division by zero
    else:
        return total_sum / length
    

with open(bwa_fileName, 'r') as file:
    for line in file:
        data = line.strip().split()
        query_numbers.append(int(data[0]))
        pm_reads.append(int(data[1]))
        bwa_reads.append(int(data[2]))


with open(blast_fileName, 'r') as file:
    for line in file:
        data = line.strip().split()
        blast_reads.append(int(data[2]))

with open(gt_fileName, 'r') as file:
    for line in file:
        data = line.strip().split()
        gt_reads.append(int(data[2]))

# Create the figure
fig, ax = plt.subplots(figsize=(10, 6))

bar_width = 0.2  # Adjust the width of the bars
index = np.arange(len(query_numbers))

# Plot PM reads as bars
ax.bar(index, pm_reads, width=bar_width, color='green', label='PARMIK', alpha=0.7)

# Plot BLAST reads as bars
ax.bar(index + bar_width, blast_reads, width=bar_width, color='red', label='BLAST', alpha=0.7)
ax.set_xlabel('Query Number', fontsize=20)
ax.set_ylabel('Matches per Query', color='black', fontsize=25)
ax.tick_params(axis='y', labelcolor='black', labelsize=15)
ax.tick_params(axis='x', labelsize=15)

ax2 = ax.twinx()
ax2.plot(query_numbers, bwa_reads, color='dodgerblue', marker=',', linewidth=0.5, label='BWA')
ax2.set_ylabel('BWA Matches per Query (0/1)', color='dodgerblue', fontsize=20)
ax2.set_yticks([0, 8])
ax2.set_yticklabels(['0', '8'], fontsize = 15)

# Adjust the y-axis ticks for BWA reads
ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
ax2.tick_params(axis='y', labelcolor='dodgerblue', labelsize=15)

ax.legend(loc='upper center', fontsize=15)

# Add legends
# ax1.legend(loc='upper left')
# ax2.legend(loc='upper right')

# plt.title('Number of Reads per Query for PM (blue) and BWA (red)')

plt.savefig(outFileName)

# plt.show()
