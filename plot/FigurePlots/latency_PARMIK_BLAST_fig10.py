import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys


def readTime(timeFilePath):
    with open(timeFilePath, 'r') as file:
        time = int(file.readline().strip())
    return time

datasets = []
metagenomicDatasets = ['SRR12432009', 'SRR14381418', 'SRR14381421', 'SRR14381422']
queryDatasets = ['NC_045512', 'SarsGenome']

tools = ['PARMIK [M=4]', 'PARMIK [M=7]', 'BLAST']

for metagenomicDataset in metagenomicDatasets:
    for queryDataset in queryDatasets:
        if queryDataset == 'NC_045512':
            datasetName = 'Wuhan-Hu-1' + '\n' + metagenomicDataset
        else:
            datasetName = queryDataset + '\n' + metagenomicDataset
        datasets.append(datasetName)

# Set up positions for equal spacing of tools across datasets
x_positions = np.arange(len(datasets) * len(tools)) * 1.5  # 1.5 unit spacing for clarity
x_ticks_tools = [x_positions[i * len(tools):(i + 1) * len(tools)] for i in range(len(datasets))]

# Create labels with tool name and dataset name
x_labels_tools = tools * len(datasets)  # Repeating tool names across all datasets

fig, ax = plt.subplots(figsize=(18, 8))  # Adjust width for a two-column layout

IKT = 0
M = [4, 7]
PI = [85, 90, 95]
T = 1
SC = 1
K = 11
R = 30
markers = ['o', 'x', 's', 'D', '^']
colors = ['blue', 'green', 'red']  # Consistent colors for each PI level
experimentPath = sys.argv[1]
# Plot the data for each dataset and PI level with new x positions
i=0
for metagenomicDataset in metagenomicDatasets:
    for queryDataset in queryDatasets:
        for j, pi in enumerate(PI):
            times = []
            for m in M:
                dir = experimentPath + '/' + metagenomicDataset + '/' + queryDataset + '/IKT' + str(IKT) + '_K' + str(K) + '_PI' + str(pi) + '_M' + str(m) + '_T' + str(T) + '_SC' + str(SC) + '_P_2484_1111_2444_2288_2848/'
                parmikTimePath = dir + 'time_parmik'
                times.append(readTime(parmikTimePath))
                if m == 4:
                    blastTimePath = dir + 'time_blast'
                    times.append(readTime(blastTimePath))
            
            # Define x positions for the tools in the current dataset
            x = x_ticks_tools[i]
            
            # Plot timing for each tool at the current PI level
            ax.plot(x, times, linestyle='', marker=markers[j], color=colors[j], 
                    label=f'PI={pi}' if i == 0 else "", linewidth=1.5, markersize=10)
        i+=1

ax.set_ylabel('Execution Time (minutes)', fontsize=18)

# Remove major ticks on the x-axis
ax.set_xticks([])  # Removes the major ticks

# Adding dataset labels in boxes at the top of the figure
dataset_tick_positions = [i * len(tools) * 1.5 + (len(tools) - 1) / 2 * 1.5 for i in range(len(datasets))]
for i, dataset in enumerate(datasets):
    # Adding text boxes above the plot
    ax.text(dataset_tick_positions[i], ax.get_ylim()[1] + 10, dataset, 
            ha='center', va='bottom', fontsize=16, bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))

plt.yticks(fontsize=14)

# Set tool names on x-axis using minor ticks
ax.set_xticks(x_positions, minor=True)  # Only tool positions as minor ticks
ax.set_xticklabels(x_labels_tools, minor=True, rotation=90, ha='center', fontsize=14)

# Adding vertical grid lines for each tool
ax.grid(True, axis='x', which='both', linestyle=':', color='grey')

# Adding vertical lines to separate datasets
for i in range(1, len(datasets)):
    plt.axvline(x=x_ticks_tools[i][0] - 0.75, color='black', linestyle='--', linewidth=1)

# Show PI level in the legend
ax.legend(loc='upper right', fontsize=14)

plt.tight_layout()
plt.show()
outputPath = sys.argv[2]
plt.savefig(outputPath + '/allTime.pdf')
