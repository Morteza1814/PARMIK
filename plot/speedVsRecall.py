import matplotlib.pyplot as plt
import sys
from matplotlib.lines import Line2D
import numpy as np

input_file = sys.argv[1]
output_file = sys.argv[2]

def read_values_from_file(file_path):
    # Initialize the variables
    percentage_identity = []
    parmik_matches_7 = []
    parmik_matches_4 = []
    blast_matches = []
    parmik_recall_7 = []
    parmik_recall_4 = []
    blast_recall = []

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            if 'percentage_identity' in line:
                percentage_identity = list(map(int, line.split('=')[1].strip().replace('[', '').replace(']', '').split(',')))
            elif 'parmik_matches_7' in line:
                parmik_matches_7 = list(map(float, line.split('=')[1].strip().replace('[', '').replace(']', '').split(',')))
            elif 'parmik_matches_4' in line:
                parmik_matches_4 = list(map(float, line.split('=')[1].strip().replace('[', '').replace(']', '').split(',')))
            elif 'blast_matches' in line:
                blast_matches = list(map(float, line.split('=')[1].strip().replace('[', '').replace(']', '').split(',')))
            elif 'parmik_recall_7' in line:
                parmik_recall_7 = list(map(float, line.split('=')[1].strip().replace('[', '').replace(']', '').split(',')))
            elif 'parmik_recall_4' in line:
                parmik_recall_4 = list(map(float, line.split('=')[1].strip().replace('[', '').replace(']', '').split(',')))
            elif 'blast_recall' in line:
                blast_recall = list(map(float, line.split('=')[1].strip().replace('[', '').replace(']', '').split(',')))
            elif 'title' in line:
                title = line.split('=')[1].strip()

    return (title, percentage_identity, parmik_matches_7, parmik_matches_4, blast_matches, 
            parmik_recall_7, parmik_recall_4, blast_recall)

# Example usage:
file_path = 'data.txt'  # Replace with the actual file path
title, percentage_identity, parmik_matches_7, parmik_matches_4, blast_matches, parmik_recall_7, parmik_recall_4, blast_recall = read_values_from_file(input_file)

# Create figure and axis objects without grids
fig, ax1 = plt.subplots()

# Increase the font size for clarity
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_xlabel('Percentage Identity (%)', fontsize=14)
ax1.set_ylabel('Speed (Matches/Second)', fontsize=14, color='tab:blue')
ax1.plot(percentage_identity, parmik_matches_7, label="PARMIK", marker='o', color='tab:blue')
ax1.plot(percentage_identity, parmik_matches_4, label="PARMIK", marker='o', color='tab:blue')
ax1.plot(percentage_identity, blast_matches, label="BLAST", marker='o', linestyle='--', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=12)

# Create a second y-axis for Recall Rate (%)
ax2 = ax1.twinx()
ax2.set_ylabel('Recall Rate (%)', fontsize=14, color='tab:red')
ax2.plot(percentage_identity, parmik_recall_7, label="PARMIK", marker='x', color='tab:red')
ax2.plot(percentage_identity, parmik_recall_4, label="PARMIK", marker='x', color='tab:red')
ax2.plot(percentage_identity, blast_recall, label="BLAST", marker='x', linestyle='--', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red', labelsize=12)

# Set the x-ticks to only show 85, 90, and 95
ax1.set_xticks(percentage_identity)
ax1.set_xticklabels([85, 90, 95], fontsize=12)

# Add the text at a position between 85 and 90 with the same slope
ax1.text(85 + 0.2, parmik_matches_7[0], 'PARMIK[M=7]', color='tab:blue', fontsize=10, ha='left',  va='bottom', style='italic')
ax1.text(85 + 0.2, parmik_matches_4[0], 'PARMIK[M=4]', color='tab:blue', fontsize=10, ha='left',  va='bottom', style='italic')
ax1.text(85 + 0.2, blast_matches[0]-120, 'BLAST', color='tab:blue', fontsize=10, ha='left',  va='bottom', style='italic')

ax2.text(95, parmik_recall_7[-1] - 10, 'PARMIK[M=7]', color='tab:red', fontsize=10, ha='right', va='top', style='italic')
ax2.text(95, parmik_recall_4[-1], 'PARMIK[M=4]', color='tab:red', fontsize=10, ha='right', va='bottom', style='italic')
ax2.text(95, blast_recall[-1] + 5, 'BLAST', color='tab:red', fontsize=10, ha='right', va='bottom', style='italic')


plt.title(title, fontsize=16)

# Adding legends with adjusted font sizes
fig.tight_layout()

# Remove the grid
ax1.grid(False)
ax2.grid(False)

# Show the plot
plt.savefig(output_file)
