import matplotlib.pyplot as plt
import sys
from matplotlib.lines import Line2D

output_file = input_file1 = sys.argv[1]
# Data from the table
percentage_identity = [85, 90, 95]
parmik_matches = [4789.330496, 4587.39375, 3986.462179]
blast_matches = [4074.732313, 3970.54058, 3272.246591]
parmik_recall = [69.27393661, 96.17265882, 99.63863281]
blast_recall = [62.81613251, 84.5726937, 83.94311319]

# Create figure and axis objects without grids
fig, ax1 = plt.subplots()

# Increase the font size for clarity
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_xlabel('Percentage Identity (%)', fontsize=14)
ax1.set_ylabel('Speed (Matches/Second)', fontsize=14, color='tab:blue')
ax1.plot(percentage_identity, parmik_matches, label="PARMIK", marker='o', color='tab:blue')
ax1.plot(percentage_identity, blast_matches, label="BLAST", marker='o', linestyle='--', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=12)

# Create a second y-axis for Recall Rate (%)
ax2 = ax1.twinx()
ax2.set_ylabel('Recall Rate (%)', fontsize=14, color='tab:red')
ax2.plot(percentage_identity, parmik_recall, label="PARMIK", marker='x', color='tab:red')
ax2.plot(percentage_identity, blast_recall, label="BLAST", marker='x', linestyle='--', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red', labelsize=12)

# Set the x-ticks to only show 85, 90, and 95
ax1.set_xticks(percentage_identity)
ax1.set_xticklabels([85, 90, 95], fontsize=12)


# Create custom lines
custom_lines = [Line2D([0], [0], color='tab:grey', linestyle='-', label='PARMIK'),
                Line2D([0], [0], color='tab:grey', linestyle='--', label='BLAST')]

# Adding the custom legend
plt.legend(custom_lines, ['PARMIK', 'BLAST'], loc='lower center', fontsize=12)
# Adding legends with adjusted font sizes
fig.tight_layout()
# ax1.legend(loc='upper left', fontsize=12)
# ax2.legend(loc='upper right', fontsize=12)

# Remove the grid
ax1.grid(False)
ax2.grid(False)

# Show the plot
plt.savefig(output_file)
