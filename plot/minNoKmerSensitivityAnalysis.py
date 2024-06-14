import matplotlib.pyplot as plt
import sys
import numpy as np

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script_name.py input_file output_file")
    sys.exit(1)

# Get input and output file names from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read data from file
M = []
FN_normalized = []
aln_time = []

with open(input_file, 'r') as file:
    # Skip the header line
    # next(file)
    for line in file:
        data = line.split()
        M.append(int(data[0]))
        FN_normalized.append(int(data[1]))
        aln_time.append(int(data[2]))

# Create a figure and a set of subplots
fig, ax1 = plt.subplots()

# Plot FN_normalized
ax1.set_xlabel('Minimum No. of Shared K-mers', fontsize=16)
ax1.set_ylabel('No. of False Negative', color='tab:blue',  fontsize=16)
ax1.plot(M, FN_normalized, color='tab:blue', marker='o', linestyle='-')
ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=12)

# Customize x-axis ticks
ax1.set_xticks(M)
ax1.set_xticklabels(M, fontsize=12)

# Create a second y-axis
ax2 = ax1.twinx()
ax2.set_ylabel('Execution Time (Minutes)', color='tab:red', fontsize=16)
ax2.plot(M, aln_time, color='tab:red', marker='x', linestyle='--')
ax2.tick_params(axis='y', labelcolor='tab:red', labelsize=12)

# Add a title
# plt.title('FN Normalized and Alignment Time vs M')
plt.tight_layout()

# Save the plot to the output file
plt.savefig(output_file)
