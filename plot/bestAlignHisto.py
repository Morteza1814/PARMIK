import matplotlib.pyplot as plt
import numpy as np
import sys
# Function to read data from a file
def read_data_from_file(filename):
    data_blast_outperform_same_read = {}
    data_blast_outperform_diff_read = {}
    data_parmik_outperfrom_same_read = {}
    data_parmik_outperfrom_diff_read = {}
    data_parmik_aln_len_blast_fn = {}

    current_section = None

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "BLAST > PARMIK (same read)":
                current_section = data_blast_outperform_same_read
            elif line == "BLAST > PARMIK (different read)":
                current_section = data_blast_outperform_diff_read
            elif line == "PARMIK > BLAST (same read)":
                current_section = data_parmik_outperfrom_same_read
            elif line == "PARMIK > BLAST (different read)":
                current_section = data_parmik_outperfrom_diff_read
            elif line == "<<<<<<<<<<<<<<<<<<<<<<<<PARMIK best aln size histo when BLAST FN>>>>>>>>>>>>>>>>>>>>>>>":
                current_section = data_parmik_aln_len_blast_fn
            elif current_section is not None and ":" in line:
                try:
                    if line:
                        key, value = map(int, line.split(":"))
                        current_section[key] = value
                except ValueError:
                    pass
    return data_blast_outperform_same_read, data_blast_outperform_diff_read, data_parmik_outperfrom_same_read, data_parmik_outperfrom_diff_read, data_parmik_aln_len_blast_fn

# Read the data from file
filename = sys.argv[1]
outputFile = sys.argv[2]
data_blast_outperform_same_read, data_blast_outperform_diff_read, data_parmik_outperfrom_same_read, data_parmik_outperfrom_diff_read, data_parmik_aln_len_blast_fn = read_data_from_file(filename)

# All keys
all_keys_extended = set(data_parmik_outperfrom_same_read.keys()).union(data_parmik_outperfrom_diff_read.keys()).union(data_blast_outperform_same_read.keys()).union(data_blast_outperform_diff_read.keys()).union(data_parmik_aln_len_blast_fn.keys())
all_keys_extended = sorted(all_keys_extended)

# Data preparation
same_read_counts_parmik = [data_parmik_outperfrom_same_read.get(key, 0) for key in all_keys_extended]
diff_read_counts_parmik = [data_parmik_outperfrom_diff_read.get(key, 0) for key in all_keys_extended]
same_read_counts_blast = [data_blast_outperform_same_read.get(key, 0) for key in all_keys_extended]
diff_read_counts_blast = [data_blast_outperform_diff_read.get(key, 0) for key in all_keys_extended]
diff_blast_fn_counts_parmik = [data_parmik_aln_len_blast_fn.get(key, 0) for key in all_keys_extended]

bar_width = 0.15
positions = np.arange(len(all_keys_extended))

# Create histogram
plt.figure(figsize=(16, 8))

# Plotting histograms side by side with log scale
plt.bar(positions, same_read_counts_parmik, width=bar_width, color='blue', label='PARMIK > BLAST (same read)')
plt.bar(positions + bar_width, diff_read_counts_parmik, width=bar_width, color='orange', label='PARMIK > BLAST (different read)')
plt.bar(positions + 2 * bar_width, same_read_counts_blast, width=bar_width, color='green', label='BLAST > PARMIK (same read)')
plt.bar(positions + 3 * bar_width, diff_read_counts_blast, width=bar_width, color='red', label='BLAST > PARMIK (different read)')
plt.bar(positions + 4 * bar_width, diff_blast_fn_counts_parmik, width=bar_width, color='purple', label='PARMIK best aln size (BLAST FN)')

# Adding titles and labels
plt.title('Histogram of Times PARMIK and BLAST Outperformed Each Other in Terms of Alignment Size (Log Scale)')
plt.xlabel('Alignment Size')
plt.ylabel('Count (Log Scale)')
plt.yscale('log')
plt.xticks(positions + 2 * bar_width, all_keys_extended, rotation=90)
plt.legend()

# Show plot
plt.savefig(outputFile)
