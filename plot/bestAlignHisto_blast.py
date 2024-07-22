import matplotlib.pyplot as plt
import numpy as np
import sys

# Function to read data from a file
def read_data_from_file(filename):
    data_blast_same_read = {}
    data_blast_diff_read = {}
    data_same_read = {}
    data_diff_read = {}
    data_parmik_diff_blast_fn = {}

    current_section = None

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line == "BLAST > PARMIK (same read)":
                current_section = data_blast_same_read
            elif line == "BLAST > PARMIK (different read)":
                current_section = data_blast_diff_read
            elif line == "PARMIK > BLAST (same read)":
                current_section = data_same_read
            elif line == "PARMIK > BLAST (different read)":
                current_section = data_diff_read
            elif "when BLAST FN" in line:
                current_section = data_parmik_diff_blast_fn
            elif current_section is not None and ":" in line:
                try:
                    key, value = map(int, line.split(":"))
                    current_section[key] = value
                except ValueError:
                    pass
    # print(data_blast_same_read)
    # print(data_blast_diff_read)
    # print(data_same_read)
    # print(data_diff_read)
    # print(data_parmik_diff_blast_fn)
    return data_blast_same_read, data_blast_diff_read, data_same_read, data_diff_read, data_parmik_diff_blast_fn

# Function to cluster data into buckets of size 10
def cluster_data(data, bucket_size=10):
    clustered_data = {}
    for key, value in data.items():
        bucket = (key // bucket_size) * bucket_size
        if bucket in clustered_data:
            clustered_data[bucket] += value
        else:
            clustered_data[bucket] = value
    return clustered_data

# Function to create histogram
def create_histogram(data_blast_same_read, data_blast_diff_read, data_same_read, data_diff_read, data_parmik_diff_blast_fn, output_filename):
    bucket_size = 10
    data_blast_same_read = cluster_data(data_blast_same_read, bucket_size)
    data_blast_diff_read = cluster_data(data_blast_diff_read, bucket_size)
    data_same_read = cluster_data(data_same_read, bucket_size)
    data_diff_read = cluster_data(data_diff_read, bucket_size)
    data_parmik_diff_blast_fn = cluster_data(data_parmik_diff_blast_fn, bucket_size)

    all_keys_extended = sorted(set(data_blast_same_read.keys()).union(
        data_blast_diff_read.keys(), data_same_read.keys(), 
        data_diff_read.keys(), data_parmik_diff_blast_fn.keys()))

    same_read_counts_parmik = [data_same_read.get(key, 0) for key in all_keys_extended]
    diff_read_counts_parmik = [data_diff_read.get(key, 0) for key in all_keys_extended]
    same_read_counts_blast = [data_blast_same_read.get(key, 0) for key in all_keys_extended]
    diff_read_counts_blast = [data_blast_diff_read.get(key, 0) for key in all_keys_extended]
    diff_blast_fn_counts_parmik = [data_parmik_diff_blast_fn.get(key, 0) for key in all_keys_extended]

    bar_width = 0.15
    positions = np.arange(len(all_keys_extended))

    plt.figure(figsize=(16, 8))

    plt.bar(positions, same_read_counts_parmik, width=bar_width, color='dodgerblue', label='PARMIK > BLAST (same read)')
    plt.bar(positions + bar_width, diff_read_counts_parmik, width=bar_width, color='darkorange', label='PARMIK > BLAST (different read)')
    plt.bar(positions + 2 * bar_width, same_read_counts_blast, width=bar_width, color='black', label='BLAST > PARMIK (same read)')
    plt.bar(positions + 3 * bar_width, diff_read_counts_blast, width=bar_width, color='grey', label='BLAST > PARMIK (different read)')
    plt.bar(positions + 4 * bar_width, diff_blast_fn_counts_parmik, width=bar_width, color='green', label='PARMIK best aln size (BLAST FN)')

    plt.yscale('log')
    # plt.title('Histogram of Times PARMIK and BLAST Outperformed Each Other in Terms of Alignment Size', fontsize=16)
    plt.xlabel('Number of Base Pairs (Ranges)', fontsize=28)
    plt.ylabel('Count (Log Scale)', fontsize=28)

    xticks_labels = [f"[{key}-{key + bucket_size - 1}]" for key in all_keys_extended]
    plt.xticks(positions + 2 * bar_width, xticks_labels, fontsize=22)
    plt.yticks(fontsize=22)
    plt.legend(fontsize=24)

    plt.tight_layout()
    plt.savefig(output_filename)
    plt.show()

# File name
filename = sys.argv[1]
output_filename = sys.argv[2]

# Read the data from file
data_blast_same_read, data_blast_diff_read, data_same_read, data_diff_read, data_parmik_diff_blast_fn = read_data_from_file(filename)

# Create histogram
create_histogram(data_blast_same_read, data_blast_diff_read, data_same_read, data_diff_read, data_parmik_diff_blast_fn, output_filename)
