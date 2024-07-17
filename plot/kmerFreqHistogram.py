import matplotlib.pyplot as plt
import sys
import numpy as np

def read_buckets_and_frequencies(file_path):
    buckets = []
    kmer_freqs_9 = []
    kmer_freqs_10 = []
    kmer_freqs_11 = []
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                parts = line.split('\t')
                boundary_part = parts[0].strip()
                bounds = boundary_part.strip('[]').split('-')
                lower_bound = int(bounds[0].replace(',', ''))
                upper_bound = int(bounds[1].replace(',', ''))
                buckets.append((lower_bound, upper_bound))
                
                kmer_freqs_9.append(int(parts[1].replace(',', '')))
                kmer_freqs_10.append(int(parts[2].replace(',', '')))
                kmer_freqs_11.append(int(parts[3].replace(',', '')))
                
    return buckets, kmer_freqs_9, kmer_freqs_10, kmer_freqs_11

def calculate_percentages(kmer_freqs, buckets):
    total_kmers = sum(kmer_freqs)
    percentages = []
    cumulative = 0
    
    for s in kmer_freqs:
        cumulative += s
        percentages.append((cumulative / total_kmers) * 100)
    
    return percentages

def format_value(value):
    if value >= 1_000_000:
        return f'{value // 1_000_000}M'
    elif value >= 1_000:
        return f'{value // 1_000}K'
    else:
        return str(value)

def format_bucket_label(lower, upper):
    return f'[{format_value(lower)}-{format_value(upper)}]'

def plot_histograms(buckets, kmer_freqs_list, labels, outputFile):
    n = len(buckets)
    ind = np.arange(n)  # the x locations for the groups
    width = 0.2  # the width of the bars

    fig, ax = plt.subplots(figsize=(14, 8))
    for i, (kmer_freqs, label) in enumerate(zip(kmer_freqs_list, labels)):
        percentages = calculate_percentages(kmer_freqs, buckets)
        rects = ax.bar(ind + i * width, percentages, width, label=label)
        annotated = False
        for rect in rects:
            height = rect.get_height()
            if height > 99 and not annotated:
                annotated = True
                color = rect.get_facecolor()
                ax.annotate(f'{height:.1f}%', 
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  
                            textcoords="offset points",
                            ha='center', va='bottom',
                            fontsize=22, color=color)
    
    # Create formatted bucket labels
    bucket_labels = [format_bucket_label(0, high) for low, high in buckets]
    
    ax.set_ylabel('Cumulative Percentage of k-mers', fontsize=26)
    ax.set_xticks(ind + width)
    ax.set_xticklabels(bucket_labels, rotation=45, fontsize=20)
    # Set the y-tick positions
    ax.set_ylim(0, 115)
    yticks = np.arange(0, 115, 25)
    ax.set_yticks(yticks)
    # Set the y-tick labels with appropriate font size
    ax.set_yticklabels(yticks, fontsize=20)
    ax.set_xlabel('Inexpensive K-mer Threshold Range', fontsize=26)
    ax.legend(fontsize=18)

    plt.tight_layout()
    plt.savefig(outputFile)
    plt.show()  # Ensure the plot is shown if running interactively


# Example usage
file_path = sys.argv[1]
outputFile = sys.argv[2]
buckets, kmer_freqs_9, kmer_freqs_10, kmer_freqs_11 = read_buckets_and_frequencies(file_path)
kmer_freqs_list = [kmer_freqs_9, kmer_freqs_10, kmer_freqs_11]
labels = ['k=9', 'k=10', 'k=11']  # Labels for different k-mer sizes

plot_histograms(buckets, kmer_freqs_list, labels, outputFile)
