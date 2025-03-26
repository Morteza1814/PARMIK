import sys
from tqdm import tqdm
import matplotlib.pyplot as plt
from collections import defaultdict


def read_fasta_to_set(file_path):
    """Reads a FASTA file and returns a set of sequences."""
    sequences = set()
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequences.add(line.strip())
    return sequences


def read_fasta_to_dict(file_path):
    """Reads a FASTA file and returns a dictionary {contig_name: sequence}"""
    sequences = {}
    with open(file_path, 'r') as f:
        contig_name = ""
        for line in f:
            if line.startswith('>'):
                contig_name = line.strip().lstrip('>')
            else:
                sequences[contig_name] = line.strip()
    return sequences


def extract_kmers(sequence, k):
    """Extracts kmers of size k from a sequence."""
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]


def analyze_contigs(contigs_file, kmers_file, classification_file, output_stats_file, k=11):
    expensive_kmers = read_fasta_to_set(kmers_file)
    contigs = read_fasta_to_dict(contigs_file)

    contig_classification = {}
    with open(classification_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            contig_classification[parts[1]] = parts[0]  # C or U

    stats = defaultdict(lambda: {'C': 0, 'U': 0})

    for contig_name, sequence in tqdm(contigs.items(), desc="Processing contigs"):
        contig_kmers = extract_kmers(sequence, k)
        expensive_count = sum(1 for kmer in contig_kmers if kmer in expensive_kmers)
        classification = contig_classification.get(contig_name, 'U')
        stats[expensive_count][classification] += 1

    with open(output_stats_file, 'w') as out:
        out.write("Expensive_Kmers\tClassified\tUnclassified\n")
        for expensive_count in sorted(stats.keys()):
            out.write(f"{expensive_count}\t{stats[expensive_count]['C']}\t{stats[expensive_count]['U']}\n")


def plot_stats(stats_file):
    labels, classified_counts, unclassified_counts = [], [], []

    with open(stats_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            labels.append(int(parts[0]))
            classified_counts.append(int(parts[1]))
            unclassified_counts.append(int(parts[2]))

    plt.bar(labels, classified_counts, label='Classified', color='skyblue')
    plt.bar(labels, unclassified_counts, bottom=classified_counts, label='Unclassified', color='lightcoral')

    plt.xlabel('Number of Expensive Kmers')
    plt.ylabel('Number of Sequences (log scale)')
    plt.yscale('log')
    plt.title('Sequence Classification by Expensive Kmers (Log Scale)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('classification_kmers_distribution.png')
    plt.show()


if __name__ == "__main__":
    # if len(sys.argv) not in [1, 5]:
    #     print("Usage:")
    #     print("Analyze and output stats: python script.py analyze <contigs_file> <expensive_kmers_file> <classification_file>")
    #     print("Plot stats: python script.py plot <output_stats_file>")
    #     sys.exit(1)

    if sys.argv[1] == "analyze" and len(sys.argv) == 5:
        contigs_file, kmers_file, classification_file = sys.argv[2], sys.argv[3], sys.argv[4]
        analyze_contigs(contigs_file, kmers_file, classification_file, 'output_stats.txt', k=11)
    elif sys.argv[1] == "plot" and len(sys.argv) == 3:
        stats_file = sys.argv[2]
        plot_stats(stats_file)
    else:
        print("Invalid arguments. Please follow usage instructions.")