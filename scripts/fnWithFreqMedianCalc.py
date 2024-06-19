import sys
import numpy as np

def read_data(filename):
    alignment_sizes = []
    frequencies = []
    with open(filename, 'r') as file:
        for line in file:
            size, frequency = map(int, line.split())
            alignment_sizes.append(size)
            frequencies.append(frequency)
    return alignment_sizes, frequencies

def weighted_median(data, weights):
    data, weights = np.array(data), np.array(weights)
    sorted_indices = np.argsort(data)
    data, weights = data[sorted_indices], weights[sorted_indices]
    cumulative_weight = np.cumsum(weights)
    midpoint = 0.5 * sum(weights)
    if any(weights > midpoint):
        return (data[weights == np.max(weights)])[0]
    else:
        median_idx = np.where(cumulative_weight >= midpoint)[0][0]
        return data[median_idx]

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    alignment_sizes, frequencies = read_data(filename)
    median_alignment_size = weighted_median(alignment_sizes, frequencies)
    print("Median alignment size:", median_alignment_size)

if __name__ == "__main__":
    main()
