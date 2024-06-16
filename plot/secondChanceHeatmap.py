import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def read_data_from_file(file_name):
    data = []
    with open(file_name, 'r') as file:
        for line in file:
            parts = line.strip().split()
            x_range = tuple(map(int, parts[0][1:-1].split(',')))
            y_range = tuple(map(int, parts[1][1:-1].split(',')))
            count = int(parts[2])
            data.append((x_range, y_range, count))
    return data

def normalize_counts(counts, min_val, max_val):
    return [(count - min_val) / (max_val - min_val) if max_val > min_val else 0 for count in counts]

def plot_heatmap(data):
    # Define the normalization ranges
    ranges = [
        (0, 10),
        (500, 2000),
        (3000, 20000),
        (800000, 2000000),
        (4000000, max(item[2] for item in data))
    ]
    
    # Define the color maps for each range
    color_maps = [
        'Blues',
        'Greens',
        'Oranges',
        'Reds',
        'Purples'
    ]
    
    # Extract the x and y ranges for plotting
    x_ranges = [(x_range[0], x_range[1]) for (x_range, _, _) in data]
    y_ranges = [(y_range[0], y_range[1]) for (_, y_range, _) in data]
    counts = [item[2] for item in data]

    # Prepare the plot
    plt.figure(figsize=(14, 10))

    # Create a dictionary to store the heatmap values
    heatmap_dict = {}
    for (x_range, y_range), count in zip(zip(x_ranges, y_ranges), counts):
        heatmap_dict[(x_range, y_range)] = count

    # Create lists of unique x and y midpoints
    x_unique = sorted(set(x_ranges))
    y_unique = sorted(set(y_ranges))

    # Create a grid for the heatmap
    # heatmap_data = np.zeros((len(y_unique), len(x_unique)))
    heatmap_rawdata = np.zeros((len(y_unique), len(x_unique)))

    # Apply normalization for each range and map to color scales
    for (min_val, max_val), cmap in zip(ranges, color_maps):
        # Filter counts within the current range
        filtered_counts = {k: v for k, v in heatmap_dict.items() if min_val <= v <= max_val}
        normalized_counts = normalize_counts(list(filtered_counts.values()), min_val, max_val)
        
        # Create a grid for the heatmap
        heatmap_data = np.zeros((len(y_unique), len(x_unique)))
        
        # Map normalized counts to the heatmap grid
        for (x_range, y_range), norm_count, count in zip(filtered_counts.keys(), normalized_counts, filtered_counts.values()):
            x_index = x_unique.index(x_range)
            y_index = y_unique.index(y_range)
            heatmap_data[y_index, x_index] = norm_count
            heatmap_rawdata[y_index, x_index] = count
        
        # Plot the heatmap
        ax = sns.heatmap(heatmap_data, 
                     xticklabels=[f'{x[0]}-{x[1]}' for x in x_unique], 
                    yticklabels=[f'{y[0]}-{y[1]}' for y in y_unique], 
                    cmap=cmap, alpha = 0.2, cbar=False)
                    # cbar_kws={'label': f'Counts {min_val}-{max_val}'})
    # Annotate the heatmap with count values
    for i in range(len(y_unique)):
        for j in range(len(x_unique)):
            ax.text(j + 0.5, i + 0.5, int(heatmap_rawdata[i, j]),
                    ha='center', va='center', color='black', fontsize=8)

    # Label the axes
    plt.xlabel('Alignment Size Range (bp')
    plt.ylabel('Alignment Size Improvement Range (bp)')
    
    # Add title
    # plt.title('Heatmap of Normalized Counts with X and Y Ranges')
    # Reverse the y-axis to have 0-9 at the bottom
    ax.invert_yaxis()
    # Show the plot
    # plt.show()
    plt.savefig('secondChanceheatmap.png')

if __name__ == "__main__":
    # Check if the filename is passed as an argument
    if len(sys.argv) < 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    file_name = sys.argv[1]
    data = read_data_from_file(file_name)
    plot_heatmap(data)
