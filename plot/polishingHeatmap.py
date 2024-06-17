import sys
import matplotlib.pyplot as plt
import numpy as np

# Define the ranges for the first and second columns
first_column_ranges = [
    (10, 29),
    (30, 59),
    (60, 99),
    (100, 129),
    (130, 159),
]

second_column_ranges = [
    (0, 0),
    (0, 9),
    (10, 29),
    (30, 59),
    (60, 99),
]

# Helper function to format numbers with 'K' and 'M'
def format_count(count):
    if count >= 1_000_000:
        return f"{count // 1_000_000}M"
    elif count >= 1_000:
        return f"{count // 1_000}K"
    else:
        return str(count)
    
# Helper function to find the range index for a given value
def find_range_index(value, ranges):
    for index, (start, end) in enumerate(ranges):
        if start <= value <= end:
            return index
    return None

# Check if the file name is provided
if len(sys.argv) < 3:
    print("Please provide the input file and output file name as a command-line argument.")
    sys.exit(1)

# Read the file name from command-line arguments
input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

# Read the data from the file
with open(input_file_name, 'r') as file:
    data_lines = file.readlines()

# Parse the data into a list of tuples
data_pairs = [tuple(map(int, line.split())) for line in data_lines]

# Initialize the 2D list (matrix) to hold counts
counts_matrix = [[0 for _ in range(len(first_column_ranges))] for _ in range(len(second_column_ranges))]

# Categorize and count the data
for first_value, second_value in data_pairs:
    first_index = find_range_index(first_value, first_column_ranges)
    second_index = find_range_index(second_value - first_value, second_column_ranges)
    
    if first_index is not None and second_index is not None:
        counts_matrix[second_index][first_index] += 1

# Print the counts for debugging purposes
for i, row in enumerate(counts_matrix):
    for j, count in enumerate(row):
        second_range = second_column_ranges[i]
        first_range = first_column_ranges[j]
        print(f"Region {first_range} & {second_range}: {count} members")

# Create the heatmap
fig, ax = plt.subplots()
cax = ax.matshow(counts_matrix, cmap='Blues')

# Add colorbar
plt.colorbar(cax)

# Set axis labels
ax.set_xticks(np.arange(len(first_column_ranges)))
ax.set_yticks(np.arange(len(second_column_ranges)))
ax.set_xticklabels([f'[{start}-{end}]' for start, end in first_column_ranges], horizontalalignment='center', fontsize=10)
ax.set_yticklabels([f'[{start}-{end}]' for start, end in second_column_ranges], rotation=90, verticalalignment='center', fontsize=10)


# Invert y-axis to maintain the correct orientation
ax.invert_yaxis()

# Move x-axis ticks to the bottom and remove from the top
ax.tick_params(bottom=True, top=False, labelbottom=True, labeltop=False)

# Annotate each cell with the formatted count value
for i in range(len(second_column_ranges)):
    for j in range(len(first_column_ranges)):
        count = counts_matrix[i][j]
        formatted_count = format_count(count)
        # Get the color value of the rectangle
        color = cax.get_cmap()(cax.norm(count))
        # Determine text color based on rectangle color intensity
        text_color = 'white' if np.mean(color[:3]) < 0.5 else 'black'
        ax.text(j, i, formatted_count, va='center', ha='center', color=text_color, fontsize=16)

# Set labels
ax.set_xlabel('Alignment (bp)', fontsize=16)
ax.set_ylabel('Polishing Improvement (bp)', fontsize=16)

# plt.tight_layout()
# Show the heatmap
plt.show()
plt.savefig(output_file_name)