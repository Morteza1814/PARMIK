import matplotlib.pyplot as plt

# Initialize variables to count the number of better alignments for each tool
tool1_better = 0
tool2_better = 0

# Initialize dictionary to store the differences in alignment sizes
alignment_differences = {}

# Prompt the user for the file address
file_address = input("Enter the file address containing alignment sizes: ")

# Prompt the user for the name of the second tool (BLAST)
second_tool_name = input("Enter the name of the second tool (BLAST): ")


# Read alignment sizes from input file
with open(file_address, 'r') as f:
    for i in range(-150, 151):
        alignment_differences[i] = 0
    for line in f:
        values = line.strip().split()
        if len(values) != 2:
            continue
        if values[0] == '-' or values[1] == '-':
            continue
        size1, size2 = map(int, values)
        difference = size1 - size2
        
        if size1 > size2:
            tool1_better += 1
        elif size2 > size1:
            tool2_better += 1

        if difference in alignment_differences:
            alignment_differences[difference] += 1
        else:
            alignment_differences[difference] = 1

# Sort alignment differences dictionary by keys
sorted_alignment_differences = dict(sorted(alignment_differences.items()))

# Plotting
x = list(sorted_alignment_differences.keys())
y = list(sorted_alignment_differences.values())

color_labels = {
    'red': 'BLAST outperfroms',
    'blue': 'Equal Performance',
    'green': 'PARMIK outperforms'
}

colors = ['blue' if val == 0 else 'red' if val <= 0 else 'green' for val in x]
labels = [color_labels[color] for color in colors]

plt.bar(x, y, color=colors, label='Same alignment size')

# Create custom legend
handles = [plt.Rectangle((0,0),1,1, color=color) for color in color_labels.keys()]
plt.legend(handles, color_labels.values())

plt.xlabel('Difference in Alignment Size or Edit Distance')
plt.ylabel('Number of Matches')
plt.title('PARMIK vs BLAST Alignment Sizes')

# Add annotations for better alignments
plt.annotate(f'{tool1_better} queries', xy=(max(x)//2, max(y)//2), xytext=(max(x)//2, max(y)//2 + 1), ha='center', color='green')
plt.annotate(f'{tool2_better} queries', xy=(-max(x)//2, max(y)//2), xytext=(-max(x)//2, max(y)//2 + 1), ha='center', color='red')

# Save the plot in a file
plt.savefig('genome_alignment_comparison_plot.png')

# Show the plot
plt.show()
