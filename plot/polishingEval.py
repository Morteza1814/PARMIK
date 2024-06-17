import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_scatter(input_file, output_file):
    # Read the data from the file using pandas for efficiency
    data = pd.read_csv(input_file, sep=' ', header=None, names=['Alignment', 'Alignment_polishing'])

    # Separate the data into the two cases, excluding the equal ones
    greater = data[data['Alignment_polishing'] > data['Alignment']]
    less = data[data['Alignment_polishing'] < data['Alignment']]
    
    print("reading the data finished!")

    # Create scatter plot
    plt.figure(figsize=(10, 6))

    # Plot each case with different colors
    if not greater.empty:
        plt.scatter(greater['Alignment'], greater['Alignment_polishing'], color='royalblue', label='Alignment+Polishing > Alignment', s=1)
    if not less.empty:
        plt.scatter(less['Alignment'], less['Alignment_polishing'], color='red', label='Alignment+Polishing < Alignment', s=1)

    print("plotting finished!")

    # Add a dotted line y = x
    max_value = max(data['Alignment'].max(), data['Alignment_polishing'].max())
    plt.plot([0, max_value], [0, max_value], 'k--', linewidth=0.5)

    # Set labels and title
    plt.xlabel('Alignment (bp)', fontsize=22)
    plt.ylabel('Alignment+Polishing (bp)', fontsize=22)
    # plt.title('Scatter Plot of Alignment vs Alignment+Second Chance')
    # plt.legend(fontsize=16)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    # Count number of points above and below y=x line
    above_line = greater.shape[0] if not greater.empty else 0
    below_line = less.shape[0] if not less.empty else 0

    # Add text annotations for counts below and above y=x line
    plt.text(max_value * 0.05, max_value * 0.95, f'Count = {above_line}', fontsize=20, color='royalblue')
    plt.text(max_value * 0.75, max_value * 0.05, f'Count = {below_line}', fontsize=20, color='red')

    # Save the plot to the output file
    plt.savefig(output_file, dpi=300)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_scatter(input_file, output_file)
