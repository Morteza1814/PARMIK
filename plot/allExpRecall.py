import sys
import matplotlib.pyplot as plt

# Function to parse the recall data from the input file
def parse_recall_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extracting experiment names (first line)
    experiments = lines[0].strip().split('\t\t')

    # Extracting recall rates (third line)
    recall_values = list(map(float, lines[2].strip().split('\t')))

    # Separating recall rates for PARMIK and BLAST
    parmik_recall = [recall_values[i] for i in range(0, len(recall_values), 2)]
    blast_recall = [recall_values[i+1] for i in range(0, len(recall_values), 2)]

    return experiments, parmik_recall, blast_recall

# Main function to run the script
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_file> <output_file>>")
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    experiments, parmik_recall, blast_recall = parse_recall_data(in_file)

    # Plotting the bar diagram
    experiments_formatted = [exp.replace(' vs. ', '\nvs.\n') for exp in experiments]
    width = 0.25  # width of the bars
    x = range(len(experiments_formatted))

    plt.figure(figsize=(10, 6))
    plt.bar([i - width/2 for i in x], parmik_recall, width=width, label='PARMIK', color='navy')
    plt.bar([i + width/2 for i in x], blast_recall, width=width, label='BLAST', color='firebrick')

    # Enlarging y-axis label and x, y ticks
    plt.ylabel('Recall Rate (%)', fontsize=20)
    plt.xticks(ticks=x, labels=experiments_formatted, fontsize=16)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)

    # Removing the figure title and grid
    plt.gca().xaxis.label.set_visible(False)
    plt.grid(False)

    plt.tight_layout()
    plt.savefig(out_file)
