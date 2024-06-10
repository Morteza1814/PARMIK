import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

def plot_stacked_bar_chart(input_file, output_image):
    # Read the file and process data
    labels = []
    parmik_vs_bwa = []
    parmik_vs_blast = []
    
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            label = ' '.join(parts[:-2]).strip('"')
            parmik_bwa_value = int(parts[-2])
            parmik_blast_value = int(parts[-1])
            
            labels.append(label)
            parmik_vs_bwa.append(parmik_bwa_value)
            parmik_vs_blast.append(parmik_blast_value)

    # Convert lists to numpy arrays for plotting
    labels = np.array(labels)
    parmik_vs_bwa = np.array(parmik_vs_bwa)
    parmik_vs_blast = np.array(parmik_vs_blast)

    sum_parmik_vs_bwa = parmik_vs_bwa.sum()
    sum_parmik_vs_blast = parmik_vs_blast.sum()
    
    ind = np.arange(2)  # two groups
    
    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Colors for each stacked part
    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in np.linspace(0, 1, len(labels))]
    
    # Plot the bars
    bottom_bwa = np.zeros(1)
    bottom_blast = np.zeros(1)
    
    bars = []
    for i in range(len(labels)):
        bar1 = ax.bar(ind[0], parmik_vs_bwa[i], width=0.35, bottom=bottom_bwa, label=labels[i], color=colors[i])
        bar2 = ax.bar(ind[1], parmik_vs_blast[i], width=0.35, bottom=bottom_blast, color=colors[i])
        bars.append(bar1)
        bottom_bwa += parmik_vs_bwa[i]
        bottom_blast += parmik_vs_blast[i]
    
    
    ax.tick_params(axis='y', labelsize=16)  # Increase font size of y-tick labels
    ax.set_ylabel('No. of Queries', fontsize=22)
    # ax.set_title('Comparison of PARMIK vs BWA and PARMIK vs BLAST')
    ax.set_xticks(ind)
    ax.set_xticklabels(['PARMIK vs. BWA', 'PARMIK vs. BLAST'], fontsize=18)
    
    # Formatter for y-ticks to display in kilo format
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{int(x/1000)}K' if x >= 1000 else f'{int(x)}'))
    
    # Add legend
    ax.legend([bar[0] for bar in bars], labels)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(output_image, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_image>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_image = sys.argv[2]
    
    plot_stacked_bar_chart(input_file, output_image)
