import matplotlib.pyplot as plt
import sys
import numpy as np

output_image = sys.argv[1]

# Data
labels = ['0 - 20', '21 - 150']
data = {
    'Wuhan-Hu-1 \nvs. \nSRR12432009': [1977411150, 19504003],
    'Wuhan-Hu-1 \nvs. \nSRR14381434': [1252240164, 5656029],
    'SarsGenome \nvs. \nSRR12432009': [1697824216, 10702742],
    'SarsGenome \nvs. \nSRR14381434': [1279545591, 3853695],
}

# Calculate percentages for each bar
percentages = {key: [size / sum(value) * 100 for size in value] for key, value in data.items()}

# Narrow stacked bar plot
fig, ax = plt.subplots()

# Plot the narrow stacked bars
width = 0.4
x = np.arange(len(data))

for i, (key, value) in enumerate(percentages.items()):
    ax.bar(x[i], value[0], width, label='0 - 20' if i == 0 else "", color='skyblue')
    ax.bar(x[i], value[1], width, bottom=value[0], label='21 - 150' if i == 0 else "", color='red')

    # Add text for percentages
    # ax.text(x[i], value[0] / 2, f'{value[0]:.2f}%', ha='center', va='center', color='black', fontsize=14)
    # ax.text(x[i], value[0] + value[1] / 2, f'{value[1]:.2f}%', ha='center', va='center', color='black', fontsize=14)

# Title and labels
ax.set_ylabel('Alignment Size (%)', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(data.keys(), fontsize=14, rotation=45, ha='center')

plt.setp(ax.get_yticklabels(), fontsize=14)
# Add legend
ax.legend(fontsize=12, bbox_to_anchor=(0.85, .75), loc='center')

plt.tight_layout()

plt.savefig(output_image)
