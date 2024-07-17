import matplotlib.pyplot as plt
import sys
import numpy as np

output_image = sys.argv[1]

# Data
labels = ['0 - 20', '21 - 150']
sizes = [1977411150, 19504003]
total = sum(sizes)

# Calculate percentages
percentages = [size / total * 100 for size in sizes]

# Narrow stacked bar plot
fig, ax = plt.subplots()

# Plot the narrow stacked bar
width = 0.01
x = [0]

ax.bar(x, percentages[0], width, label='0 - 20', color='skyblue')
ax.bar(x, percentages[1], width, bottom=percentages[0], label='21 - 150', color='lightgrey')

# Add text for percentages
ax.text(x[0], percentages[0] / 2, f'{percentages[0]:.2f}%', ha='center', va='center', color='black', fontsize=14)
ax.text(x[0], percentages[0] + percentages[1] / 2, f'{percentages[1]:.2f}%', ha='center', va='center', color='black', fontsize=14) 

# Title and labels
# ax.set_title('Percentage of Total by Range')
ax.set_ylabel('Alignment Size Percentage (%)', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(['Wuhan-Hu-1 vs. SRR12432009'], fontsize=14)

plt.setp(ax.get_yticklabels(), fontsize=14)
# Add legend
ax.legend(fontsize=14)

plt.savefig(output_image)
