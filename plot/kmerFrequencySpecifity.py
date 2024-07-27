import matplotlib.pyplot as plt
import pandas as pd

# Data provided
data = {
    'Range': ['[0-5K]', '[5K-50K]', '[50K-500K]', '[500K-5M]', '[5M-50M]', '[50M-500M]'],
    'k=11': [4167750, 23467, 2561, 513, 11, 1],
    'k=14': [135873777, 22180, 2512, 421, 11, 1],
    'k=17': [262161198, 21527, 2456, 373, 10, 0]
}

# Creating DataFrame
df = pd.DataFrame(data)

# Define bar width and positions
bar_width = 0.25
positions = list(range(len(df['Range'])))

# Adjusting the plot to make the legends and axis ticks larger

fig, ax = plt.subplots(figsize=(12, 8))

# Define different shades of blue for each k-mer size
colors = ['#1f77b4', '#2ca02c', '#ff7f0e']

# Plotting bars for each k-mer size with logarithmic scale
plt.bar([p - bar_width for p in positions], df['k=11'], width=bar_width, align='center', label='k=11', color=colors[0])
plt.bar(positions, df['k=14'], width=bar_width, align='center', label='k=14', color=colors[1])
plt.bar([p + bar_width for p in positions], df['k=17'], width=bar_width, align='center', label='k=17', color=colors[2])

# Adding labels and title
ax.set_xlabel('K-mer Frequency Ranges', fontsize=24)
ax.set_ylabel('Frequency (Log Scale)', fontsize=24)
# ax.set_title('K-mer Frequency Distribution', fontsize=16)
ax.set_xticks(positions)
ax.set_xticklabels(df['Range'], fontsize=12)
ax.set_yscale('log')

# Enlarging the legend
plt.legend(fontsize=22)

# Enlarging the ticks
plt.xticks(fontsize=20, rotation=45)
plt.yticks(fontsize=20)

plt.tight_layout()
# Display the plot
plt.show()

plt.savefig('/p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/kmerSizeSpecifityAnalysis/kmerFrequencySpecifity.png')
