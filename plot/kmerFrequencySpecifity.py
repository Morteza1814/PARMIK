import matplotlib.pyplot as plt
import pandas as pd

# Data provided
data = {
    # 'Range': ['[0-5K]', '[5K-50K]', '[50K-500K]', '[500K-5M]', '[5M-50M]', '[50M-500M]'
    'Range': ['[0-25]', '[26-50]', '[51-100]', '[101-150]', '[151-200]', '[201-250]', '[251-500]', '[501-750]', '[751-1000]', '[1001-1500]', '[1501-2000]'],
    'k=11': [110035, 457116, 1079833, 769764, 494205, 319901, 609191, 156160, 59664, 51632, 24939],
    'k=14': [131569770, 2227792, 1060533, 357538, 167659, 101340, 210784, 61759, 30754, 36609, 17939],
    'k=17': [258453785, 1762322, 952797, 343622, 164668, 100141, 203463, 60978, 31907, 37652, 17704],
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

plt.savefig('/p/lava/parmik/partialMatcherProject/tools/parmik/CK_INDEX/SRR12432009/kmerSizeSpecifityAnalysis/kmerFrequencySpecifity.pdf')
