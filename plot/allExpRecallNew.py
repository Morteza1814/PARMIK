import pandas as pd
import matplotlib.pyplot as plt
import sys

output_file = sys.argv[1]
# Data preparation
datasets = [
    "Wuhan-Hu-1 vs. SRR12432009",
    "Wuhan-Hu-1 vs. SRR14381434",
    "SarsGenome vs. SRR12432009",
    "SarsGenome vs. SRR14381434"
]

# Data for each R value
data_R20 = {
    "datasets": datasets,
    "IKT5000 (M=4)": [52.80095746, 37.25102095, 53.07612692, 34.49456604],
    "IKT5000 (M=3)": [62.54326518, 49.45897665, 62.84850645, 50.39463211],
    "IKT5000 (M=2)": [0, 0, 0, 70.36842839],
    "BLAST": [38.66163316, 6.161500533, 23.05702906, 6.8264615]
}

data_R25 = {
    "datasets": datasets,
    "IKT5000 (M=4)": [89.70509073, 73.65415134, 95.08374426, 31.6410791],
    "IKT5000 (M=3)": [90.29823846, 74.93916996, 95.69130276, 58.68833887],
    "IKT5000 (M=2)": [0, 0, 0, 93.59675425],
    "BLAST": [87.53722366, 14.98524763, 81.39740709, 53.90948672]
}

data_R30 = {
    "datasets": datasets,
    "IKT5000 (M=4)": [93.98951691, 67.06122773, 98.22773907, 27.9798201],
    "IKT5000 (M=3)": [94.26968699, 69.0871133, 98.62382188, 52.79704134],
    "IKT5000 (M=2)": [0, 0, 0, 98.13914856],
    "BLAST": [92.31669621, 19.64922576, 87.6407173, 89.12914712]
}

df_R20 = pd.DataFrame(data_R20)
df_R25 = pd.DataFrame(data_R25)
df_R30 = pd.DataFrame(data_R30)

# Plotting with specified colors for BLAST and PARMIK
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Colors for the bars
colors_R20 = ['skyblue', 'dodgerblue', 'navy', 'firebrick']
colors_R25 = ['skyblue', 'dodgerblue', 'navy', 'firebrick']
colors_R30 = ['skyblue', 'dodgerblue', 'navy', 'firebrick']

# Plot for R=20
df_R20.plot(x='datasets', kind='bar', ax=axs[0], color=colors_R20)
axs[0].set_title('R=20', fontsize=18)
axs[0].set_xlabel('', fontsize=12)
axs[0].set_ylabel('Recall Rate (%)', fontsize=16)
axs[0].tick_params(axis='both', which='major', labelsize=14)
axs[0].legend(fontsize=14, ncol=1)
axs[0].set_ylim(0, 100)
axs[0].set_yticks(range(0, 101, 20))


# Plot for R=25
df_R25.plot(x='datasets', kind='bar', ax=axs[1], color=colors_R25)
axs[1].set_title('R=25', fontsize=18)
axs[1].set_xlabel('', fontsize=12)
axs[1].set_ylabel('Recall Rate (%)', fontsize=16)
axs[1].tick_params(axis='both', which='major', labelsize=14)
# axs[1].legend(fontsize=14, ncol=1)
axs[1].get_legend().remove()
axs[1].set_ylim(0, 100)
axs[1].set_yticks(range(0, 101, 20))

# Plot for R=30
df_R30.plot(x='datasets', kind='bar', ax=axs[2], color=colors_R30)
axs[2].set_title('R=30', fontsize=18)
axs[2].set_xlabel('', fontsize=12)
axs[2].set_ylabel('Recall Rate (%)', fontsize=16)
axs[2].tick_params(axis='both', which='major', labelsize=14)
# axs[2].legend(fontsize=14, ncol=1)
axs[2].get_legend().remove()
axs[2].set_ylim(0, 100)
axs[2].set_yticks(range(0, 101, 20))

# Adjusting the layout
fig.tight_layout()

# Adjust x-tick labels to be on two lines
for ax in axs:
    labels = ax.get_xticklabels()
    new_labels = ['\n'.join(label.get_text().split(' vs. ')) for label in labels]
    ax.set_xticklabels(new_labels, fontsize=14, rotation=0)

plt.tight_layout()
plt.savefig(output_file)
