import pandas as pd
import matplotlib.pyplot as plt
import sys

output_file = sys.argv[1]
# Data preparation
datasets = [
    "Wuhan-Hu-1 vs. SRR12432009",
    "SarsGenome vs. SRR12432009",
    "Wuhan-Hu-1 vs. SRR14381418",
    "SarsGenome vs. SRR14381418",
    "Wuhan-Hu-1 vs. SRR14381422",
    "SarsGenome vs. SRR14381422",
    "Wuhan-Hu-1 vs. SRR14381425",
    "SarsGenome vs. SRR14381425",
]

# Data for each R value
data_R20 = {
    "datasets": datasets,
    "PARMIK IKT0+Polishing": [55.75292153, 54.11686794, 58.56529963, 52.59595631, 32.78385713, 0, 45.03342278, 31.39805896],
    "PARMIK IKT5000+Polishing": [49.64803852, 50.11007477, 54.66315307, 50.92667077, 28.14434258, 0, 22.36153924, 17.61697452],
    "BLAST": [35.36943983, 20.10555709, 20.55594172, 12.33344455, 1.631199671, 0, 0.7196505234, 0.4904444938]
}

data_R25 = {
    "datasets": datasets,
    "PARMIK IKT0+Polishing": [98.38093052, 96.50427376, 95.36952138, 94.74096345, 80.94908302, 0, 98.03892284, 96.4476342],
    "PARMIK IKT5000+Polishing": [89.02840505, 93.86112245, 94.51615055, 94.33153348, 69.29804234, 0, 3.720593268, 14.19277189],
    "BLAST": [84.5726937, 66.69167501, 60.03873669, 46.91984642, 29.97161647, 0, 50.41945759, 55.19129307]
}

data_R30 = {
    "datasets": datasets,
    "PARMIK IKT0+Polishing": [99.14989974, 97.69927326, 97.06250289, 96.95812982, 86.1822551, 0, 99.68246641, 96.1198819],
    "PARMIK IKT5000+Polishing": [93.57272954, 97.28001911, 96.65432061, 96.78509413, 79.15427524, 0, 2.598533445, 31.08621366],
    "BLAST": [89.44510742, 71.4221907, 64.43928323, 49.76253448, 36.25012331, 0, 66.59704421, 82.73610604]
}

df_R20 = pd.DataFrame(data_R20)
df_R25 = pd.DataFrame(data_R25)
df_R30 = pd.DataFrame(data_R30)

# Plotting with specified colors for BLAST and PARMIK
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Colors for the bars
colors_R20 = ['navy', 'dodgerblue', 'firebrick']
colors_R25 = ['navy', 'dodgerblue', 'firebrick']
colors_R30 = ['navy', 'dodgerblue', 'firebrick']

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
    ax.set_xticklabels(new_labels, fontsize=14, rotation=45)

plt.tight_layout()
plt.savefig(output_file)
