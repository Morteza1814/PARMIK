import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

def read_first_and_last_column(file_path):
    try:
        df = pd.read_csv(file_path, sep='\s+', header=None)
        first_column = df.iloc[:, 0].to_numpy()
        last_column = df.iloc[:, -1].to_numpy()
        return first_column, last_column
    except Exception as e:
        print(f"An error occurred: {e}")

def calculate_recall_rate(tpArray, fnArray):
    if tpArray.shape != fnArray.shape:
        raise ValueError("Both arrays must be of the same size.")
    tpSum = np.sum(tpArray)
    fnSum = np.sum(fnArray)
    return (tpSum / (tpSum + fnSum)) * 100 if (tpSum + fnSum) > 0 else 0

def getRecallOutOfTpFn(tpFilePath, fnFilePath, r):
    tp = read_first_and_last_column(tpFilePath)
    fn = read_first_and_last_column(fnFilePath)
    if tp is not None and fn is not None:
        tp_array = tp[1]
        fn_array = fn[1]
        if r - 1 < len(tp_array) and r - 1 < len(fn_array):
            return calculate_recall_rate(tp_array[r-1:], fn_array[r-1:])

def getTime(timeFilePath):
    with open(timeFilePath, 'r') as file:
        return int(file.readline().strip())

def getSpeedOutOfTp(tpFilePath, timeFilePath, r):
    with open(timeFilePath, 'r') as file:
        time = int(file.readline().strip())
    tp = read_first_and_last_column(tpFilePath)
    if tp is not None:
        tp_array = tp[1]
        if r - 1 < len(tp_array) and time > 0:
            return np.sum(tp_array[r-1:]) / (60 * time)

# Experiment settings
metagenomicDatasets = ['SRR12432009', 'SRR14381418', 'SRR14381421', 'SRR14381422']
queryDatasets = ['NC_045512', 'SarsGenome']
IKT = 0
M = [7, 4]
PI = [85, 90, 95]
T = 1
SC = 1
K = 11
R = 30

# Pair map for labeling dots
pair_list = list(product(queryDatasets, metagenomicDatasets))
pair_map = {f"{q}-{m}": i+1 for i, (q, m) in enumerate(pair_list)}

# Plot setup
rows = 1
cols = 3
width = 10
height = 3
fig, axs = plt.subplots(rows, cols, figsize=(width, height))

experimentPath = sys.argv[1]

for metagenomicDataset in metagenomicDatasets:
    for queryDataset in queryDatasets:
        pair_id = f"{queryDataset}-{metagenomicDataset}"
        pair_number = pair_map[pair_id]
        blast_speed, blast_time, blast_recall = [], [], []
        bowtie_speed, bowtie_time, bowtie_recall = [], [], []

        for m in M:
            parmik_recall, parmik_speed, parmik_time = [], [], []
            for i, pi in enumerate(PI):
                ax = axs[i]
                dir = f"{experimentPath}/{metagenomicDataset}/{queryDataset}/IKT{IKT}_K{K}_PI{pi}_M{m}_T{T}_SC{SC}_P_2484_1111_2444_2288_2848/"
                parmikTpFilePath = dir + 'cmp_Baseline_parmik_AlnSz_tp.txt'
                parmikFnFilePath = dir + 'cmp_Baseline_parmik_AlnSz_fn.txt'
                parmikTimePath = dir + 'time_parmik'
                recall_val = getRecallOutOfTpFn(parmikTpFilePath, parmikFnFilePath, R)
                time_val = getTime(parmikTimePath)
                parmik_recall.append(recall_val)
                parmik_speed.append(getSpeedOutOfTp(parmikTpFilePath, parmikTimePath, R))
                parmik_time.append(time_val)
                if m == 7:
                    dx = 0
                    dy = 0
                    if pi == 90:
                        if pair_number == 2:
                            dx = 2
                        if pair_number == 6:
                            dx = -2
                        if pair_number == 7:
                            dx = -3
                            dy=10
                    if pi == 95:
                        if pair_number == 3 or pair_number == 4:
                            dx = -8
                            dy = 10
                        if pair_number == 7:
                            dy = 15
                            dx = -3
                        if pair_number == 6:
                            dy = -25
                        if pair_number == 5:
                            dy = 5
                        if pair_number == 1:
                            dy = 10

                    ax.scatter(recall_val + dx, time_val + dy, label='PARMIK [M=7]', color='dodgerblue', s=90)
                    ax.text(recall_val + dx, time_val + dy, str(pair_number), fontsize=9, ha='center', va='center', color='white')
                elif m == 4:
                    dx = 0
                    dy = 0
                    if pi == 85:
                        if pair_number == 2 or pair_number == 7 or pair_number == 8:
                            dy = -40
                        if pair_number == 6:
                            dy = -70
                    if pi == 90:
                        if pair_number == 2:
                            dy = -30
                        if pair_number == 6 or pair_number == 8:
                            dy = -60
                        if pair_number == 7:
                            dx = -5
                    if pi == 95:
                        if pair_number == 3:
                            dy = 80
                        if pair_number == 4:
                            dy = 20
                        if pair_number == 2:
                            dy = -80
                        if pair_number == 7 or pair_number == 8:
                            dy = -120
                        if pair_number == 5:
                            dy = 50
                        if pair_number == 1:
                            dy = 20

                        # if pair_number == 6 or pair_number == 8:
                        #     dy = -60
                        # if pair_number == 7:
                        #     dx = -5
                    ax.scatter(recall_val + dx, time_val + dy, label='PARMIK [M=4]', color='navy', s=90)
                    ax.text(recall_val + dx, time_val + dy, str(pair_number), fontsize=9, ha='center', va='center', color='white')

                    blastTpFilePath = dir + 'cmp_Baseline_blast_AlnSz_tp.txt'
                    blastFnFilePath = dir + 'cmp_Baseline_blast_AlnSz_fn.txt'
                    blastTimePath = dir + 'time_blast'
                    recall_val = getRecallOutOfTpFn(blastTpFilePath, blastFnFilePath, R)
                    time_val = getTime(blastTimePath)
                    blast_recall.append(recall_val)
                    blast_speed.append(getSpeedOutOfTp(blastTpFilePath, blastTimePath, R))
                    blast_time.append(time_val)
                    dx = 0
                    dy = 0
                    if pi == 85: 
                        if pair_number == 3 or pair_number == 7:
                            dy = 20   
                    if pi == 95:
                        if pair_number == 1 or pair_number == 5:
                            dx = -4       
                    if pi == 90:
                        if pair_number == 8:
                            dy = -10        
                    ax.scatter(recall_val + dx, time_val + dy, label='BLAST', color='firebrick', s=90)
                    ax.text(recall_val + dx, time_val + dy, str(pair_number), fontsize=9, ha='center', va='center', color='white')

                    bowtieTpFilePath = dir + 'cmp_Baseline_bowtie_AlnSz_tp.txt'
                    bowtieFnFilePath = dir + 'cmp_Baseline_bowtie_AlnSz_fn.txt'
                    bowtieTimePath = dir + 'time_bowtie'
                    recall_val = getRecallOutOfTpFn(bowtieTpFilePath, bowtieFnFilePath, R)
                    time_val = getTime(bowtieTimePath)
                    bowtie_recall.append(recall_val)
                    bowtie_speed.append(getSpeedOutOfTp(bowtieTpFilePath, bowtieTimePath, R))
                    bowtie_time.append(time_val)
                    dx = 0
                    dy = 0
                    if pi == 85: 
                        if pair_number == 5:
                            dy = 10
                    if pi == 95:
                        if pair_number == 3:
                            dy = 80
                        if pair_number == 4:
                            dy = 20
                    ax.scatter(recall_val + dx, time_val + dy, label='Bowtie', color='darkgreen', s=90)
                    ax.text(recall_val + dx, time_val + dy, str(pair_number), fontsize=9, ha='center', va='center', color='white')

                ax.set_title(f'PI = {PI[i]}%')

for ax in axs:
    ax.set_xlabel('Recall Rate (%)')
    ax.set_ylabel('Time (Minutes)')
    ax.set_yscale('log')
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

handles, labels = axs[-1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[-1].legend(by_label.values(), by_label.keys(), frameon=False, loc='upper left', bbox_to_anchor=(1, 0.3), fontsize=9)

# Mapping: number -> (query, metagenome)
mapping_lines = [
    "(1) Wuhan-Hu-1 vs. SRR12432009",
    "(2) Wuhan-Hu-1 vs. SRR14381418",
    "(3) Wuhan-Hu-1 vs. SRR14381421",
    "(4) Wuhan-Hu-1 vs. SRR14381422",
    "(5) Sars-1 vs. SRR12432009",
    "(6) Sars-1 vs. SRR14381418",
    "(7) Sars-1 vs. SRR14381421",
    "(8) Sars-1 vs. SRR14381422"
]

for i, line in enumerate(mapping_lines):
    y_pos = 0.92 - i * 0.075  # increased spacing from 0.04 → 0.055

    # Extract number and label
    number = line.split(')')[0].replace('(', '')  # e.g., "1"
    label = line.split(')')[1].strip()            # e.g., "Wuhan-Hu-1 vs. SRR12432009"

    # Draw black circle with white number inside
    axs[-1].text(1.14, y_pos, number, transform=axs[-1].transAxes,
                 fontsize=9, ha='center', va='center',
                 color='white', bbox=dict(boxstyle='circle,pad=0.1', fc='black', ec='none'))

    # Draw the label next to the circle
    axs[-1].text(1.20, y_pos, label, transform=axs[-1].transAxes,
                 fontsize=9, ha='left', va='center', color='black')

plt.tight_layout()
# Uncomment below and set valid experimentPath when ready to save
plt.savefig(experimentPath + '/allExpFigures/timeVSrecall.pdf')
