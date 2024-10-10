import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def read_first_and_last_column(file_path):
    try:
        # Read the space-separated file into a DataFrame
        df = pd.read_csv(file_path, sep='\s+', header=None)  # header=None if there is no header row

        # Get the first and last columns
        first_column = df.iloc[:, 0].to_numpy()
        last_column = df.iloc[:, -1].to_numpy()

        return first_column, last_column

    except Exception as e:
        print(f"An error occurred: {e}")

def calculate_recall_rate(tpArray, fnArray):
    # Check if both arrays are of the same size
    if tpArray.shape != fnArray.shape:
        raise ValueError("Both arrays must be of the same size.")

    # Create an array to store the recall rates
    tpSum = np.sum(tpArray)
    fnSum = np.sum(fnArray)
    
    denominator = tpSum + fnSum
    if denominator > 0:
        recall_rate_array = (tpSum / denominator) * 100
    else:
        recall_rate_array = 0  # Define recall rate as 0 if both are zero

    return recall_rate_array

def getRecallOutOfTpFn(tpFilePath, fnFilePath, r):
    tp = read_first_and_last_column(tpFilePath,)
    fn = read_first_and_last_column(fnFilePath)
    if tp is not None and fn is not None: 
        tp_array = tp[1]  # Second array from tp tuple
        fn_array = fn[1]  # Second array from fn tuple
        if r - 1 < len(tp_array) and r - 1 < len(fn_array):
            recall = calculate_recall_rate(tp_array[r-1:], fn_array[r-1:])
            return recall
        
def getTime(timeFilePath):
    with open(timeFilePath, 'r') as file:
        time = int(file.readline().strip())
        return time

def getSpeedOutOfTp(tpFilePath, timeFilePath, r):
    with open(timeFilePath, 'r') as file:
        time = int(file.readline().strip())
    tp = read_first_and_last_column(tpFilePath)
    if tp is not None: 
        tp_array = tp[1]  # Second array from tp tuple
        if r - 1 < len(tp_array) and time > 0:
            speed = np.sum(tp_array[r-1:]) / (60 * time)
            return speed

metagenomicDatasets = ['SRR12432009', 'SRR14381418', 'SRR14381421', 'SRR14381422']
queryDatasets = ['NC_045512', 'SarsGenome']
IKT = 0
M = [7,4]
PI = [85, 90, 95]
T = 1
SC = 1
K = 11
R = 30

rows = 1
cols = 3
width = 10
height = 3
fig, axs = plt.subplots(rows, cols, figsize=(width, height))

experimentPath = sys.argv[1]
markers = ['o', 'x', 's', 'D', '^']
firstRound=True
for metagenomicDataset in metagenomicDatasets:
    for queryDataset in queryDatasets:
        
        if queryDataset == 'NC_045512':
            datasetName = 'Wuhan-Hu-1' + ' vs. ' + metagenomicDataset
        else:
            datasetName = queryDataset + ' vs. ' + metagenomicDataset
        print('=========================================')
        print(datasetName)
        blast_speed = []  
        blast_time = []     
        blast_recall = []
        print('PI:', PI)
        for m in M:
            parmik_recall = [] 
            parmik_speed = []
            parmik_time = []
            for i, pi in enumerate(PI):
                ax = axs[i]
                dir = experimentPath + '/' + metagenomicDataset + '/' + queryDataset + '/IKT' + str(IKT) + '_K' + str(K) + '_PI' + str(pi) + '_M' + str(m) + '_T' + str(T) + '_SC' + str(SC) + '_P_2484_1111_2444_2288_2848/'
                parmikTpFilePath = dir + 'cmp_Baseline_parmik_AlnSz_tp.txt'
                parmikFnFilePath = dir + 'cmp_Baseline_parmik_AlnSz_fn.txt'
                parmikTimePath = dir + 'time_parmik'
                parmik_recall.append(getRecallOutOfTpFn(parmikTpFilePath, parmikFnFilePath, R))
                parmik_speed.append(getSpeedOutOfTp(parmikTpFilePath, parmikTimePath, R))
                parmik_time.append(getTime(parmikTimePath))
                if m == 7:
                    ax.scatter(parmik_recall[i], parmik_time[i], label='PARMIK [M=' + str(m) + ']', color='green')
                elif m == 4:
                    ax.scatter(parmik_recall[i], parmik_time[i], label='PARMIK [M=' + str(m) + ']', color='navy')
                    blastTpFilePath = dir + 'cmp_Baseline_blast_AlnSz_tp.txt'
                    blastFnFilePath = dir + 'cmp_Baseline_blast_AlnSz_fn.txt'
                    blastTimePath = dir + 'time_blast'
                    blast_recall.append(getRecallOutOfTpFn(blastTpFilePath, blastFnFilePath, R))
                    blast_speed.append(getSpeedOutOfTp(blastTpFilePath, blastTimePath, R))
                    blast_time.append(getTime(blastTimePath))
                    ax.scatter(blast_recall[i], blast_time[i], label='BLAST', color='firebrick')
                ax.title.set_text(f'PI = {PI[i]}%')

            print('PARMIK M =', m) 
            print('parmik_recall: ', parmik_recall)
            print('parmik_time (Minutes): ', parmik_time)
            print('parmik_speed (Matches/Sec): ', parmik_speed)
            print('---------------------------------------')
        print('BLAST:')
        print('blast_recall: ', blast_recall)
        print('blast_time (Matches/Sec): ', blast_time)
        print('blast_speed: ', blast_speed)
        

for ax in axs:
    ax.set_xlabel('Recall Rate (%)')
    ax.set_ylabel('Time (Minutes)')
    ax.set_xticks([0, 25, 50, 75, 100])  # Set x-axis ticks to 25, 50, 75, and 100
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

ax = axs[-1]
# keep only unique lables in legend
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()
plt.savefig(experimentPath + '/allExpFigures/timeVSrecall.pdf')