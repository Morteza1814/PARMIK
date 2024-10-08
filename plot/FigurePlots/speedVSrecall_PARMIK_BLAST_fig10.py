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
M = [4, 7]
PI = [85, 90, 95]
T = 1
SC = 1
K = 11
R = 30

experimentPath = sys.argv[1]
markers = ['o', 'x', 's', 'D', '^']
firstRound=True
for metagenomicDataset in metagenomicDatasets:
    for queryDataset in queryDatasets:
        # Create figure and axis objects without grids
        fig, ax1 = plt.subplots()

        # Increase the font size for clarity
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        ax1.set_xlabel('Percentage Identity (%)', fontsize=20)
        ax1.set_ylabel('Speed (Matches/Second)', fontsize=20, color='tab:blue')
        # Create a second y-axis for Recall Rate (%)
        ax2 = ax1.twinx()
        ax2.set_ylabel('Recall Rate (%)', fontsize=20, color='tab:red')

        if queryDataset == 'NC_045512':
            datasetName = 'Wuhan-Hu-1' + ' vs. ' + metagenomicDataset
        else:
            datasetName = queryDataset + ' vs. ' + metagenomicDataset
        print('=========================================')
        print(datasetName)
        blast_speed = []     
        blast_recall = []
        markerInd = 0
        for m in M:
            parmik_recall = [] 
            parmik_speed = []
            for pi in PI:                
                dir = experimentPath + '/' + metagenomicDataset + '/' + queryDataset + '/IKT' + str(IKT) + '_K' + str(K) + '_PI' + str(pi) + '_M' + str(m) + '_T' + str(T) + '_SC' + str(SC) + '_P_2484_1111_2444_2288_2848/'
                parmikTpFilePath = dir + 'cmp_Baseline_parmik_AlnSz_tp.txt'
                parmikFnFilePath = dir + 'cmp_Baseline_parmik_AlnSz_fn.txt'
                parmikTimePath = dir + 'time_parmik'
                parmik_recall.append(getRecallOutOfTpFn(parmikTpFilePath, parmikFnFilePath, R))
                parmik_speed.append(getSpeedOutOfTp(parmikTpFilePath, parmikTimePath, R))
                if m == 4:
                    blastTpFilePath = dir + 'cmp_Baseline_blast_AlnSz_tp.txt'
                    blastFnFilePath = dir + 'cmp_Baseline_blast_AlnSz_fn.txt'
                    blastTimePath = dir + 'time_blast'
                    blast_recall.append(getRecallOutOfTpFn(blastTpFilePath, blastFnFilePath, R))
                    blast_speed.append(getSpeedOutOfTp(blastTpFilePath, blastTimePath, R))

            ax1.plot(PI, parmik_speed, label=f'PARMIK [M={m}]', marker=markers[markerInd], color='tab:blue')
            ax2.plot(PI, parmik_recall, label=f'PARMIK [M={m}]', marker=markers[markerInd], color='tab:red')
            print('=====>>>>', dir)        
            print('parmik_recall: ', parmik_recall)
            print('parmik_speed: ', parmik_speed)
            markerInd+=1
        print('blast_recall: ', blast_recall)
        print('blast_speed: ', blast_speed)
        ax1.plot(PI, blast_speed, label="BLAST", marker='^', linestyle='--', color='tab:blue')
        ax2.plot(PI, blast_recall, label="BLAST", marker='^', linestyle='--', color='tab:red')
        ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=14)
        ax2.tick_params(axis='y', labelcolor='tab:red', labelsize=14)
        # Set the x-ticks to only show 85, 90, and 95
        ax1.set_xticks(PI)
        ax1.set_xticklabels([85, 90, 95], fontsize=14)
        plt.title(datasetName, fontsize=20)
        # Combining legends
        if firstRound:
            firstRound=False    
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines + lines2, labels + labels2, loc='center left', fontsize=14)

        # Adding legends with adjusted font sizes
        fig.tight_layout()

        # Remove the grid
        ax1.grid(False)
        ax2.grid(False)
        ax2.set_ylim(20, 105)


        # Show the plot
        plt.savefig(experimentPath + '/' + metagenomicDataset + '/' + queryDataset + '/SvR_' + datasetName + '_' + str(R)+ '_ylim.pdf')