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
        
def plot_recall_rate(data_collection, output_file, r_list):
    data_collection_size = len(data_collection)
    fig, axs = plt.subplots(data_collection_size, 1, figsize=(10, data_collection_size * 5))
    if not isinstance(axs, np.ndarray):
        axs = [axs] 
    i = 0
    for data in data_collection:
        df = pd.DataFrame(data)

        # Plotting with specified colors for BLAST and PARMIK

        # Colors for the bars
        colors = ['navy', 'dodgerblue', 'firebrick']

        # Plot for R=20
        df.plot(x='datasets', kind='bar', ax=axs[i], color=colors)
        # axs[i].set_title('R=' + str(R[i]), fontsize=18)
        axs[i].set_xlabel('', fontsize=12)
        axs[i].set_ylabel('Recall Rate (%)', fontsize=16)
        axs[i].tick_params(axis='both', which='major', labelsize=14)
        if i ==0:
            axs[i].legend(fontsize=14)
        else:
            axs[i].legend().set_visible(False)
        axs[i].set_ylim(0, 100)
        axs[i].set_yticks(range(0, 101, 20))
        i+=1

    # Adjusting the layout
    fig.tight_layout()

    # Adjust x-tick labels to be on two lines
    for ax in axs:
        labels = ax.get_xticklabels()
        new_labels = ['\n'.join(label.get_text().split(' vs. ')) for label in labels]
        ax.set_xticklabels(new_labels, fontsize=14, rotation=45)

    plt.tight_layout()
    plt.savefig(output_file)


experimentPath = sys.argv[1]

metagenomicDatasets = ['SRR12432009', 'SRR14381418', 'SRR14381421', 'SRR14381422']
queryDatasets = ['NC_045512', 'SarsGenome']
IKT = [0,5000]
M = 4
PI = 90
T = 1
SC = 1
K = 11
R = [30]

data_collection = []

for r in R:
    datasets = []
    parmik_ikt_all_polishing = []
    parmik_ikt_5000_polishing = []
    blast = []
    for metagenomicDataset in metagenomicDatasets:
        for queryDataset in queryDatasets:
            if queryDataset == 'NC_045512':
                datasets.append('Wuhan-Hu-1' + ' vs. ' + metagenomicDataset)
            else:
                datasets.append(queryDataset + ' vs. ' + metagenomicDataset)
            for ikt in IKT:
                dir = experimentPath + '/' + metagenomicDataset + '/' + queryDataset + '/IKT' + str(ikt) + '_K' + str(K) + '_PI' + str(PI) + '_M' + str(M) + '_T' + str(T) + '_SC' + str(SC) + '_P_2484_1111_2444_2288_2848/'
                parmikTpFilePath = dir + 'cmp_Baseline_parmik_AlnSz_tp.txt'
                parmikFnFilePath = dir + 'cmp_Baseline_parmik_AlnSz_fn.txt'
                if ikt == 5000: 
                    parmik_ikt_5000_polishing.append(getRecallOutOfTpFn(parmikTpFilePath, parmikFnFilePath, r))
                else:
                    parmik_ikt_all_polishing.append(getRecallOutOfTpFn(parmikTpFilePath, parmikFnFilePath, r))
                    blastTpFilePath = dir + 'cmp_Baseline_blast_AlnSz_tp.txt'
                    blastFnFilePath = dir + 'cmp_Baseline_blast_AlnSz_fn.txt'
                    blast.append(getRecallOutOfTpFn(blastTpFilePath, blastFnFilePath, r))
    data = {
        "datasets": datasets,
        "PARMIK [IKT=all, Polishing]": parmik_ikt_all_polishing,
        "PARMIK [IKT=5000, Polishing]": parmik_ikt_5000_polishing,
        "BLAST": blast
    }   
    data_collection.append(data)

print(data_collection)
output_file = experimentPath + '/allExpFigures/allExpRecall_30.pdf'
plot_recall_rate(data_collection, output_file, R)

