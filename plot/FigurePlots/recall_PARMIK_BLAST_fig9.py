import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def plotFig(x, y2, y3, output_file):
    second_ax_name = "PARMIK"
    third_ax_name = "BLAST"

    # Create a figure and axes for both plots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    bar_width = 0.2

    # Plot using dots
    ax1.plot(x, y2, '-', label=second_ax_name, color='navy', markersize=10)
    ax1.plot(x, y3, '-', label=third_ax_name, color='firebrick', markersize=10)

    ax1.legend(fontsize=18, loc='lower right')
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.set_ylabel('Recall Rate (%)', fontsize=20)
    ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

    # Specify the range for the zoomed-in plot
    zoom_start = int(21 - x[0])
    zoom_end = int(35 - x[0])

    print("Zoom Start:", zoom_start, " Zoom End:", zoom_end)

    # Plot the zoomed-in portion as a bar chart in the bottom subplot
    ax2.bar(x[zoom_start:zoom_end], y2[zoom_start:zoom_end], label=second_ax_name, color='navy', width=bar_width)
    ax2.bar([i + bar_width for i in x[zoom_start:zoom_end]], y3[zoom_start:zoom_end], label=third_ax_name, color='firebrick', width=bar_width)  # Shift x values for Dataset 3
    ax2.legend(fontsize=18, loc='lower right')
    ax2.tick_params(axis='both', which='major', labelsize=18)

    # Set y-axis limits to 50-100 for the zoomed-in plot
    ax2.set_ylim(0, 100)

    # Format the y-axis ticks
    # formatter = ticker.FuncFormatter(lambda x, pos: '0' if x == 0 else '{:.0f}k'.format(x*1e-3))
    # ax2.yaxis.set_major_formatter(formatter)
    ax2.set_ylabel('Recall Rate (%)', fontsize=20)
    ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

    plt.xlabel('Minimum Alignment Size (R)', fontsize=20)

    # Show the plot
    plt.tight_layout()
    plt.show()

    plt.savefig(output_file)

def cumulative_sum_from_each_element(arr):
    # Create an array to store the cumulative sums
    cumulative_sum_array = np.zeros_like(arr)

    # Calculate cumulative sums
    for i in range(len(arr)):
        cumulative_sum_array[i] = np.sum(arr[i:])

    return cumulative_sum_array

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
    recall_rate_array = np.zeros_like(tpArray, dtype=float)

    # Calculate recall rates
    for i in range(len(tpArray)):
        # Avoid division by zero
        denominator = tpArray[i] + fnArray[i]
        if denominator > 0:
            recall_rate_array[i] = (tpArray[i] / denominator) * 100
        else:
            recall_rate_array[i] = 0  # Define recall rate as 0 if both are zero

    return recall_rate_array

def getRecallOutOfTpFn(tpFilePath, fnFilePath):
    aln_len, tp = read_first_and_last_column(tpFilePath)
    aln_len, fn = read_first_and_last_column(fnFilePath)
    tp_cum = cumulative_sum_from_each_element(tp)
    fn_cum = cumulative_sum_from_each_element(fn)
    recall = calculate_recall_rate(tp_cum, fn_cum)
    np.set_printoptions(suppress=True)
    # print(recall)
    return aln_len, recall


def trim_arrays_within_range(array1, lower_limit, upper_limit):
    # Create a mask for the elements that fall within the specified range
    mask = (array1 >= lower_limit) & (array1 <= upper_limit)

    # Trim the arrays using the mask
    trimmed_array1 = array1[mask]

    return trimmed_array1

experimentPath = sys.argv[1]
parmikTpFilePath = experimentPath + '/cmp_Baseline_parmik_AlnSz_tp.txt'
parmikFnFilePath = experimentPath + '/cmp_Baseline_parmik_AlnSz_fn.txt'
blastTpFilePath = experimentPath + '/cmp_Baseline_blast_AlnSz_tp.txt'
blastFnFilePath = experimentPath + '/cmp_Baseline_blast_AlnSz_fn.txt'
output_file = experimentPath + '/PARMIKvsBLAST_recall.pdf'

aln_len, parmikRecall = getRecallOutOfTpFn(parmikTpFilePath, parmikFnFilePath)
aln_len, blastRecall = getRecallOutOfTpFn(blastTpFilePath, blastFnFilePath)

plotFig(aln_len[10:150], parmikRecall[10:150], blastRecall[10:150], output_file)
for a, p, b in zip(aln_len[10:150], parmikRecall[10:150], blastRecall[10:150]):
    print(f"{a:<10} {p:<15.2f} {b:<15.2f}")

