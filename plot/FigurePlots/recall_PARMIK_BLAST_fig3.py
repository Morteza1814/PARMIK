import pandas as pd
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def plotFig(x, y2, y3, y4, output_file, minimap2_recall=None, show=True):
    second_ax_name = "PARMIK"
    third_ax_name = "BLAST"
    fourth_ax_name = "BOWTIE"

    # Create a figure and axes for both plots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    bar_width = 0.2

    # Plot using dots
    ax1.plot(x, y2, '-', label=second_ax_name, color='navy', markersize=10)
    ax1.plot(x, y3, '-', label=third_ax_name, color='firebrick', markersize=10)
    ax1.plot(x, y4, '-', label=fourth_ax_name, color='darkgreen', markersize=10)

    if minimap2_recall is not None:
        ax1.plot(x, minimap2_recall, '-', label='MINIMAP2', color='darkorange', markersize=10)

    ax1.legend(fontsize=18, loc='lower right')
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.set_ylabel('Recall Rate (%)', fontsize=20)
    ax1.text(-0.1, 1.1, '(A)', transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

    zoom = (x >= 21) & (x < 35)
    x_zoom = x[zoom]

    # Plot the zoomed-in portion as a bar chart in the bottom subplot
    ax2.bar(x_zoom, y2[zoom], label=second_ax_name, color='navy', width=bar_width)
    ax2.bar([i + bar_width for i in x_zoom], y3[zoom], label=third_ax_name, color='firebrick', width=bar_width)
    ax2.bar([i + 2 * bar_width for i in x_zoom], y4[zoom], label=fourth_ax_name, color='darkgreen', width=bar_width)

    if minimap2_recall is not None:
        ax2.bar(x_zoom + 3 * bar_width, minimap2_recall[zoom],
                label='MINIMAP2', color='darkorange', width=bar_width)

    ax2.legend(fontsize=18, loc='lower right')
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax2.set_ylim(0, 100)
    ax2.set_ylabel('Recall Rate (%)', fontsize=20)
    ax2.text(-0.1, 1.1, '(B)', transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top', ha='right', color='blue')

    plt.xlabel('Minimum Alignment Size (R)', fontsize=20)
    plt.tight_layout()
    plt.savefig(output_file)
    if show:
        plt.show()
    plt.close(fig)
    return fig

def cumulative_sum_from_each_element(arr):
    return np.cumsum(arr[::-1])[::-1]

def read_first_and_last_column(file_path):
    df = pd.read_csv(file_path, sep=r'\s+', header=None)
    if df.empty or df.shape[1] != 3:
        raise ValueError(f"{file_path}: expected length, baseline count, tool count")
    values = df.to_numpy(dtype=float)
    if (not np.isfinite(values).all() or (values < 0).any()
            or (values != np.floor(values)).any()):
        raise ValueError(f"{file_path}: histogram values must be nonnegative integers")
    lengths = values[:, 0].astype(np.int64)
    if (lengths < 1).any() or (np.diff(lengths) <= 0).any():
        raise ValueError(f"{file_path}: lengths must be positive and strictly increasing")
    return lengths, values[:, -1].astype(np.int64)


def calculate_recall_rate(tpArray, fnArray):
    if tpArray.shape != fnArray.shape:
        raise ValueError("Both arrays must be of the same size.")
    recall_rate_array = np.zeros_like(tpArray, dtype=float)
    for i in range(len(tpArray)):
        denominator = tpArray[i] + fnArray[i]
        if denominator > 0:
            recall_rate_array[i] = (tpArray[i] / denominator) * 100
        else:
            recall_rate_array[i] = 0
    return recall_rate_array

def getRecallOutOfTpFn(tpFilePath, fnFilePath):
    aln_len, tp = read_first_and_last_column(tpFilePath)
    fn_len, fn = read_first_and_last_column(fnFilePath)
    if not np.array_equal(aln_len, fn_len):
        raise ValueError("TP and FN alignment lengths must match")
    tp_cum = cumulative_sum_from_each_element(tp)
    fn_cum = cumulative_sum_from_each_element(fn)
    recall = calculate_recall_rate(tp_cum, fn_cum)
    np.set_printoptions(suppress=True)
    return aln_len, recall

def trim_arrays_within_range(array1, lower_limit, upper_limit):
    mask = (array1 >= lower_limit) & (array1 <= upper_limit)
    trimmed_array1 = array1[mask]
    return trimmed_array1

def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot Figure 3 from saved TP/FN histograms")
    parser.add_argument('experiment_path', type=Path)
    parser.add_argument('--include-minimap2', action='store_true',
                        help='require minimap2 counts; otherwise include them when present')
    parser.add_argument('--output', type=Path, help='PDF destination (default: input directory)')
    parser.add_argument('--no-show', action='store_true', help='save without opening a display')
    args = parser.parse_args(argv)
    tools = ['parmik', 'blast', 'bowtie']
    minimap2_files = [args.experiment_path / f'cmp_Baseline_minimap2_AlnSz_{kind}.txt'
                      for kind in ('tp', 'fn')]
    if args.include_minimap2 or any(path.exists() for path in minimap2_files):
        if not all(path.is_file() for path in minimap2_files):
            parser.error('minimap2 requires both cmp_Baseline_minimap2_AlnSz_tp.txt '
                         'and cmp_Baseline_minimap2_AlnSz_fn.txt')
        tools.append('minimap2')

    recalls = []
    lengths = None
    for tool in tools:
        paths = [args.experiment_path / f'cmp_Baseline_{tool}_AlnSz_{kind}.txt'
                 for kind in ('tp', 'fn')]
        tool_lengths, recall = getRecallOutOfTpFn(*paths)
        if lengths is not None and not np.array_equal(lengths, tool_lengths):
            raise ValueError(f"{tool}: alignment lengths differ from PARMIK")
        lengths = tool_lengths
        recalls.append(recall)

    mask = (lengths >= 11) & (lengths <= 150)
    if not mask.any():
        raise ValueError('No histogram lengths in the Figure 3 range (11 to 150)')
    x = lengths[mask]
    plotted = [recall[mask] for recall in recalls]
    filename = ('PARMIKvsBLASTvsBOWTIEvsMINIMAP2_recall.pdf' if len(tools) == 4
                else 'PARMIKvsBLASTvsBOWTIE_recall.pdf')
    output = args.output or args.experiment_path / filename
    plotFig(x, *plotted[:3], output,
            minimap2_recall=plotted[3] if len(tools) == 4 else None,
            show=not args.no_show)
    print('R          ' + ' '.join(f'{tool.upper():<15}' for tool in tools))
    for row in zip(x, *plotted):
        print(f"{row[0]:<10} " + ' '.join(f'{value:<15.2f}' for value in row[1:]))


if __name__ == '__main__':
    main()
