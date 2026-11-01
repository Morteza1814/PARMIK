import pandas as pd
import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
STYLES = {
    'parmik': ('PARMIK', 'navy'),
    'blast': ('BLAST', 'firebrick'),
    'bowtie': ('BOWTIE2', 'darkgreen'),
    'minimap2': ('MINIMAP2', 'darkorange'),
}


def tool_label(tool):
    return STYLES.get(tool, (tool.upper(), None))[0]


def plotFig(x, recalls, output_file, show=False):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    bar_width = min(0.2, 0.8 / len(recalls))
    zoom = (x >= 21) & (x < 35)
    for index, (tool, recall) in enumerate(recalls.items()):
        label, color = STYLES.get(tool, (tool.upper(), None))
        line, = ax1.plot(x, recall, '-', label=label, color=color)
        ax2.bar(x[zoom] + index * bar_width, recall[zoom],
                label=label, color=line.get_color(), width=bar_width)
    for ax, panel in ((ax1, '(A)'), (ax2, '(B)')):
        ax.legend(fontsize=18, loc='lower right')
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.set_ylabel('Recall Rate (%)', fontsize=20)
        ax.text(-0.1, 1.1, panel, transform=ax.transAxes, fontsize=24,
                fontweight='bold', va='top', ha='right', color='blue')
    ax2.set_ylim(0, 100)
    ax2.set_xlabel('Minimum Alignment Size (R)', fontsize=20)
    fig.tight_layout()
    fig.savefig(output_file)
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
    parser = argparse.ArgumentParser(description="Save Figure 3 and recall numbers for every available tool")
    parser.add_argument('experiment_path', type=Path, help='directory containing TP/FN count files')
    parser.add_argument('--output', type=Path,
                        default=REPO / 'figure-results/figure3/Figure3.pdf')
    parser.add_argument('--show', action='store_true', help='also open a display')
    parser.add_argument('--no-show', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--include-minimap2', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if not args.experiment_path.is_dir():
        parser.error(f'Data directory does not exist: {args.experiment_path}')
    prefix = 'cmp_Baseline_'
    available = set()
    for kind in ('tp', 'fn'):
        suffix = f'_AlnSz_{kind}.txt'
        available.update(p.name[len(prefix):-len(suffix)]
                         for p in args.experiment_path.glob(f'{prefix}*{suffix}'))
    tools = [tool for tool in STYLES if tool in available]
    tools += sorted(available - set(STYLES))
    recalls = {}
    lengths = None
    for tool in tools:
        paths = [args.experiment_path / f'{prefix}{tool}_AlnSz_{kind}.txt'
                 for kind in ('tp', 'fn')]
        if not all(path.is_file() for path in paths):
            print(f'Skipping {tool}: both TP and FN files are needed.', file=sys.stderr)
            continue
        tool_lengths, recall = getRecallOutOfTpFn(*paths)
        if lengths is not None and not np.array_equal(lengths, tool_lengths):
            raise ValueError(f'{tool}: alignment lengths differ between tools')
        lengths = tool_lengths
        recalls[tool] = recall
    if not recalls:
        parser.error('No complete TP/FN histogram pairs found')
    mask = (lengths >= 11) & (lengths <= 150)
    if not mask.any():
        parser.error('No histogram lengths in the Figure 3 range (11 to 150)')
    plotted = {tool: recall[mask] for tool, recall in recalls.items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plotFig(lengths[mask], plotted, args.output, show=args.show and not args.no_show)
    lines = ['R          ' + ' '.join(f'{tool_label(tool):<15}' for tool in plotted)]
    for row in zip(lengths[mask], *plotted.values()):
        lines.append(f"{row[0]:<10} " + ' '.join(f'{value:<15.2f}' for value in row[1:]))
    table = '\n'.join(lines) + '\n'
    numbers = args.output.with_name(args.output.stem + '_recall.txt')
    numbers.write_text(table)
    print(table, end='')
    print(f'Saved plot: {args.output}\nSaved numbers: {numbers}', file=sys.stderr)


if __name__ == '__main__':
    main()
