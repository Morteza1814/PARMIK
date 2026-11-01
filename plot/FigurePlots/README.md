# Figure commands

Run from the repository root. Replace `<...>` with your own values.
Run only the steps you need; existing count files can go straight to plotting.

## 1. Build

```bash
make -C src
```

For plotting, activate a Python environment with NumPy, pandas and Matplotlib.
On this server, the environment is already installed:

```bash
source .venv-plot-miniforge/bin/activate
```

For a new environment on UVA CS (once):

```bash
module load miniforge/26.7.2
python3 -m venv .venv-plot-miniforge
source .venv-plot-miniforge/bin/activate
python3 -m pip install numpy pandas matplotlib
```

## 2. Generate data for Figure 3

Generate TP/FN counts from existing tool alignments and a matching
Smith–Waterman baseline. Run once per tool you want to evaluate:

```bash
./src/parmik -a 4 -l <tool> \
  -s 48 -c 150 -m 14 -t 0 -k 11 -d 90 -x \
  -i <read_count> -j <query_count> \
  -r <reads_file> -q <queries_file> -f <index_directory> \
  -b <tool_alignment_file> -z <baseline_directory> \
  -p <penalty_file> -o <results_directory>
```

Use `parmik`, `blast`, `bowtie`, or `minimap2` for `<tool>`.
The original Figure 3 uses 17,000,000 reads and up to 30,000 queries;
use the same processed inputs and penalties as the saved baseline.
This step can require substantial memory and overwrites existing comparison counts.

Counts are written in the `experiment dir` printed by PARMIK, beneath
`<results_directory>`, as `cmp_Baseline_<tool>_AlnSz_tp.txt` and
`cmp_Baseline_<tool>_AlnSz_fn.txt`. For `parmik`, the program reads
`parmikAln.txt` from that experiment directory instead of the supplied `-b` path.
This step evaluates saved alignments; it does not rerun the aligners.

## 3. Generate Figure 3 and recall numbers

```bash
python3 plot/FigurePlots/recall_PARMIK_BLAST_fig3.py <counts_directory> --output <figure.pdf>
```

Use the experiment directory containing the count files as `<counts_directory>`.
Every tool with both TP and FN files is plotted; missing pairs are skipped.
Invalid data produces an error. The command saves `<figure.pdf>` and
`<figure>_recall.txt` beside it, creates output directories automatically,
and prints recall percentages. No graphical display is needed.
