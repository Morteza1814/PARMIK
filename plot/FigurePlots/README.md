# Figure commands

Run from the repository root. Replace `<...>` with your own values.
Run only the steps you need; current baseline recall count files can go straight
to plotting. Rebuild and regenerate comparisons when changing the recall definition.

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
Comparison queries must have unique, consecutive numeric IDs starting at 0,
ending at the loaded query count minus 1. The comparison rejects unsupported IDs
instead of silently skipping queries.

For parallel comparisons in mode `-a 4`, append `--comparison-threads 8` to
the command above. This option defaults to 1 and must be a positive integer.
It controls only comparison workers; it does not change the experiment
directory or the PARMIK alignment input path. Keep `-w` at the value used for
the saved experiment, since `-w` still participates in its directory name.

Append `--comparison-summary-only` to omit per-query details from
`cmp_Baseline_<tool>.txt` while retaining its overall summary and all count
reports. Detailed reporting remains enabled by default. Baseline data is loaded
one 1,000-query chunk at a time, then queries in that chunk are distributed
among workers. Counts are merged and detailed reports retain query order.
Loading each chunk before processing its queries also fixes the previous case
where an `N` in query 0, 1000, etc. prevented loading the baseline for the rest
of that chunk. Results for such inputs can therefore differ from older runs.

When benchmarking, compare 1, 2, 4, and 8 workers on the same inputs with the
same reporting option. Input parsing and disk I/O can limit scaling; use the
worker count that improves elapsed time on your machine. Each run overwrites
the comparison reports in the same experiment directory, so run these trials
sequentially and save reports separately if you want to compare them.

### Existing Figure 3 Slurm jobs

| Purpose | Script (click for its repository location) |
| --- | --- |
| PARMIK counts | [cmp_parmik_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_parmik_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm) |
| BLAST counts | [cmp_blast_K11_PI90.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_blast_K11_PI90.slurm) |
| Bowtie2 counts | [cmp_bowtie_K11_PI90.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_bowtie_K11_PI90.slurm) |
| Minimap2 counts | [cmp_minamp2_K11_PI90.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_minamp2_K11_PI90.slurm) |
| PARMIK alignments | [aln_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/aln_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm) |
| BLAST alignments | [aln_ws11_pi90_mts50M.slurm](../../experiments/blast/SRR12432009/NC_045512/mts_50M/aln_ws11_pi90_mts50M.slurm) |
| Bowtie2 alignments | [SRR12432009_NC_045512.slurm](../../experiments/bwotie2/SRR12432009/NC_045512/SRR12432009_NC_045512.slurm) |

Update paths and scheduler settings for your machine, then submit the desired job:

```bash
sbatch <script_path>
```

Wait for alignment jobs to finish before submitting their comparison jobs;
plot once the counts are ready. `cmp_minamp2_K11_PI90.slurm` retains its original
filename typo; its input now points to `minimap2.sam`. Counts from the previous
version used Bowtie2 output and must be regenerated before plotting minimap2.

List all available Slurm scripts from the repository root:

```bash
rg --files -g '*.slurm'
```

## 3. Generate Figure 3 and recall numbers

```bash
python3 plot/FigurePlots/recall_PARMIK_BLAST_fig3.py <counts_directory> --output <figure.pdf>
```

Use the experiment directory containing the count files as `<counts_directory>`.
Rebuild PARMIK and rerun the baseline comparisons for **all tools** before plotting.
The plot now requires `cmp_Baseline_<tool>_BaselineRecallByThreshold.txt` and
rejects older comparison reports for tools missing this file. It also rejects
different threshold grids or different eligible-baseline counts between tools,
which catches incompatible count sets. All runs must still use the same actual
baseline and query set; equal count totals alone do not prove input identity.

`BaselineRecallByThreshold.txt` has three columns: `R`, `TP_baseline(R)`,
`FN_baseline(R)`. These are direct counts, not histogram bins, and must not be
cumulatively summed. For each R, the eligible set contains accepted baseline
query/read pairs whose baseline alignment length B is at least R. A pair is TP
when its best accepted tool alignment has length L at least R; otherwise it is
FN. Thus `TP_baseline(R) + FN_baseline(R)` is the same baseline denominator for
every tool. An absent or rejected tool alignment uses L=0, and a missed pair
contributes FN when `L < R <= B`. Queries with no tool candidates contribute all
their eligible baseline pairs. The plotted percentage is
`100 * TP_baseline(R) / (TP_baseline(R) + FN_baseline(R))`; a zero denominator is
shown as 0%.

Two separate reports preserve the other information:

- `cmp_Baseline_<tool>_AdditionalHitsByThreshold.txt` has two columns: R and the
  number of accepted tool hits of length at least R outside the eligible baseline
  set. This includes pairs absent from the accepted baseline and pairs whose
  baseline alignment is shorter than R. These hits do not increase baseline recall.
- `cmp_Baseline_<tool>_RecallByThreshold.txt` retains the previous three-column
  counts using all accepted tool hits as TP. Its `TP / (TP + FN)` ratio is a
  separate metric, rather than recall against a fixed baseline-positive set.

The older `AlnSz_tp.txt` and `AlnSz_fn.txt` histograms remain available for
alignment-length summaries, but cannot reconstruct the threshold relationship.
Invalid data produces an error. The command saves `<figure.pdf>` and
`<figure>_recall.txt` beside it, creates output directories automatically,
and prints recall percentages. No graphical display is needed.
