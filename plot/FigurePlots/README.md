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

### Existing Figure 3 Slurm jobs

| Purpose | Script (click for its repository location) |
| --- | --- |
| PARMIK counts | [cmp_parmik_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_parmik_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm) |
| BLAST counts | [cmp_blast_K11_PI90.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_blast_K11_PI90.slurm) |
| Bowtie2 counts | [cmp_bowtie_K11_PI90.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_bowtie_K11_PI90.slurm) |
| Minimap2 counts (edit input first) | [cmp_minamp2_K11_PI90.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/cmp_minamp2_K11_PI90.slurm) |
| PARMIK alignments | [aln_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm](../../experiments/parmik/IK_SSW_PF/SRR12432009/NC_045512/aln_IKT0_K11_PI90_M4_T1_SC1_P_2484_1111_2444_2288_2848.slurm) |
| BLAST alignments | [aln_ws11_pi90_mts50M.slurm](../../experiments/blast/SRR12432009/NC_045512/mts_50M/aln_ws11_pi90_mts50M.slurm) |
| Bowtie2 alignments | [SRR12432009_NC_045512.slurm](../../experiments/bwotie2/SRR12432009/NC_045512/SRR12432009_NC_045512.slurm) |

Update paths and scheduler settings for your machine, then submit the desired job:

```bash
sbatch <script_path>
```

Wait for alignment jobs to finish before submitting their comparison jobs;
plot once the counts are ready. `cmp_minamp2_K11_PI90.slurm` retains its original
filename typo and currently points `-b` at Bowtie2 output: replace that argument
with your minimap2 SAM file before using it.

List all available Slurm scripts from the repository root:

```bash
rg --files -g '*.slurm'
```

## 3. Generate Figure 3 and recall numbers

```bash
python3 plot/FigurePlots/recall_PARMIK_BLAST_fig3.py <counts_directory> --output <figure.pdf>
```

Use the experiment directory containing the count files as `<counts_directory>`.
Every tool with both TP and FN files is plotted; missing pairs are skipped.
Invalid data produces an error. The command saves `<figure.pdf>` and
`<figure>_recall.txt` beside it, creates output directories automatically,
and prints recall percentages. No graphical display is needed.
