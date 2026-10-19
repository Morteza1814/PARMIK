# PARMIK project discussion and tool comparison

This README consolidates the shareable project context and recommendations from the PARMIK discussion. It covers the paper and available data, the earlier tool comparisons, Minimap2 and other suggested comparators, and revision points raised in prior reviews. Internal system/tool instructions and hidden content are not included.

Prepared 4 October 2026. Literature and software details below are a research summary, not results from a new benchmark. Options and defaults can change; pin versions and retain each tool's full help output before running comparisons.

## Project context from the discussion

PARMIK means “PArtial Read Matching using Inexpensive K-mers.” The intended task is to search many short nucleotide query sequences (around 150 nt) against a large collection of individual short metagenomic reads, retrieve qualifying local matches, and report the matched target read, alignment coordinates, and edit information under minimum alignment-length (`R`) and percent-identity (`PI`) criteria.

The current manuscript supplied in the project is:

`/u/rgq5aw/GIT/PARMIK/PARMIK__PArtial_Read_Matching_using_Inexpensive_K_mers__PLOS_ONE_.pdf`

The earlier conversation also identified the GitHub repository as <https://github.com/Morteza1814/PARMIK>. A prior submission introduced the project to Jake and Ram, who offered to add a Minimap2 comparison. The proposed collaboration was to retain the existing comparisons, add Minimap2 first, and discuss the interpretation and manuscript changes with Morteza. Morteza noted that he was not very familiar with Minimap2 and wanted to discuss it further after the group read more about it.

### Data and experiment context

The previously described query sources were SARS-CoV-2 (`NC_045512`) and SARS-CoV-1 (`NC_004718`) reference genomes. Queries were described as 150-nt sliding windows with a one-base step (approximately 29,754 and 29,602 windows, respectively). This is a set of overlapping windows, not a set of independent simulated reads; do not describe them as error-simulated unless a script confirms added noise.

The metagenomic accessions discussed were:

- `SRR12432009`
- `SRR14381418`
- `SRR14381421`
- `SRR14381422`

The project data location previously given was `/p/lava/parmik/partialMatcherProject/datasets/metagenomicDB/<SRR>/<SRR>_reads`. The processed target files were described as retaining reads of exactly 150 nt, removing ambiguous-base reads when configured, and removing duplicate sequences. The earlier description says that the processed metagenomic reads were supplied as the reference/database and the viral 150-nt fragments as queries. The exact preprocessing options and counts should be taken from the scripts and run logs for each experiment.

The earlier conversation identified experiment artifacts under `/p/lava/parmik/experiments` and PARMIK output under `/p/lava/parmik/partialMatcherProject/tools/parmik/IK_SSW_PF`. Paths on that cluster may not be available to collaborators outside the environment; share or transfer files as needed.

## What was compared in the paper/code

The previously examined experiments included PARMIK, BLASTN, Bowtie2, and BWA. In the described workflow, the recall analysis compared PARMIK/BLAST/Bowtie2 with a Smith–Waterman baseline; the best-match-per-query analysis included BWA, BLAST, Bowtie2, and PARMIK. The proposed Minimap2 addition was initially scoped to the best-match analysis, while broad Minimap2 secondary alignments could be investigated separately for recall. This distinction is a study-design proposal, not a claim that Minimap2 cannot report secondary hits.

The parameter values previously associated with several PARMIK runs included `K=11`, `M=4` (or `M=7` for some experiments), `PI=90%`, `R=30` in several summaries, an inexpensive-k-mer threshold that was either unrestricted or 5,000, and polishing where indicated. Treat these as reported examples, not as universal settings for every figure; recover the exact per-figure commands from scripts and logs.

### Existing command details that need careful checking

- A BLASTN script was found using word size 11, a large `max_target_seqs` value (50 million in one script), and output columns including query/target IDs, coordinates, alignment length, mismatches, strand, aligned sequences, identity, gaps, and E-value. Word size alone does not make BLAST exhaustive or ensure it retains every alignment satisfying PARMIK's `R`/`PI` rules.
- A Bowtie2 script was found using local mode, `-a`, and seed length 11 (`--local -a -L 11`). `-a` affects reporting limits, but Bowtie2's local-score threshold and heuristic search still matter.
- A BWA-MEM invocation using `-a` was present for a best-hit comparison. BWA was not part of the described recall-baseline comparison. This is an experimental-design distinction, not a statement that BWA is incapable of reporting multiple alignments.
- PARMIK's `-a` modes were described in the repository as: `0` index, `1` align, `2` compare best hits, `3` baseline, and `4` compare against baseline. Confirm these against the exact code revision used for the paper.
- In the inspected CLI, `-m` refers to a minimum exact-match length, while the paper's `M` denotes a k-mer-count threshold. Do not conflate similarly named variables. The command-line conversion between them should be described and tested against the version used.
- PARMIK's `-t` was described as the inexpensive-k-mer threshold (with zero meaning all), whereas `-w` controls threads in the inspected CLI. Archive the program help with each run to avoid ambiguity.

### Two high-priority benchmark audits

1. **Bowtie2 local score cutoff.** Bowtie2's documented default local minimum score is `20 + 8 ln(query_length)`, around 60 for a 150-nt query. With default FASTA scoring, a 30-column alignment with 27 matches and 3 substitutions scores 36 (`27×2 − 3×6`): it meets 90% identity but falls below that score cutoff. Thus `--local -a -L 11` alone does not establish a fair search for short `R=30`, `PI=90%` alignments. Recover the actual `--score-min` and scoring options from the run scripts/logs and test them on known fixtures.
2. **SAM edit-distance parsing.** The inspected SAM reader appears to use the SAM `NM` edit-distance tag as substitutions and then add CIGAR indels, potentially double-counting indels and miscounting matches. Audit the parser using synthetic SAM records with substitutions, insertions, deletions, soft clipping, and reverse-strand alignments before relying on published identity/edit-distance comparisons. This is a code-audit concern, not a conclusion that every historical result is wrong.

## Tool comparison: direct local alignment and mapping

“Can report coordinates” means the tool can produce an alignment trace/coordinates in an appropriate output mode. It does not mean the tool returns every possible qualifying interval, uses PARMIK's exact `R`/`PI` rule, or has been benchmarked on millions of separate short target records.

| Tool | Role and useful controls | Alignment detail / useful output | Main limitation for this study |
|---|---|---|---|
| **PARMIK** | Inexpensive-k-mer index/filter followed by alignment and polishing. Parameters include `K`, matched-k-mer threshold `M`, `PI`, minimum length `R`, IKT, scoring schemes, and second-chance/polishing behavior. | Per-read match output and alignment information are central to the method; validate coordinate and edit calculations with independent fixtures. | Novelty must be stated as the demonstrated combination and trade-offs, not as inventing k-mer indexes, frequent-k-mer filtering, local alignment, or inverted files. Ablate each distinctive step. |
| **Minimap2** | `-x sr` for short reads; `-k` k-mer length; `-w` minimizer window; `-f`, `-U`, `--q-occ-frac` affect seed-frequency filtering; `-n` minimum minimizers per chain; `-m` minimum chain score; `-s` minimum DP score; `--secondary`, `-N`, `-p` control secondary alignments; `-A`, `-B`, `-O`, `-E`, `-z` and `--end-bonus` affect scoring/extension; `-t` threads. | `-a` SAM or `-c` CIGAR-bearing PAF; `--eqx`, `--MD`, `--cs` preserve edit detail. | `sr` disables secondary alignments by default. Even with secondaries enabled, seed/chaining is heuristic and does not promise all `R`/`PI`-qualifying local hits. Several score/chain parameters are not equivalent to `R` or `PI`. Changing indexed seed parameters may require rebuilding the index. |
| **BWA-MEM / BWA-MEM2** | `mem`; `-k` seed length; `-r` reseeding; `-c` occurrence cap; `-w` band; `-d` Z-drop; `-A/-B/-O/-E/-L` scoring/clipping; `-T` minimum output score; `-a` report all found alignments for single-end/unpaired reads; `-t` threads. | SAM CIGAR/MD and NM tags; reconstruct edit positions from trace plus sequences. | “All found” is not exhaustive enumeration. Mapping is primarily a best-alignment task and scoring/output thresholds can omit short partial matches. Keep as a best-hit comparator, and separately define the recall objective. |
| **Bowtie2** | `--local`; `--score-min`; `--ma`, `--mp`, `--rdg`, `--rfg`; `-L`, `-N`, `-i`, `-D`, `-R`; `-k`/`-a`; `-p`; `--xeq`. | SAM positions, CIGAR and optional `=/X` operations; MD can help reconstruct edits. | It has no direct `R`/`PI` acceptance rule. Default local score can exclude matches satisfying short length/identity thresholds. `-a` does not remove heuristic search limits. |
| **BLASTN** | `-task`, `-word_size`, `-perc_identity`, `-evalue`, `-qcov_hsp_perc`, reward/penalty and gap scores, X-drop settings, masking options, strand, `-max_target_seqs`, `-max_hsps`, `-num_threads`. | Tabular output can include IDs, coordinates, strand, identity, gaps, aligned query/target strings, BTOP, and E-value. | E-value/score filters optimize statistical significance, not the exact retrieval of all `R`/`PI`-eligible read pairs. `blastn-short` is intended for much shorter queries than 150 nt. Candidate/output limits and masking require explicit audit. |
| **MMseqs2 nucleotide search** | Use the native nucleotide workflow (`--search-type 3`); relevant controls include `--min-seq-id`, `--seq-id-mode`, `--min-aln-len`, `-e`, coverage controls, `-k`, exact/spaced-kmer options, masking, `--max-seqs`, `--max-accept`, `--max-rejected`, gap scores, Z-drop, and threads. | Use traceback (`-a`) and `convertalis` fields for coordinates, aligned strings/CIGAR, identity, length, and IDs. | Not exhaustive by virtue of a shared identity threshold; prefilter/candidate caps and stopping parameters matter. Native nucleotide defaults differ from generic/protein search documentation; pin release and workflow. |
| **Lambda3** | Current release offers nucleotide indexing/search commands (`mkindexn`, `searchn`); paper describes fast/default/sensitive modes and filtering by identity, E-value, and bit score. | Supports BLAST-like and SAM/BAM-style reporting as described in its paper. | Strong modern nucleotide-search candidate, but verify exact release-specific flags, hit caps, threading, and output fields with `--full-help`; do not borrow Lambda2 options blindly. Published results on assembled reference databases do not establish performance on millions of 150-nt records. |
| **LAST** | `lastal`; `-T` local alignment; `-e` score threshold; `-m/-l/-L/-k/-W` seed controls; scoring/extension parameters `-r/-q/-a/-b/-z/-x/-y/-d`; strand/threads; pruning and hit limits such as `-K/-C/-N/-n`. | MAF can retain coordinates, strand, and gapped sequences for common post-filtering. | `-k` is not simply k-mer length. Repeat masking, pruning, statistical cutoffs, and hit caps can remove otherwise eligible matches. Useful optional comparator, but older than the main modern-tool shortlist. |
| **STELLAR / Stellar3** | Error threshold `-e/--epsilon`; minimum length `-l/--minLength`; seed k; abundance/repeat filters; verification and X-drop; result limits. | GFF can include coordinates, identity, strand, query interval, CIGAR, and mutation positions. | Particularly relevant related work because it targets local matches with length/error bounds. Its output is maximal ε-matches under its definition, not every overlapping qualifying interval. Original STELLAR defaults included a per query-target-pair cap; verify current Stellar3 flags/defaults. |
| **SigAlign** | Rust API offers `Local` and `SemiGlobal`; local constructor sets mismatch/gap penalties, minimum length, and maximum penalty per length. Reference builder accepts multiple FASTA targets. | Per-target alignments can include endpoints and match/substitution/insertion/deletion operations. | Maximum penalty per length is not percent identity. Basic local mode applies overlap/output rules; limited mode may stop early. Needs a benchmark driver and memory/thread characterization for this workload. |
| **Strobealign** | Short-read mapper; `-r` sets mean read length profile, and seeding/profile controls include `-k`, `-s`, `-l`, `-u`, `-m`; noisy profile available. | SAM plus `--eqx`; `-N` bounds secondary output; mapping-only PAF does not preserve base-level trace. | A good modern short-read mapping comparison, but ordinary mapper outputs and secondary caps do not imply complete local-hit recovery. Check meanings against the pinned version's `--help`. |
| **Minibwa** | Minimap2-derived short-read mapper; `-x sr`; seeds/chaining (`-k/-c/-m/-w/-W/-g`); candidate retention `-N`, emitted secondaries `--outn`, score ratios; SAM/PAF, `-b cs|ds|MD`, `--eqx`; threading/index-memory controls. | Traceback-capable SAM/PAF modes. | Newer relevant candidate; score/candidate caps still differ from exhaustive `R`/`PI` search. Avoid chain-only mode for edit-position comparisons; use native interface rather than assuming legacy BWA flags behave identically. |
| **PuffAligner / Pufferfish** | Pufferfish index options include reference, output, k, sparsity, and build threads; aligner accepts index and reads. | CIGAR output is documented on a `cigar-strings` branch. | Originally oriented to transcript quantification/end-to-end mapping; verify local-match behavior, single-end flags, hit reporting, and output in the exact build. `align -p` is binary PAM output, not a thread setting. |
| **X-Mapper** | Java mapper with configurable maximum penalty, substitution/indel scores, near-best span, and thread/cache controls. | SAM output; verify needed trace tags or reconstruct from sequences. | Paper describes end-to-end alignment; no verified arbitrary internal-local `R` mode. Maximum penalty per length is not identical to PI. `--max-num-matches` may reject a query on excess matches rather than simply truncate output. |
| **DREAM-Stellar** | Parallel, space-efficient exact local alignment; error/length controls and indexing/binning parameters. | Can report local matches with coordinates, CIGAR, and mutation details. | Very relevant objective-matched comparator. Understand its probabilistic filter and repeat-bin behavior; test whether the output policy returns all qualifying target pairs. Not the same as a generic best-hit mapper. |

### Other direct tools considered

| Tool | Assessment |
|---|---|
| **BLAT** | Historical local nucleotide search; options include tile/step size, minimum score/identity, repeat-seed cap, and fine search. Useful background or optional comparator, but PSL blocks alone may not fully specify every edit. |
| **USEARCH local** | Local nucleotide-search mode with identity/E-value/coverage and acceptance/hit limits. Potential older comparator; verify edition, licensing, and stopping behavior. |
| **VSEARCH** | Its documented feature set does not provide the same arbitrary local alignment task; global identity/coverage is not a substitute for partial matches. |

## Tool comparison: k-mer retrieval, pseudoalignment, and classification

These tools can retrieve candidate sequences efficiently, but usually do not by themselves provide the same per-read local alignment trace as PARMIK. A strong reviewer-facing experiment is to use one as a candidate generator, then run a common aligner on retrieved reads and compare recall, coordinates, and cost.

| Tool | Main controls/capabilities | What it does not establish by itself | Fit for a follow-up |
|---|---|---|---|
| **Themisto** | Colored de Bruijn graph; sequence/sample colors; threshold-union style query matching. | A matched color or k-mer threshold is not an alignment interval, PI, edit distance, or traceback for a particular read. | Candidate generation followed by common local alignment. Preserve per-read identity through colors/input design. |
| **COBS** | Compact k-mer index; configurable term size and false-positive rate; threshold search and result limits. Supports multifasta handling for separate target records. | Approximate k-mer membership is not a local alignment or edit trace; verify false-positive/negative effects and target-record conventions. | Candidate filter followed by exact read-level alignment. Good reviewer response if indexed with one document per read. |
| **Fulgor** | Colored compacted de Bruijn graph; can return color candidates and optionally query k-mer positions. | Query positions do not automatically provide target coordinates or a local gapped alignment. | Candidate stage followed by alignment; document sequence-to-color mapping and any lost read multiplicity. |
| **Raptor** | Interleaved Bloom filter index; seed/bin configuration and threshold-based candidate lookup. | Candidate/bin hits are not a local alignment proof or per-base edit trace. | Candidate generation plus common alignment; report false candidate burden and sensitivity. |
| **Mantis** | Compact exact k-mer membership/count queries using a counting quotient filter. | K-mer membership alone has no local alignment coordinates or PI. | Background or candidate prefilter, if the specific workflow returns read identities. |
| **ganon / ganon2** | Taxonomic classification/database workflows; custom database construction can operate on sequences/reads with settings affecting sources, minimizers, and classification. | Taxonomic assignment is not a per-read local alignment with coordinates and edit positions. | Related-work context or upstream classifier; not a direct replacement without a read-level retrieval plus alignment stage. |
| **Kraken2** | k-mer-based taxonomic classification; database and confidence/threshold options. | Usually assigns taxonomy rather than returning every qualifying target read and local traceback. | Related work and possible candidate source only if a read-level mapping back to original targets is designed. |
| **MetaGraph** | De Bruijn graph plus annotations; build k, annotate headers/coordinates, query signatures/counts/coordinates, and align reads. Can preserve coordinate-to-header mapping with sidecar metadata. | Graph alignment is not automatically equivalent to independently aligning every original read; graph paths and annotation coordinates must be mapped back carefully. | Potential candidate/alignment system. Build the index so each original read is recoverable, preserve coordinate mapping, and validate output IDs, strand, and CIGAR. |
| **LOCALE** | Learned local-alignment embeddings retrieve top candidate sequences/accessions; search quality depends on index/model and candidate settings. | Candidate ranking is not a guarantee of all threshold-qualified read pairs, nor an edit trace. | Emerging literature candidate for screening; follow with common exact local alignment if testing. Confirm availability and exact short-read behavior. |

### Key retrieval distinctions

- Indexing query or target k-mers, returning candidate read IDs, computing a local alignment, and emitting all qualifying hits are four separate steps. Measure and describe each.
- A threshold-union or shared-k-mer count is a candidate criterion. It is not equivalent to local alignment length `R` or percent identity `PI`.
- Graph compaction may merge sequence paths. A graph method can still be useful if annotations preserve the original read IDs and coordinates, but those mappings must be verified.
- If a method cannot emit edit positions directly, preserve the original query/target sequences and compute a common traceback on candidates.

## How to define and compare matches

Before comparing tools, specify the match rule in one place and apply it to every output. A clear starting definition for an alignment with `M` matches, `X` substitutions, `I` inserted query bases, and `D` deleted query bases is:

`alignment_columns = M + X + I + D`

`PI = 100 × M / alignment_columns`

`edit_count = X + I + D`

State whether terminal soft clips are excluded, how ambiguous bases are treated, and whether a gap of length `g` contributes `g` errors or one gap event. Report aligned query and target spans separately. Do not assume a tool's “identity,” “edit distance,” or “alignment length” field uses this same convention.

Keep three evaluation tasks distinct:

1. **Pair recall:** for each query-target read pair with at least one qualifying alignment, did the tool recover that pair?
2. **Best match per query:** did the tool choose a preferred target under a defined ranking rule?
3. **Coordinate/trace agreement:** do start/end coordinates, strand, CIGAR/edit positions, and alignment score agree under a defined truth set?

A best-scoring Smith–Waterman alignment is not necessarily the longest subalignment that meets `R` and `PI`. A baseline that uses the same alignment implementation/scoring choices as PARMIK is useful for a controlled comparison but is not independent biological truth. Label it accurately and add independent synthetic/known-truth tests.

## Recommendations from the reviewer comments

### Highest priority: correctness and fair evaluation

1. **Audit every comparator's eligibility filters.** For BLASTN, record task, word size, E-value, masking, scoring/gap parameters, target/HSP caps, and output format. For Bowtie2, record local score minimum, scoring, seed settings, search effort, and reporting mode. Test fixtures where the exact `R`/`PI` outcome is known.
2. **Validate the SAM parser and every post-filter.** Create small alignments with substitutions, insertions, deletions, clipping, multiple matches, and reverse-strand hits. Compare computed lengths, PI, edit counts, and coordinates against manually checkable truth.
3. **Build independent ground truth.** Use controlled synthetic query-target pairs with known substitutions/indels and known nonmatches. Include composition-matched shuffled or unrelated controls, but do not assume unrelated species are guaranteed negatives. Report pair-level precision and recall, not just agreement with a PARMIK-like baseline.
4. **Define the biological meaning and limitations of a hit.** A short high-identity local match may occur by chance or in a conserved/repetitive region; it does not by itself prove organism presence or support a clinical decision. Analyze match uniqueness, complexity, coverage, and genomic context, and explicitly limit claims.
5. **Separate performance measures.** Report index build time, load time, query/alignment time, postprocessing time, peak RAM, index disk size, and output size. Report both queries/second and qualifying pairs/second, because tools emitting different hit counts make a single speed number misleading.

### High priority: tool coverage and novelty

6. **Retain BWA, Bowtie2, and BLASTN, then add Minimap2.** Do not assume Minimap2 must outperform a short-read mapper because it was first known for long reads; evaluate its short-read preset and task-specific seeding. Decide whether it belongs in best-hit, pair-recall, or both analyses based on the actual secondary-hit output and validation.
7. **Prioritize additional direct comparators:** Lambda3 and native nucleotide MMseqs2 are strong candidates; DREAM-Stellar and STELLAR are especially relevant to the bounded local-match objective; SigAlign is useful for direct local traceback. LAST is an optional additional local-search baseline. Pilot on a manageable subset before expanding the matrix.
8. **Address the index-plus-realignment reviewer suggestion directly.** Test a k-mer retrieval system (for example Themisto, COBS, Fulgor, Raptor, or MetaGraph) as candidate generation, then apply the same alignment and filtering step to candidates. Compare candidate recall, final pair recall, coordinate trace, and build/query/RAM costs.
9. **Narrow novelty claims.** Inverted indexes, k-mer filtering, local alignment, and alignment coordinates are established ideas. State precisely what PARMIK contributes as a combination, how its IKT and candidate filters differ operationally, and which measured trade-offs are better on the tested workload. Cite the read-alignment landscape and closely related local-search/filtering work.
10. **Ablate the claimed contributions.** Test one versus five scoring schemes; with/without IKT; multiple IKT values; with/without polishing/second chance; and alternate trimming choices. Keep postfilters and data fixed so each change has an interpretable effect.

### High priority: clarify method and dataset

11. Define “partial match” as a local alignment and clarify how that differs from a full-query mapping. Give a concise high-level statement early: many 150-nt queries are searched against many individual 150-nt target reads; the goal is to retrieve read-level local alignments that satisfy specified `R` and `PI` thresholds.
12. Define “inexpensive k-mer” versus “cheap k-mer,” state how frequent k-mers are identified, and explain the trade-off. A frequent k-mer can still be informative; filtering it may lower index/query cost but can also lose sensitivity. Avoid implying that common means uninformative in every dataset.
13. Explain that IKT is an absolute occurrence threshold in the current implementation and report the range tested. The claim that 5,000–10,000 worked across the tested datasets is empirical for those datasets, not a universal guarantee. Evaluate normalized/adaptive alternatives or state them as future work with a portability limitation.
14. Define `K`, the matched-k-mer threshold (`M`), `R`, `PI`, and IKT once with general notation. Distinguish `M` in the manuscript from any CLI variable such as minimum exact-match length. Replace `%` when it could be confused with percent identity; use `mod` for modulo.
15. Describe the index concretely (e.g. mapping k-mers to target/read identifiers and positions), whether it is in memory or on disk, whether construction time is included, and what exact reads/records and query counts were used in each result.
16. Explain scoring and gap conventions: match reward, substitution penalty, gap-open and gap-extension definitions, and whether gap-open includes the first gap base. State how a local alignment is selected, its score/length cutoff, and what the “second chance”/polishing step changes.
17. Explain the five scoring models and evaluate whether taking the best result across models increases false positives. Report the contribution of each model and define tie-breaking.
18. Clarify whether `R` means alignment-column length or exact-match length. The paper and code/comments appear to use “exact” and “alignment” length in places that may not mean the same thing; make the acceptance rule unambiguous.

### Figures, results, and writing

19. Put the 150-nt read/query scale and dataset sizes in the abstract/introduction and repeat exact query/read counts, filtering rules, and accession identifiers in Methods or a compact data table.
20. Include construction time in a separate metric if search time excludes it. State whether timing includes index creation/loading, output writing, SAM parsing, and common postfiltering.
21. Report full benchmark results and modern comparisons, including index size and peak memory as reviewers requested. Include exact versions, commands, parameters, hardware, threads, and reproducible scripts in supplement/repository.
22. Redesign figures for legibility: larger labels, consistent subfigure lettering/colors/fonts, distinguish zero from missing data, use consistent plot types, clarify zoomed panels, and reduce unrelated panels per figure. Consider replacing sparse bars with line, point-density, or heatmap plots where appropriate.
23. Explain that BLAST E-value estimates expected chance alignments for a score and database size; it is not directly the same acceptance rule as a minimum local `R`/`PI` threshold. Compare the actual sequences/alignments where tools disagree and report which differences satisfy independent truth.
24. Clarify simulation/benchmark inputs: how many queries, their lengths, whether mutations or sequencing errors were introduced, which strands were searched, and the number of target reads after each preprocessing step. A sliding-window fragment set should be described as such.
25. Discuss limitations and biological interpretation. More reported local matches do not automatically mean more true organism detections. Explain whether additional matches are longer, unique, plausible, repetitive, low-complexity, or potentially chance matches.
26. Clean manuscript wording and references: introduce PARMIK before using the acronym; explain terms such as “directly searchable” and “robust metagenomic sampling”; complete the cut-off sentence; use primary Smith–Waterman and inverted-index references rather than Wikipedia; fix citation years/links and typos; normalize PI capitalization and subfigure punctuation.
27. Address algorithm presentation: give pseudocode/numbered steps, relate the approach to seed-chain-extend, explain where simple k-mer counting substitutes for chaining and where it does not, and state complexity in terms of query count, candidate count, and alignment lengths—not as an all-to-all Smith–Waterman operation if filtering avoids most pairs.
28. Explain the reviewer questions about why increasing `M` can improve speed in some settings: requiring more matching k-mers can reduce the number of candidates passed to expensive alignment, even though computing/counting more evidence can itself cost time. Show measured trade-offs; do not state a monotonic speed rule without data.
29. Revisit the reported short-match/polishing example where a shorter pre-polish match became longer after polishing. Explain the algorithm's selection behavior and whether the longer exact segment was on the same strand/trace; use the code and figure source to avoid implying that the method simply overlooked a contained exact match.
30. Consider moving implementation-heavy parameter sweeps and complete command tables to supplementary material, while making the main figures tell one question each.

## Suggested comparison plan

| Stage | Work | Decision/output |
|---|---|---|
| 1. Reproduce | Recover per-figure scripts, exact commit, data preprocessing, tool versions, and logs. | A manifest linking each figure/table to its command, inputs, outputs, and plotting script. |
| 2. Validate | Audit parser and thresholds with small known-truth fixtures; independently calculate R/PI/edit traces. | Corrected, tested common evaluation layer. |
| 3. Pilot | Run existing tools plus Minimap2, Lambda3, MMseqs2 nucleotide, and one objective-matched aligner (DREAM-Stellar or STELLAR) on a subset. | Feasibility and evidence of which output modes/parameters retrieve relevant pairs. |
| 4. Retrieval pipeline | Test at least one sequence/k-mer index with common realignment (e.g. COBS or Themisto plus traceback aligner). | Reviewer-responsive comparison of candidate retrieval versus read-level alignment. |
| 5. Expand | Run selected parameters on the full accession/query matrix with resource accounting and independent controls. | Final tables/figures with pair recall, best-hit quality, precision controls, build/query costs, RAM, and index size. |
| 6. Revise | Rewrite claims/methods/results/discussion and archive commands/data manifests. | Reproducible manuscript and supplement. |

## Prioritization

**Do first:** retrieve exact per-figure commands and versions; inspect Bowtie2 score thresholds; test SAM parsing; make a common R/PI/coordinate definition; establish independent ground truth and negative controls.

**Then pilot:** Minimap2 short-read mode; Lambda3; native nucleotide MMseqs2; DREAM-Stellar or STELLAR; a candidate-index-plus-realignment workflow (COBS or Themisto are practical examples).

**Then revise:** strengthen related work and novelty framing; give full input sizes and reproducible methods; add ablations and cost metrics; redesign figures; expand discussion of chance matches and biological limits.

## Sources for tool details and key methodological statements

- Minimap2 manual: <https://lh3.github.io/minimap2/minimap2.html>
- BWA manual: <https://bio-bwa.sourceforge.net/bwa.shtml>
- BWA-MEM2: <https://github.com/bwa-mem2/bwa-mem2>
- Bowtie2 manual: <https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>
- BLASTN options and output: <https://www.ncbi.nlm.nih.gov/books/NBK279684/>
- MMseqs2 parameters: <https://github.com/soedinglab/MMseqs2/blob/master/src/commons/Parameters.cpp>; nucleotide workflow: <https://github.com/soedinglab/MMseqs2/blob/master/src/workflow/Search.cpp>
- Lambda3 repository and paper: <https://github.com/seqan/lambda>; <https://pmc.ncbi.nlm.nih.gov/articles/PMC10955267/>
- STELLAR: <https://github.com/seqan/seqan/tree/main/apps/stellar>; <https://pmc.ncbi.nlm.nih.gov/articles/PMC3283304/>
- DREAM-Stellar: <https://link.springer.com/article/10.1186/s12859-026-06389-0>; <https://github.com/seqan/dream-stellar>
- SigAlign API and paper: <https://docs.rs/sigalign/latest/sigalign/>; <https://academic.oup.com/nar/article/52/15/8717/7714453>
- LAST manual: <https://gitlab.com/mcfrith/last/-/raw/main/doc/lastal.rst>
- Strobealign: <https://github.com/ksahlin/strobealign>
- Minibwa: <https://github.com/lh3/minibwa>; <https://arxiv.org/abs/2606.15357>
- Pufferfish/PuffAligner: <https://github.com/COMBINE-lab/pufferfish>; <https://doi.org/10.1093/bioinformatics/btab408>
- X-Mapper: <https://github.com/mathjeff/Mapper>; <https://link.springer.com/article/10.1186/s13059-024-03473-7>
- Themisto: <https://github.com/algbio/themisto>; <https://pmc.ncbi.nlm.nih.gov/articles/PMC10311346/>
- COBS: <https://cobs.readthedocs.io/en/latest/tutorial.html>
- Fulgor: <https://github.com/jermp/fulgor>
- Raptor: <https://docs.seqan.de/raptor/main/tutorial_search.html>
- MetaGraph quick start and sequence search: <https://metagraph.ethz.ch/static/docs/quick_start.html>; <https://metagraph.ethz.ch/static/docs/sequence_search.html>
- ganon documentation and paper: <https://pirovc.github.io/ganon/custom_databases/>; <https://pmc.ncbi.nlm.nih.gov/articles/PMC12267982/>
- Kraken2 manual: <https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown>
- Mantis: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10964368/>
- Alser et al., read-alignment survey: <https://pubmed.ncbi.nlm.nih.gov/34446078/>
- Smith–Waterman (1981): <https://pubmed.ncbi.nlm.nih.gov/7265238/>
- Parasail SIMD alignment library: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4748600/>
- Inverted files for text search engines (Zobel & Moffat): <https://www.inf.ed.ac.uk/teaching/courses/ad/lectures16/ZobelMoffat.pdf>

## Conversation scope

### Conversation history captured here

The project discussion began with locating the PARMIK paper, finding the repository and identifying where scripts/results might be stored. The user supplied a local PDF path and later asked about finding the paper online and sharing a newer manuscript. The discussion then moved to restoring the project after a pause, locating the metagenomic read datasets, and adding Minimap2 with help from Jake and Ram.

The email-drafting discussion covered: thanking the students; describing the initial Minimap2 comparison; retaining BWA while adding Minimap2; explaining that Minimap2 should not be assumed better or worse for short reads; potentially limiting the first comparison to best match per query; keeping recall as a separate experiment if broad secondary hits can be configured and validated; noting that MMseqs2 might be considered later; describing the viral sliding-window queries and metagenomic read database; and proposing a Friday meeting at 3 p.m. Mountain Time. The user also wanted to offer to share datasets and mention that Professor Skadron was copied for reference. Exact final email wording should be taken from the email actually sent, if one exists.

There was also a practical compression exchange about the processed read files. `gzip -9 -k FILE` writes `FILE.gz` in the current directory and keeps the original. `tar` can combine multiple files into one archive; a compressed tar archive is conventionally `.tar.gz`, while `gzip` alone compresses one input stream/file. Faster compression generally uses a lower gzip level, trading file size for speed. These are command-line notes from the conversation, not part of the scientific evaluation.

The reviewer feedback supplied by the user included two prior journal reviews, with concerns about terminology and definitions, acronym order, local versus partial matching, seed/k-mer novelty, frequent-k-mer filtering, IKT portability, algorithm detail, score cutoffs, input sizes, comparison coverage, baseline independence, precision and negative controls, output-coordinate interpretation, runtime/index construction, RAM/index size, figure clarity, and biological interpretation of additional matches. The recommendations above consolidate those points into an action plan rather than copying the long review text verbatim.

This file is therefore a readable consolidated record of the project-relevant discussion, not a line-by-line transcript. Private/internal system instructions, tool orchestration, and non-shareable content are omitted.
