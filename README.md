# AASV.jl

Julia implementation of the [aasv](https://github.com/alexology/aasv) pipeline
for the functional analysis of COI ASVs. The scientific formulas are compatible
with R version 0.1.0, commit `9960347948da3142714deb656fcf38124b4e38d8`;
the architecture uses Julia types, dispatch and conventions.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

Requires Julia 1.10 or later. Does not require R, Biostrings or DECIPHER.

```powershell
julia --project=. setup.jl
julia -t auto --project=. test/runtests.jl
julia -t auto --project=. examples/offline_pipeline.jl
```

The local example creates FASTA files, codon tables, Excel reports and
`results/demo/analysis/classification.csv`.

To install the version published on GitHub (ProfileAlignments.jl is not
registered, so add it first; `setup.jl` does the same for a local checkout):

```julia
using Pkg
Pkg.develop(url="https://github.com/alexology/ProfileAlignments.jl")
Pkg.develop(url="https://github.com/alexology/AASV.jl")
using AASV, BioSequences

record = Sequence("COI", dna"ATGCCAGCA")
protein = best_translation(record.sequence)
```

## Pipeline with real data

```julia
using AASV, DataFrames, CSV

taxonomy = CSV.read("taxonomy.csv", DataFrame)
# Required columns: family, genus, species.
align_species_seq!(taxonomy; alignment_dir="references",
                   min_length=500, max_length=1000, genetic_code="5")

for family in unique(skipmissing(taxonomy.family))
    problems = check_alignment(family; alignment_dir="references", tax_lev=SpeciesRank)
    improve_alignment!(problems; log_dir="logs")
    dir = joinpath("references", "alignments", family, "species")
    paths = filter(p -> endswith(p, ".fasta"), readdir(dir; join=true))
    filter_divergent_seq!(paths; log_dir="logs", max_distance=0.10)
end

align_genus_seq!(; alignment_dir="references")
align_family_seq!(; alignment_dir="references")
# Records whose insertions open gap columns in most of the merged profile
# are removed and listed in references/logs/gap_causing_references.

asvs = CSV.read("asvs.csv", DataFrame)
# Required columns: family, genus, species, ASV (DNA sequence).
# Additional numeric columns: abundances in the samples.
revised = asv_functional_structure!(asvs;
    alignment_dir="references", output_dir="analysis", genetic_code="5")
verdicts = classify_asv("analysis")
CSV.write("classification.csv", verdicts)
```

`align_genus_seq!` and `align_family_seq!` remove, after each merge, the
records whose insertions open columns that are gaps in most references
(`remove_gap_causing=true`). R aasv merges without this step and relies on
`check_alignment` and `improve_alignment` at each level; `false` reproduces
its merges. A few such records otherwise make every query an indel at that
level: on 9,768 Scarabaeidae references, 44 records (0.45%) made all 34
tested ASVs `INVALID_ORF` at family level; with the removal the results
matched those obtained with the R references.

The package precompiles a small offline analysis, so the first call of
`asv_functional_structure!` in a session does not wait for compilation.

Julia options are passed after `;`. Missing statistical values use `missing`;
an absent result or option uses `nothing`. CSV, Excel and `DataFrame` are
input/output adapters, while the core uses concrete structures. Sequences are
vectors of `Sequence(id, sequence)` with unique names. The `alignment_dir` and
`output_dir` arguments accept Windows, Linux and macOS paths. Use
`overwrite=true` to rebuild references that already exist.

## Julia API

`Sequence` is the same type exported by `ProfileAlignments`: alignments no
longer copy each record between two domain representations. The typed
translation returns a single result type:

```julia
translation = find_translation(dna"ATGCCAGCA")
if !isnothing(translation)
    translation.amino_acids
    translation.frame
end
```

`assess_structure` and `assess_asv` return `StructureAssessment` and
`ASVAssessment` respectively; amino acid differences are `MutationObservation`
values. These APIs compose the computation in memory without creating
intermediate tables. `check_asv_structure` and `calculate_asv_distance` keep
the tabular output for the existing reports.

Taxonomic options accept `SpeciesRank`, `GenusRank` and `FamilyRank`.
A reusable FI policy is represented by an immutable value:

```julia
policy = ClassificationConfig(quality_threshold=0.995,
                              include_hydro=false)
verdicts = classify_asv("analysis", policy)
```

Functions that modify files or directories end with `!`:
`align_species_seq!`, `align_genus_seq!`, `align_family_seq!`, `re_align!`,
`improve_alignment!`, `filter_divergent_seq!` and
`asv_functional_structure!`. The previous names remain as adapters for
compatibility with existing code.

`filter_divergent_seq!` reports progress over files and, every 5 seconds,
over the sequence comparisons in the current file. To show it even with
`print_messages=false`, pass `show_progress=true`; the frequency is set with
`progress_interval` (in seconds).
The distance computation works directly on the bytes of the ASCII sequences,
without copying rows for each pair, and precomputes the gap-free ends. Identical
sequences share their comparisons while keeping all the original identifiers.
Distances and removal rules remain exact, including cases with no overlap or
with gaps only. The gain grows with the number of duplicates; when all
sequences differ, the maximum number of comparisons remains quadratic.

To process all families without writing loops in the script:

```julia
species_bad = check_alignments!(; alignment_dir=reference_dir,
    tax_lev=SpeciesRank)
isempty(species_bad) || improve_alignment!(species_bad)

# After the fix, only the modified families are checked again.
species_still_bad = check_alignments!(; alignment_dir=reference_dir,
    tax_lev=SpeciesRank)

filter_divergent_references!(; alignment_dir=reference_dir,
    tax_lev=SpeciesRank, max_distance=0.10)
```

Both functions discover the families and skip directories without FASTA files.
`families=["FamilyA", "FamilyB"]` restricts the selection. `resume=true`
(default) reuses the results saved in `alignment_dir/checkpoints` only if the
file names and SHA-256 contents, the parameters and the algorithm version all
match. New downloads, added or removed files and realignments also invalidate
the results when the content changes. `resume=false` forces a new run.

On a cluster, `align_species_seq!`, `align_genus_seq!`, `align_family_seq!`,
`check_alignments!` and `filter_divergent_references!` accept
`shard_index` and `shard_count`. Families are sorted and partitioned
deterministically, without overlaps. For example, an 8-element job array
uses `shard_index=task_id, shard_count=8`. Jobs may share `alignment_dir`,
because each shard writes distinct families and checkpoints; downloads or
unsharded runs must not be started in parallel on the same files.

`asv_functional_structure!` instead uses dynamic tasks within the node and
accepts `threads=:auto` or an integer limit. ASV computations run in parallel
in batches bounded by the number of threads, while the aggregated reports are
assembled in the same order as the input.
Each reference profile and its nucleotide references are read once per run
and released after the last ASV of the family; reference codons are
translated only at novel residues.

`asv_functional_structure!` also reads reference sets built by earlier aasv
versions. Exact repeats of a record (same ID and sequence) are ignored; an ID
repeated with a different sequence is an error. Profile records without a
nucleotide reference in the taxon cannot be checked codon by codon and are
excluded, together with the columns left with gaps only. This removes, for
example, the `Canthonidia` species that an earlier R version merged into the
`Canthon` genus profile. Both cases are reported once per file with a warning. The other functions
that read references (`check_alignment`, `improve_alignment!`, `re_align!`,
`filter_divergent_seq!` and the genus/family merges) also ignore exact
repeats, and a family merge keeps a species found in two genus profiles once.
`check_reference_records(; alignment_dir)` lists these records in a table
before an analysis:

```julia
issues = check_reference_records(; alignment_dir="references")
```

Genus and family profiles can be rebuilt with `overwrite=true`. Species
alignments with altered IDs need `align_species_seq!(overwrite=true)`:
realigning them from `raw_sequences` would restore records removed by earlier
curation.
`check_alignments!` also returns the problematic paths of reused families;
it additionally supports `GenusRank`, `FamilyRank` and `type="raw_sequences"`.

The filter records the file fingerprints **after** modification and stores the
summary of removed sequences in the checkpoints. Per-family reports are
written to `alignment_dir/logs/divergent_references/<level>/<family>.tsv`
(customisable with `log_dir`, or disabled with `log_dir=nothing`).
Resumed runs do not duplicate report rows. Outputs are prepared before the
FASTA files are replaced: an interrupted replacement can be completed on the
next run without repeating the filter. If the interruption happens during the
computation, the family in progress is recomputed; completed families are
reused. Staging files remain available in case of interruption. The filter
works on the current files: a less restrictive threshold does not restore
previously removed sequences; to do that, regenerate the references.
`show_progress=false` disables the progress messages.

By default, `align_species_seq!` reports the families processed or reused,
progress through BOLD pages, filtering, and the species being aligned with its
number of sequences. Pass `show_progress=false` to disable these messages;
`progress_interval` controls updates during filtering. Alignment messages are
emitted before each species.
Independent species of the same family are aligned in batches with tasks from
the separate `ProfileAlignments` package; each batch is limited to the number
of available threads to bound memory use and avoid nested parallelism.
Alignment algorithms, presets and numerical structures live in the
`ProfileAlignments` package; AASV handles downloads, taxonomy and files.

`align_genus_seq!` and `align_family_seq!` report the families and taxa
processed or reused, the number of sequences in the profiles to merge, the
completed merges with the percentage for the current taxon, and the final time.
Updates are emitted before and after each merge; `show_progress=false`
disables them.

## Ported functions

| Function | Operation |
|---|---|
| `aa_similarity_matrix` | 4096 codon pairs with Hamming similarity |
| `align_species_seq!` | BOLD COI-5P download, filters, translation, alignments |
| `align_genus_seq!`, `align_family_seq!` | progressive profile merging |
| `check_alignment` | paths of alignments with internal gaps |
| `check_alignments!` | check of all families with persistent results |
| `filter_divergent_references!` | per-family filter with resume and persistent reports |
| `re_align!` | realignment, optionally from the original FASTA files |
| `improve_alignment!` | iterative removal of problematic sequences, with logs |
| `filter_divergent_seq!` | nearest-neighbour distance filter |
| `asv_functional_structure!` | ASV analysis at the three taxonomic levels |
| `classify_asv` | FI index, class, evidence and intermediate metrics |

Also available: `read_fasta`, `write_fasta`, `best_translation`, `aa_table`,
`align_sequences`, `align_profiles`, `calculate_asv_distance`,
`get_hydro_metrics` and `check_asv_structure`. `?function_name` shows the
documentation in Julia.

## Compatibility and differences

- The FI formula, the weights (0.35/0.30/0.35), the 0.90 threshold, the quality
  filter, the normalised entropy and the conservation discount
  `1-exp(-n_ref/20)` follow the R code. Grantham distances use the original
  table, not an approximation. Hydrophobicity uses the Kyte–Doolittle scale.
- Alignment uses the independent
  [ProfileAlignments.jl](../ProfileAlignments.jl/README.md) package, installed
  as a sibling project and written entirely in Julia:
  BLOSUM62, affine gaps (`gap_open=-10`, `gap_extend=-1`), k-mer distances,
  guide tree, two progressive passes and iterative refinement.
  Terminal gaps are free; ASV queries are inserted into the profile with
  `mode=:query`. `anchors=:fft` enables the approximate FFT acceleration,
  while the default uses full dynamic programming.
  The previous linear-penalty engine has been replaced.
  It does not reproduce DECIPHER or MAFFT: coordinates, gaps, mutations and FI
  may differ. Biological validation on real benchmarks is still required.
- `Sequence`, `best_translation` and the protein metrics accept BioSequences
  types directly (`LongDNA`, `LongAA` and the `dna"..."`, `aa"..."` literals).
  The pipeline converts these inputs once into the internal ASCII encoding,
  keeping compatibility with the files and with the string API.
- DECIPHER-specific arguments (for example `processors`, `terminalGap`) do not
  carry over directly; use the Julia parameters listed above. The previous
  arguments `match`, `mismatch` and `gap` are still accepted: `gap` sets a
  linear penalty.
  The filter distance is the fraction of differing residues in the
  intersection of the covered regions, excluding terminal gaps.
- Built-in genetic codes: the 20 codes of the Biostrings table used in the R
  comparison (1–6, 9–16, 21–26), without reinterpreting start codons. For
  other codes, pass a complete codon→amino acid dictionary.
  IUPAC bases that translate unambiguously produce the corresponding residue;
  otherwise `X`. Translation searches only the three forward frames, ignores
  the incomplete terminal codon and does not search the reverse complement.
- Identical profiles are deduplicated for the structural analysis, as in R.
  When restoring an entirely problematic alignment without an external
  reference, the first sequence is kept deterministically, whereas R picks
  one at random.
- BOLD is read through the paginated API, even for small families. Retries
  keep the backoff state; permanent HTTP 4xx errors stop the download
  immediately. A `.complete` marker prevents a directory created during a
  failed download from being treated as complete. Downloaded pages are saved
  to `alignment_dir/bold_downloads/<family>` only once fully received. With
  `overwrite=false`, a rerun reuses the local pages and requests only the
  missing pages of the same query; `overwrite=true` starts a new query without
  mixing pages from different downloads. The cache is kept after completion.
  If BOLD no longer serves the saved query, restart with `overwrite=true`.
  For local data, pass `records_provider=family -> records` to
  `align_species_seq`.
- The analysis also records comparisons without mutations or whose
  hydrophobicity cannot be computed, so no ASV is lost from the reports. When
  the hydrophobicity metric is missing, its FI contribution is zero, as in R;
  the report contains `UNAVAILABLE`. A stop codon or an incompatible gap
  produces the sentinel `aa_pos=-1` and an FI of zero.
- Reports keep the main names and columns of the R version.
  `classify_asv` reads Excel or TSV. Re-running the pipeline replaces the
  previous reports in the `aa_structure_results` directory.

The FI index measures compatibility with the evolution of functional COI;
it is neither a probability nor a diagnosis of NUMT/mitochondrial origin.

## Read quality with VSEARCH

```julia
asv_functional_structure(asvs;
    alignment_dir="references", output_dir="analysis",
    vsearch_path=raw"C:\tools\vsearch.exe",
    trimmed_folder="trimmed", ee=40)
```

Samples are selected from the numeric columns with abundance >0, by looking
for the sample name in the FASTQ file name. Uncompressed FASTQ files
(`.fastq`, `.fq`) with Phred+33 qualities are supported. Error probabilities
are `10^(-Q/10)`; the `ee` filter uses their sum, and codon quality is
`1-sum(errors)` clamped to [0,1]. VSEARCH errors stop the analysis. The
probability computation is tested locally; integration with a real VSEARCH
executable requires installing it.

## Local checks

The tests cover FASTA, translation, gaps, alignment cleaning, distances, Phred
quality, FI, simulated BOLD retries and the complete offline pipeline.
The optional comparison with R checks every numeric column and the class
across eight scenarios, as well as translation, hydrophobicity and genetic
codes:

```powershell
julia --project=. test/parity.jl 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
```

This check requires R with `dplyr`, `readxl`, `Biostrings`, `Peptides` and the
original sources under `reference/aasv` (or an alternative path as the second
argument). Normal use and the Julia tests do not require R.
The included CI runs the tests on Windows and Linux, with Julia 1.10 and the
stable release. The alignment package also has its own suite, which can be
run on its own:

```powershell
julia --project=../ProfileAlignments.jl -e 'using Pkg; Pkg.test()'
```

## Provenance

This port derives from [aasv](https://github.com/alexology/aasv), by Alex
Laini, Samuele Voyron and Laura Gruppuso, distributed under GPL >=2. This port
is also distributed under GPL >=2. The Grantham table comes from the
[grantham](https://github.com/patterninstitute/grantham) package, MIT licence;
see `THIRD_PARTY_NOTICES.md`. The copies of the original repositories in
`reference/` are local comparison material and are excluded from Git.
