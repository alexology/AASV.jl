# Local R / Julia comparison

To compare only the distance computation of the divergent filter with the
previous implementation:

```powershell
julia --project=. benchmark/divergent_filter.jl
```

The benchmark first checks that the distances are equal, then times one call
after warm-up, reporting time and allocated bytes. It includes 400 sequences
that are all different and 1,000 sequences with 50 distinct sequences, each 220
positions long. The data are synthetic: the measurements exclude FASTA
reading/writing, checkpoints and reports, and do not represent the time of the
whole pipeline.
The local measurements are in `benchmark/DIVERGENT_FILTER_RESULTS.md`.

Run sequentially, avoiding other CPU load:

```powershell
julia --project=. benchmark/compare.jl
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' benchmark/compare.R results/benchmark reference/aasv
```

The first script creates identical fixtures for both languages under
`results/benchmark` and times Julia; the second reads the fixtures and times
the original R code at the commit given in the project README.

Scenarios: multiple alignment of 20 and 80 synthetic proteins of 150 residues,
and functional analysis + classification of 10 ASVs (390 nt) against 20
references at species, genus and family level. The pipeline profiles are
already aligned and identical for both versions; the pipeline time includes
FASTA/Excel I/O. It does not include BOLD downloads or VSEARCH.

The ASVs are fragments covering positions 13–142 of the reference proteins,
which have 150 residues. In the original R source, a query with zero leading
gaps causes an error in the subset `audit_results[-seq_len(0)]`; the pipeline
comparison uses fragments to avoid that case without modifying the R code.
The results therefore do not represent all possible queries.
`--pipeline-only` repeats the pipeline while keeping the alignment
measurements, both in the Julia script and as the third argument of the
R script.

Each scenario contains the first call and three subsequent calls; report the
median of the three subsequent calls. The first Julia call includes any
function compilation, but **not** process start-up or package loading, which
happen before the measurement.

The command above starts Julia with a single thread, so DECIPHER uses
`processors=1` and ProfileAlignments uses one thread. To also measure Julia's
parallel mode, run `julia -t auto --project=. benchmark/compare.jl`
separately. The comparison uses the normal settings of the two engines, which
use different algorithms, weights and parameters. The timings measure these
implementations on these inputs, not the intrinsic performance of R versus
Julia. They are not benchmarks of biological accuracy.

Raw results: `results/benchmark/julia_times.csv` and `r_times.csv`.
The summary of the latest local run is in `benchmark/RESULTS.md`.
