# Local measurements, 14 October 2026

Windows, Julia 1.12.2 and R 4.5.0, on the same computer. DECIPHER with
`processors=1`; the two benchmarks were never run at the same time.
Median of the three calls following the first one, in seconds:

| Operation | R | Julia |
|---|---:|---:|
| Alignment of 20 proteins × 150 aa | 0.280 | 0.089 |
| Alignment of 80 proteins × 150 aa | 1.000 | 0.527 |
| Functional analysis + classification of 10 ASVs at three levels | 11.600 | 0.284 |

The pipeline uses 20 pre-aligned references and 390 nt ASVs. It includes
FASTA and Excel input/output; it excludes BOLD and VSEARCH. On these inputs,
FI and classes match for all 30 ASV/level combinations.

First pipeline call, excluding process start-up and package import:
Julia 13.838 s, R 12.660 s. Julia's initial compilation therefore weighs
heavily compared with the subsequent calls.

On these inputs, the Julia engine is about 3.1 times faster with 20 proteins
and 1.9 times faster with 80; the functional and classification pipeline is
also faster in steady state. These data do not establish which complete
workflow is faster on a real dataset: the alignment engines differ and the
data here are synthetic.

Profiling of the previous version attributed almost all the time of the
80-sequence case to refinement: 32.66 s out of about 33.9 s. Scoring now walks
the ASCII strings directly, and refinement evaluates only the pairs crossed by
each bipartition, caching the others. The 80-sequence median dropped from
33.913 s to 0.527 s. UPGMA selection also uses a min-queue of distances, with
cost O(N² log N) instead of O(N³).

The original R code fails on a query without leading gaps; the pipeline is
measured on fragments shorter than the references to avoid that case without
modifying the original source. See README.md.

Reproducible scripts: `compare.jl`, `compare.R`, `summarize.jl`.
Raw measurements: `../results/benchmark/julia_times.csv`, `r_times.csv` and
`summary.csv` in the same results directory.
