# Divergent filter: local measurements, 15 October 2026

Windows, Julia 1.12.2. Comparison of the distance computation before and after
optimisation, on synthetic sequences of 220 positions. One call after warm-up,
with exact equality of the results checked before timing.

| Sequences | Distinct | Before (s) | After (s) | Speed-up | Allocations before | Allocations after |
|---|---:|---:|---:|---:|---:|---:|
| 400 | 400 | 0.1255123 | 0.0200584 | 6.26× | 452,541,550 bytes | 58,148 bytes |
| 1,000 | 50 | 0.6708134 | 0.0003065 | 2,189× | 2,831,735,566 bytes | 36,467 bytes |

In the second scenario each distinct sequence appears 20 times: every record
therefore has a duplicate at distance zero and no pair scans are needed.
This is a particularly favourable case, not a prediction for real data.
Allocations are the cumulative total, not the peak resident memory.
Excluded from the timings: start-up, package loading, FASTA reading/writing,
checkpoints and reports. Reproduce with:
`julia --project=. benchmark/divergent_filter.jl`.
