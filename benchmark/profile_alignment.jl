using ProfileAlignments, Statistics

const PA = ProfileAlignments
const fixture = joinpath(@__DIR__, "..", "results", "benchmark", "proteins_80.fasta")
isfile(fixture) || error("Run benchmark/compare.jl once to create the benchmark fixture")
records = read_fasta(fixture)

function measured(f, label)
    GC.gc()
    result = nothing
    elapsed = @elapsed result = f()
    println(rpad(label, 30), round(elapsed; digits=3), " s")
    result
end

# Compile the relevant paths before reporting steady-state timings.
warmup = align_sequences(records[1:4]; retree=2, max_iterations=0, return_trace=true)
PA.refine_alignment(warmup.alignment; tree=warmup.tree, max_iterations=2,
                    return_trace=true)
align_sequences(records)

base = measured("progressive (retree=1)") do
    align_sequences(records; retree=1, max_iterations=0, return_trace=true)
end
retreed = measured("progressive + retree") do
    align_sequences(records; retree=2, max_iterations=0, return_trace=true)
end
refined = measured("refinement only") do
    PA.refine_alignment(retreed.alignment; tree=retreed.tree, max_iterations=2,
                        return_trace=true)
end
complete = measured("complete defaults") do
    align_sequences(records; return_trace=true)
end

println("representatives: ", length(base.representatives))
println("refinement accepted/cycles: ", refined.accepted, "/", refined.cycles)
println("complete accepted/cycles:   ", complete.accepted, "/", complete.cycles)
println("same alignment: ", getfield.(refined.alignment.records, :sequence) ==
                            getfield.(complete.alignment.records, :sequence))
