using AASV, DataFrames, CSV, Statistics
root=joinpath(@__DIR__,"..","results","benchmark")
times=vcat(CSV.read(joinpath(root,"julia_times.csv"),DataFrame),CSV.read(joinpath(root,"r_times.csv"),DataFrame))
summary=combine(groupby(times[times.phase.=="warm",:],[:language,:operation,:n]),:seconds=>median=>:seconds)
sort!(summary,[:operation,:n,:language])
CSV.write(joinpath(root,"summary.csv"),summary)
show(stdout,summary;allrows=true)
println()
julia=sort(classify_asv(joinpath(root,"julia_output")),[:ASV_id,:tax_lev])
r=sort(CSV.read(joinpath(root,"r_verdicts.csv"),DataFrame),[:ASV_id,:tax_lev])
@assert nrow(julia)==nrow(r)==30
@assert julia.ASV_id==r.ASV_id && julia.tax_lev==r.tax_lev
println("Maximum FI difference on these synthetic inputs: ",maximum(abs.(julia.FI.-r.FI)))
println("Classes agreeing: ",count(julia.Class.==r.Class),"/",nrow(r))
println("Julia ",VERSION,"; OS ",Sys.KERNEL,"; hardware threads ",Sys.CPU_THREADS)
