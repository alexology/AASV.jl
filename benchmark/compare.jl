using AASV, DataFrames, CSV, Statistics, Random

root=joinpath(@__DIR__,"..","results","benchmark")
mkpath(root)
rng=MersenneTwister(812)
alphabet=collect("ACDEFGHIKLMNPQRSTVWY")
base=rand(rng,alphabet,150)
base[1]='M'; base[30]='A'; base[60]='A'
codons=Dict(aa=>first(sort([c for (c,a) in AASV.genetic_table("5") if a==aa])) for aa in alphabet)
to_dna(aa)=join(codons[a] for a in aa)

fixtures=Dict{Int,Vector{Sequence}}()
for n in (20,80)
    records=Sequence[]
    for i in 1:n
        seq=copy(base)
        positions=randperm(rng,150)[1:4]
        for p in positions
            p in (1,30,60) && continue
            seq[p]=rand(rng,alphabet)
        end
        push!(records,Sequence("Bench_species__REF$i",String(seq)))
    end
    fixtures[n]=records
    write_fasta(records,joinpath(root,"proteins_$n.fasta"))
end

references=joinpath(root,"references")
records=fixtures[20]
for (level,name) in (("species","Bench_species"),("genus","Bench"),("family","BenchFamily"))
    write_fasta(records,joinpath(references,"alignments","BenchFamily",level,name*".fasta"))
end
write_fasta([Sequence(r.id,to_dna(r.sequence)) for r in records],
    joinpath(references,"raw_sequences_nt","BenchFamily","species","Bench_species.fasta"))
asvs=DataFrame(family=fill("BenchFamily",10),genus=fill("Bench",10),species=fill("Bench species",10),
    ASV=[begin
        aa=copy(base)
        isodd(i) && (aa[30]='D')
        i%3==0 && (aa[60]='K')
        to_dna(aa[13:142])
    end for i in 1:10],sample1=fill(10,10))
CSV.write(joinpath(root,"asvs.csv"),asvs)

timings=DataFrame(language=String[],operation=String[],n=Int[],phase=String[],seconds=Float64[])
pipeline_only="--pipeline-only" in ARGS
if pipeline_only && isfile(joinpath(root,"julia_times.csv"))
    previous=CSV.read(joinpath(root,"julia_times.csv"),DataFrame)
    append!(timings,previous[previous.operation.!="functional_pipeline",:])
end
function measure(operation,n,fn)
    time=@elapsed fn()
    push!(timings,("Julia",operation,n,"first",time))
    for repeat in 1:3
        GC.gc()
        time=@elapsed fn()
        push!(timings,("Julia",operation,n,"warm",time))
    end
    println(operation," ",n,": warm median ",median(timings.seconds[(timings.operation.==operation) .& (timings.n.==n) .& (timings.phase.=="warm")]))
    flush(stdout)
end
if !pipeline_only
    for n in (20,80)
        measure("multiple_alignment",n,()->align_sequences(fixtures[n]))
    end
end
out=joinpath(root,"julia_output")
measure("functional_pipeline",10,()->begin
    asv_functional_structure!(asvs;output_dir=out,alignment_dir=references)
    classify_asv(out)
end)
CSV.write(joinpath(root,"julia_times.csv"),timings)
println("Benchmark fixtures and timings: ",abspath(root))
