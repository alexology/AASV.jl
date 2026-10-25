function alignment_matrix(records::AbstractVector{Sequence})
    isempty(records) && return Matrix{Char}(undef,0,0)
    n = length(first(records).sequence)
    all(r->length(r.sequence)==n,records) || throw(ArgumentError("Aligned sequences must have equal lengths"))
    mat = Matrix{Char}(undef,length(records),n)
    for (i,r) in enumerate(records)
        mat[i,:] = collect(r.sequence)
    end
    mat
end

function inner_gaps(mat::AbstractMatrix)
    [begin
        inds = findall(!=('-'), row)
        length(inds)>1 && any(==('-'), row[first(inds):last(inds)])
    end for row in eachrow(mat)]
end

function causes_internal_gaps(mat::AbstractMatrix)
    spanning = falses(size(mat))
    for i in axes(mat,1)
        inds = findall(!=('-'),mat[i,:])
        isempty(inds) || (spanning[i,first(inds):last(inds)] .= true)
    end
    flagged = [begin
        n = sum(spanning[:,j])
        n > 1 && sum(spanning[:,j] .& (mat[:,j] .== '-')) / n > 0.5
    end for j in axes(mat,2)]
    [any((mat[i,:] .!= '-') .& flagged) for i in axes(mat,1)]
end

function trim_columns(records)
    mat = alignment_matrix(records)
    keep = [any(!=('-'),col) for col in eachcol(mat)]
    [Sequence(r.id, String(mat[i,keep])) for (i,r) in enumerate(records)]
end

engine_records(records)=records
aasv_records(alignment)=alignment.records

function alignment_model(;model=nothing,alphabet=:protein,match=nothing,mismatch=nothing,
                         gap=nothing,gap_open=nothing,gap_extend=nothing,terminal_gap=0.0,mode=:overlap)
    !isnothing(model) && return model
    if !isnothing(gap)
        (!isnothing(gap_open) || !isnothing(gap_extend)) && throw(ArgumentError("Choose gap or affine gap_open/gap_extend"))
        gap_open=gap_extend=gap
    end
    ProfileAlignments.ScoringModel(;alphabet,match,mismatch,
        gap_open=isnothing(gap_open) ? -10.0 : gap_open,
        gap_extend=isnothing(gap_extend) ? -1.0 : gap_extend,terminal_gap,mode)
end

"""Align existing profiles through the standalone Julia ProfileAlignments engine.
Defaults: BLOSUM62, affine gap_open=-10/gap_extend=-1, free terminal gaps.
`mode=:query` fits the first profile into the second. FFT anchors are optional.
Legacy `gap` selects a linear penalty; `match`/`mismatch` select scalar scores.
"""
function align_profiles(a::AbstractVector{Sequence},b::AbstractVector{Sequence};
                        model=nothing,alphabet=:protein,match=nothing,mismatch=nothing,
                        gap=nothing,gap_open=nothing,gap_extend=nothing,terminal_gap=0.0,
                        mode=:overlap,kwargs...)
    scoring=alignment_model(;model,alphabet,match,mismatch,gap,gap_open,gap_extend,terminal_gap,mode)
    result=ProfileAlignments.merge_profiles(engine_records(a),engine_records(b);model=scoring,kwargs...)
    aasv_records(result)
end

"""Align sequences using guide-tree progressive alignment and iterative refinement.
Uses the independent Julia ProfileAlignments package; preserves IDs/order/duplicates.
"""
function align_sequences(records::AbstractVector{Sequence};model=nothing,alphabet=:protein,
                         match=nothing,mismatch=nothing,gap=nothing,gap_open=nothing,
                         gap_extend=nothing,terminal_gap=0.0,mode=:overlap,kwargs...)
    scoring=alignment_model(;model,alphabet,match,mismatch,gap,gap_open,gap_extend,terminal_gap,mode)
    aasv_records(ProfileAlignments.align_sequences(engine_records(records);model=scoring,kwargs...))
end

"""Return paths of FASTA alignments containing internal gaps."""
function check_alignment(taxon; alignment_dir=pwd(), type="alignments", tax_lev=SpeciesRank)
    tax_lev=rank_name(parse_rank(tax_lev))
    type in ("alignments","raw_sequences") || throw(ArgumentError("Invalid type"))
    files = fasta_files(joinpath(alignment_dir,type,safe_name(taxon),tax_lev))
    isempty(files) && throw(ArgumentError("No FASTA files in specified directory"))
    filter(files) do path
        records = read_reference_records(path)
        length(records)>1 && any(inner_gaps(alignment_matrix(records)))
    end
end

"""Realign FASTA alignments in place, without their gaps. With `raw_sequences=true`,
species alignments are rebuilt from `raw_sequences`, as in R aasv: this resets
the curation, restoring records removed by `improve_alignment!` or
`filter_divergent_seq!`. To rebuild species references, use
`align_species_seq!(overwrite=true)`.
"""
function re_align!(paths; raw_sequences=false, verbose=true, kwargs...)
    for path in paths
        if raw_sequences
            parts = splitpath(normpath(path))
            "species" in parts || throw(ArgumentError("raw_sequences=true requires species-level paths"))
            source = joinpath((p=="alignments" ? "raw_sequences" : p for p in parts)...)
        else
            source = path
        end
        write_fasta(align_sequences(read_reference_records(source);kwargs...),path)
        verbose && @info "Realigned" path
    end
    paths
end
re_align(args...;kwargs...)=re_align!(args...;kwargs...)

function append_log(log_dir, filename, rows)
    isnothing(log_dir) && return
    isempty(rows) && return
    mkpath(log_dir)
    path = joinpath(log_dir,filename)
    append = isfile(path) && filesize(path)>0
    CSV.write(path,DataFrame(rows); delim='\t', append, header=!append)
end

function align_blocks(records; kwargs...)
    groups = Dict{String,Vector{Sequence}}()
    for r in records
        push!(get!(groups,first(split(r.id,"__")),Sequence[]),r)
    end
    blocks = [align_sequences(groups[k];kwargs...) for k in sort(collect(keys(groups)))]
    sort!(blocks;by=length,rev=true,alg=MergeSort)
    reduce((a,b)->align_profiles(a,b;kwargs...), blocks)
end

"""Remove gap-causing representatives first, then internally gapped ones.
Duplicates are restored. If all candidates are flagged, a sibling species
reference breaks the tie; without one, retain the first candidate deterministically.
"""
function improve_alignment!(paths; print_messages=true, log_dir=nothing, kwargs...)
    for path in paths
        original = read_reference_records(path)
        isempty(original) && continue
        raw = [Sequence(r.id,replace(r.sequence,"-"=>"")) for r in original]
        reps = unique(r->r.sequence,raw)
        aligned = align_blocks(reps;kwargs...)
        while length(aligned)>1
            mat = alignment_matrix(aligned)
            bad = causes_internal_gaps(mat)
            any(bad) || (bad = inner_gaps(mat))
            any(bad) || break
            if all(bad)
                genus = first(split(basename(path),'_'))
                siblings = filter(p->p!=path && startswith(basename(p),genus*"_"),fasta_files(dirname(path)))
                winner = 1
                if !isempty(siblings)
                    refs = read_reference_records(first(siblings))
                    if !isempty(refs)
                        ref = Sequence("__external_reference__",replace(first(refs).sequence,"-"=>""))
                        counts = [count(==('-'), first(align_sequences([Sequence(r.id,replace(r.sequence,"-"=>"")),ref];kwargs...)).sequence) for r in aligned]
                        winner = argmin(counts)
                    end
                end
                bad[winner]=false
            end
            aligned = trim_columns(aligned[.!bad])
            if any(inner_gaps(alignment_matrix(aligned))) || any(causes_internal_gaps(alignment_matrix(aligned)))
                aligned = align_blocks([Sequence(r.id,replace(r.sequence,"-"=>"")) for r in aligned];kwargs...)
            end
        end
        lookup = Dict(replace(r.sequence,"-"=>"")=>r.sequence for r in aligned)
        out = [Sequence(r.id,lookup[r.sequence]) for r in raw if haskey(lookup,r.sequence)]
        write_fasta(out,path)
        kept = Set(r.id for r in out)
        append_log(log_dir,"removed_sequences.tsv",[(file=basename(path),taxon=splitext(basename(path))[1],sequence=r.id,
            n_start=length(raw),n_retained=length(out)) for r in raw if !(r.id in kept)])
        print_messages && @info "Improved alignment" path retained=length(out) original=length(raw)
    end
    paths
end
improve_alignment(args...;kwargs...)=improve_alignment!(args...;kwargs...)

function aligned_distance(a,b)
    ia,ib = findall(!=('-'),a),findall(!=('-'),b)
    (isempty(ia) || isempty(ib)) && return Inf
    lo,hi = max(first(ia),first(ib)),min(last(ia),last(ib))
    lo>hi && return Inf
    compared=different=0
    @inbounds for j in lo:hi
        a[j]=='-' && b[j]=='-' && continue
        compared+=1
        different+=a[j]!=b[j]
    end
    iszero(compared) ? Inf : different/compared
end

function bounded_aligned_distance(a,b,first_a,last_a,first_b,last_b,gap)
    lo,hi=max(first_a,first_b),min(last_a,last_b)
    lo>hi && return Inf
    compared=different=0
    @inbounds for k in lo:hi
        x,y=a[k],b[k]
        x==gap && y==gap && continue
        compared+=1
        different+=x!=y
    end
    iszero(compared) ? Inf : different/compared
end

# A function barrier specializes the distance loop on byte sequences.
function nearest_unique_distances(rows,counts,progress)
    n=length(rows)
    gap=convert(eltype(first(rows)),'-')
    starts=Vector{Int}(undef,n)
    stops=Vector{Int}(undef,n)
    nearest=fill(Inf,n)
    for i in eachindex(rows)
        a=rows[i]
        lo=findfirst(!=(gap),a)
        hi=findlast(!=(gap),a)
        starts[i]=isnothing(lo) ? 1 : lo
        stops[i]=isnothing(hi) ? 0 : hi
        # An all-gap duplicate still has no comparable positions.
        counts[i]>1 && !isnothing(lo) && (nearest[i]=0.0)
    end
    comparisons=0
    total=n*(n-1)÷2
    for i in 1:n-1
        for j in i+1:n
            # Neither exact minimum can improve below zero.
            if !(iszero(nearest[i]) && iszero(nearest[j]))
                d=bounded_aligned_distance(rows[i],rows[j],starts[i],stops[i],starts[j],stops[j],gap)
                nearest[i]=min(nearest[i],d)
                nearest[j]=min(nearest[j],d)
            end
            comparisons+=1
            if !isnothing(progress) && comparisons%1024==0
                progress(comparisons,total,n)
            end
        end
    end
    !isnothing(progress) && progress(comparisons,total,n)
    nearest
end

"""Exact nearest distances, sharing identical sequences and avoiding per-pair copies."""
function nearest_aligned_distances(records::AbstractVector{Sequence};progress=nothing)
    isempty(records) && return Float64[]
    width=length(first(records).sequence)
    all(r->length(r.sequence)==width,records) ||
        throw(ArgumentError("Aligned sequences must have equal lengths"))
    sequences=String[]
    counts=Int[]
    membership=Int[]
    lookup=Dict{String,Int}()
    for r in records
        index=get(lookup,r.sequence,0)
        if iszero(index)
            push!(sequences,r.sequence)
            push!(counts,0)
            index=length(sequences)
            lookup[r.sequence]=index
        end
        counts[index]+=1
        push!(membership,index)
    end
    # Sequence validates ASCII input; code units reference its strings without copies.
    rows=codeunits.(sequences)
    nearest=nearest_unique_distances(rows,counts,progress)
    nearest[membership]
end

"""Remove singleton outliers using nearest-neighbour uncorrected amino-acid distance.
Terminal gaps are excluded; if every sequence is flagged, the file is retained.
Identical sequences share comparisons; duplicates with nongap residues have distance zero.
Distances remain exact, including those written to removal logs.
`show_progress` reports file counts and pairwise comparisons every
`progress_interval` seconds (default 5), plus start/completion messages.
It defaults to `print_messages`, but can be enabled independently.
"""
function filter_divergent_seq!(paths; max_distance=0.10, min_seqs=4, print_messages=true,
                               log_dir=nothing, show_progress=print_messages, progress_interval=5.0)
    0<=max_distance<=1 || throw(ArgumentError("max_distance must be in [0,1]"))
    isfinite(progress_interval) && progress_interval>0 || throw(ArgumentError("progress_interval must be finite and positive"))
    total = length(paths)
    started = last_report = time()
    show_progress && @info "Filtering divergent sequences: started" files=total
    for (file_index,path) in enumerate(paths)
        if show_progress && (file_index==1 || file_index%100==0 || time()-last_report>=progress_interval)
            @info "Divergent sequence filter: file" file=file_index total_files=total path
            last_report = time()
        end
        records = read_reference_records(path)
        n = length(records)
        n<max(2,min_seqs) && continue
        progress=function(comparisons,total_comparisons,unique_sequences)
            if time()-last_report>=progress_interval
                percent=iszero(total_comparisons) ? 100.0 : round(100*comparisons/total_comparisons;digits=1)
                @info "Divergent sequence filter: comparisons" file=file_index total_files=total path sequences=n unique_sequences comparisons total_comparisons percent
                last_report=time()
            end
        end
        nearest=nearest_aligned_distances(records;progress=show_progress ? progress : nothing)
        bad = nearest .> max_distance
        any(bad) || continue
        if all(bad)
            @warn "Every sequence exceeds max_distance; retaining alignment" path
            continue
        end
        out = trim_columns(records[.!bad])
        write_fasta(out,path)
        append_log(log_dir,"removed_divergent.tsv",[(file=basename(path),taxon=splitext(basename(path))[1],
            sequence=records[i].id,distance=round(nearest[i];digits=4),n_start=n,n_retained=length(out)) for i in findall(bad)])
        print_messages && @info "Filtered divergent sequences" path retained=length(out) original=n
    end
    show_progress && @info "Filtering divergent sequences: completed" files=total elapsed_seconds=round(time()-started;digits=1)
    paths
end
filter_divergent_seq(args...;kwargs...)=filter_divergent_seq!(args...;kwargs...)
