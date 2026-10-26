# Bump these versions when the corresponding analysis algorithm changes.
const ALIGNMENT_CHECK_VERSION = 1
const DIVERGENT_FILTER_VERSION = 1

file_digest(path) = open(sha256,path) |> bytes2hex
reference_fingerprint(files) = Dict(basename(path)=>file_digest(path) for path in files)

function reference_families(alignment_dir, type, families)
    root=joinpath(alignment_dir,type)
    selected=isnothing(families) ? (isdir(root) ? readdir(root) : String[]) :
        families isa AbstractString ? [families] : collect(families)
    sort!(unique!([safe_name(string(f)) for f in selected if isdir(joinpath(root,string(f)))]))
end

function read_reference_checkpoint(path)
    isfile(path) || return nothing
    try
        JSON3.read(read(path,String))
    catch e
        e isa InterruptException && rethrow()
        @warn "Unreadable reference checkpoint; repeating analysis" path exception=e
        nothing
    end
end

function checkpoint_matches(state, version, parameters, fingerprint)
    !isnothing(state) && get(state,:version,0)==version &&
        JSON3.write(get(state,:parameters,nothing))==JSON3.write(parameters) &&
        Dict(String(k)=>String(v) for (k,v) in pairs(state.fingerprint))==fingerprint
end

"""Check all reference families and return problematic FASTA paths, including cached results.
Checkpoints in `alignment_dir/checkpoints/check_alignment` are reused only when
file names, SHA-256 contents, rank, type and algorithm version match. Set
`resume=false` to repeat checks. Missing/empty rank directories are skipped.
"""
function check_alignments!(; alignment_dir=pwd(), tax_lev=SpeciesRank,
                            type="alignments", families=nothing, resume=true,
                            show_progress=true, shard_index=1, shard_count=1)
    type in ("alignments","raw_sequences") || throw(ArgumentError("Invalid type"))
    rank=rank_name(parse_rank(tax_lev))
    parameters=(type=type,rank=rank)
    selected=select_shard(reference_families(alignment_dir,type,families);shard_index,shard_count)
    bad=String[]
    for (index,family) in enumerate(selected)
        files=fasta_files(joinpath(alignment_dir,type,family,rank))
        isempty(files) && continue
        checkpoint=joinpath(alignment_dir,"checkpoints","check_alignment",type,rank,family*".json")
        fingerprint=reference_fingerprint(files)
        state=resume ? read_reference_checkpoint(checkpoint) : nothing
        reused=checkpoint_matches(state,ALIGNMENT_CHECK_VERSION,parameters,fingerprint)
        if reused
            append!(bad,[joinpath(dirname(first(files)),String(name)) for name in state.bad])
        else
            result=check_alignment(family;alignment_dir,type,tax_lev)
            write_json_checkpoint(checkpoint,(version=ALIGNMENT_CHECK_VERSION,
                parameters=parameters,fingerprint=fingerprint,bad=basename.(result)))
            append!(bad,result)
        end
        show_progress && @info "Reference alignment check: family completed" family index total=length(selected) reused problems=length(bad)
    end
    bad
end

# A prepared filter transaction retains its outputs and report until every file
# has been installed. After interruption, each file may have its input or output
# hash; either is safe to resume without filtering the same records twice.
function install_reference_filter!(pending, directory, stage)
    before=Dict(String(k)=>String(v) for (k,v) in pairs(pending.before))
    after=Dict(String(k)=>String(v) for (k,v) in pairs(pending.fingerprint))
    files=fasta_files(directory)
    current=reference_fingerprint(files)
    issubset(Set(keys(current)),Set(keys(after))) || return false
    all(current[name] in (before[name],after[name]) for name in keys(current)) || return false
    all(isfile(joinpath(stage,name)) && file_digest(joinpath(stage,name))==after[name]
        for name in keys(after)) || return false
    for name in sort!(collect(keys(after)))
        get(current,name,nothing)==after[name] && continue
        path=joinpath(directory,name)
        temporary,io=mktemp(directory)
        close(io)
        try
            cp(joinpath(stage,name),temporary;force=true)
            # Atomic replacement: on failure the previous file stays in place.
            rename_atomic(temporary,path)
        finally
            isfile(temporary) && rm(temporary)
        end
    end
    true
end

function write_reference_filter_report(path, rows)
    # A report belongs to one completed family, so resumes never append duplicates.
    table=DataFrame(file=String[String(row.file) for row in rows],
        taxon=String[String(row.taxon) for row in rows],
        sequence=String[String(row.sequence) for row in rows],
        distance=Float64[Float64(row.distance) for row in rows],
        n_start=Int[Int(row.n_start) for row in rows],
        n_retained=Int[Int(row.n_retained) for row in rows])
    mkpath(dirname(path))
    CSV.write(path,table;delim='\t')
end

"""Filter all reference families with persistent SHA-256 checkpoints.
Returns the current FASTA paths. `resume=true` reuses completed families only
when their output files, `max_distance`, `min_seqs` and algorithm version match.
Each family's removal report is retained in its checkpoint and written to
`log_dir/divergent_references/<rank>/<family>.tsv`. Filtering operates on the
current files; changing thresholds cannot restore previously removed records.
Prepared outputs allow interrupted installation to resume without filtering twice.
"""
function filter_divergent_references!(; alignment_dir=pwd(), tax_lev=SpeciesRank,
        families=nothing, max_distance=0.10, min_seqs=4, resume=true,
        log_dir=joinpath(alignment_dir,"logs"), print_messages=false,
        show_progress=true, progress_interval=5.0, shard_index=1, shard_count=1)
    0<=max_distance<=1 || throw(ArgumentError("max_distance must be in [0,1]"))
    isfinite(progress_interval) && progress_interval>0 || throw(ArgumentError("progress_interval must be finite and positive"))
    rank=rank_name(parse_rank(tax_lev))
    parameters=(rank=rank,max_distance=Float64(max_distance),min_seqs=max(2,min_seqs))
    selected=select_shard(reference_families(alignment_dir,"alignments",families);shard_index,shard_count)
    result=String[]
    for (index,family) in enumerate(selected)
        directory=joinpath(alignment_dir,"alignments",family,rank)
        files=fasta_files(directory)
        cache=joinpath(alignment_dir,"checkpoints","filter_divergent",rank,family)
        checkpoint=joinpath(cache,"complete.json")
        pending_path=joinpath(cache,"pending.json")
        isempty(files) && !isfile(pending_path) && continue
        fingerprint=reference_fingerprint(files)
        state=resume ? read_reference_checkpoint(checkpoint) : nothing
        reused=checkpoint_matches(state,DIVERGENT_FILTER_VERSION,parameters,fingerprint)
        if !reused
            pending=resume ? read_reference_checkpoint(pending_path) : nothing
            prepared=!isnothing(pending) && get(pending,:version,0)==DIVERGENT_FILTER_VERSION &&
                JSON3.write(pending.parameters)==JSON3.write(parameters)
            stage=prepared ? joinpath(cache,safe_name(String(pending.stage))) : ""
            installed=prepared && install_reference_filter!(pending,directory,stage)
            if !installed
                # Never leave an old completion marker valid while replacing outputs.
                rm(checkpoint;force=true)
                rm(pending_path;force=true)
                mkpath(cache)
                stage=mktempdir(cache;prefix="prepared_",cleanup=false)
                for path in files
                    cp(path,joinpath(stage,basename(path)))
                end
                filter_divergent_seq!(fasta_files(stage);max_distance,min_seqs,
                    log_dir=stage,print_messages,show_progress,progress_interval)
                report_path=joinpath(stage,"removed_divergent.tsv")
                rows=isfile(report_path) ? NamedTuple.(eachrow(CSV.read(report_path,DataFrame;delim='\t'))) : NamedTuple[]
                pending=(version=DIVERGENT_FILTER_VERSION,parameters=parameters,
                    before=fingerprint,fingerprint=reference_fingerprint(fasta_files(stage)),
                    removed=rows,stage=basename(stage))
                write_json_checkpoint(pending_path,pending)
                install_reference_filter!(pending,directory,stage) ||
                    error("Reference files changed during filtering of $family; rerun the analysis")
            end
            state=pending
            write_json_checkpoint(checkpoint,state)
            rm(pending_path;force=true)
            # Only remove the known files in this generated staging directory.
            for path in fasta_files(stage)
                rm(path)
            end
            rm(joinpath(stage,"removed_divergent.tsv");force=true)
            isempty(readdir(stage)) && rm(stage)
        end
        if !isnothing(log_dir)
            write_reference_filter_report(joinpath(log_dir,"divergent_references",rank,family*".tsv"),state.removed)
        end
        append!(result,fasta_files(directory))
        show_progress && @info "Reference divergence filter: family completed" family index total=length(selected) reused
    end
    result
end

"""Report reference records that `asv_functional_structure!` ignores or excludes.
Returns a table with `family`, `tax_lev`, `file`, `record` and `issue`:
`"repeated"` (exact repeat of a record), `"conflicting"` (ID repeated with
another sequence; analyses stop on it) or `"no_nucleotide_reference"` (profile
record without a nucleotide reference in the taxon, such as another genus in a
genus profile or an altered ID written by an earlier aasv version).
Genus and family profiles can be rebuilt with `overwrite=true`. Realigning
species files from the raw sequences would restore records removed by earlier
curation; rebuild them with `align_species_seq!(overwrite=true)` instead.
"""
function check_reference_records(; alignment_dir=pwd(), families=nothing,
                                 tax_lev=(SpeciesRank,GenusRank,FamilyRank))
    ranks=tax_lev isa Union{TaxonomicRank,Symbol,AbstractString} ? [tax_lev] : collect(tax_lev)
    levels=[rank_name(parse_rank(rank)) for rank in ranks]
    out=DataFrame(family=String[],tax_lev=String[],file=String[],record=String[],issue=String[])
    for family in reference_families(alignment_dir,"alignments",families)
        nt_dir=joinpath(alignment_dir,"raw_sequences_nt",family,"species")
        nt_ids=Dict(path=>Set(first.(fasta_entries(path))) for path in fasta_files(nt_dir))
        for level in levels, path in fasta_files(joinpath(alignment_dir,"alignments",family,level))
            name=splitext(basename(path))[1]
            ids=reduce(union,(nt_ids[p] for p in taxon_nt_paths(nt_dir,level,name));init=Set{String}())
            seen=Dict{String,String}()
            for (id,sequence) in fasta_entries(path)
                previous=get(seen,id,nothing)
                issue=if isnothing(previous)
                    seen[id]=sequence
                    id in ids ? nothing : "no_nucleotide_reference"
                else
                    previous==sequence ? "repeated" : "conflicting"
                end
                isnothing(issue) || push!(out,(family,level,path,id,issue))
            end
        end
    end
    out
end
