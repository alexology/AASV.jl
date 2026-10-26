"""HTTP requester with shared exponential backoff. Permanent 4xx errors fail immediately.
`max_attempts=0` retries transient errors indefinitely, as in the R package.
"""
function make_retry_performer(; initial_wait=60.0, max_wait=7200.0, max_attempts=0,
                              transport=HTTP.get, sleeper=sleep)
    initial_wait>=0 && max_wait>=0 && max_attempts>=0 || throw(ArgumentError("Invalid retry parameters"))
    level = Ref(0)
    function perform(url; timeout=120)
        attempt=0
        while true
            attempt+=1
            level[]+=1
            try
                response = transport(url; readtimeout=timeout, connect_timeout=60, retry=false, status_exception=false)
                status=response.status
                if 200<=status<300
                    level[]=max(0,level[]-1)
                    return JSON3.read(String(response.body))
                end
                if 400<=status<500 && !(status in (408,429))
                    throw(ArgumentError("BOLD HTTP $status for $url"))
                end
                error("Transient HTTP $status")
            catch e
                e isa InterruptException && rethrow()
                e isa ArgumentError && rethrow()
                max_attempts>0 && attempt>=max_attempts && rethrow()
                pause=min(initial_wait*2.0^min(level[]-1,40),max_wait)
                @warn "BOLD request failed; retrying" attempt pause exception=e
                sleeper(pause)
            end
        end
    end
end

# Unlike Filesystem.rename on Julia 1.10, never fall back to copying/removing.
function rename_atomic(source::AbstractString, destination::AbstractString)
    err = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), source, destination)
    err < 0 && Base.uv_error("rename($(repr(source)), $(repr(destination)))", err)
    nothing
end

"""Write a checkpoint only after its entire JSON payload is on disk."""
function write_json_checkpoint(path, value)
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        JSON3.write(io,value)
        close(io)
        # A single rename replaces the file atomically; mv(force=true) could
        # remove the old checkpoint and then fail to move the new one.
        rename_atomic(temporary,path)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

"""Select one deterministic, one-based shard from a collection.
Useful for scheduler job arrays; every distinct value belongs to exactly one shard.
"""
function select_shard(values; shard_index::Integer=1, shard_count::Integer=1)
    shard_count >= 1 || throw(ArgumentError("shard_count must be positive"))
    1 <= shard_index <= shard_count ||
        throw(ArgumentError("shard_index must be between 1 and shard_count"))
    ordered=sort!(unique(collect(values));by=string)
    ordered[shard_index:shard_count:end]
end

"""Retrieve COI-5P BOLD records for a family using the paginated portal API.
`requester` and `sleeper` can be replaced for offline/reproducible workflows.
`cache_dir` retains the query metadata and complete pages across interruptions.
`overwrite=true` starts a fresh download; cached pages are never mixed across queries.
"""
function fetch_bold_records(family; timeout=7200, page_size=50000,
                            requester=make_retry_performer(), sleeper=sleep,
                            show_progress=false, cache_dir=nothing, overwrite=false)
    page_size>0 || throw(ArgumentError("page_size must be positive"))
    state_path=isnothing(cache_dir) ? nothing : joinpath(cache_dir,"download.json")
    overwrite && !isnothing(state_path) && rm(state_path;force=true)
    state=(!overwrite && !isnothing(state_path) && isfile(state_path)) ?
        JSON3.read(read(state_path,String)) : nothing
    if !isnothing(state) && (get(state,:version,0)!=1 || get(state,:family,nothing)!=family ||
                            get(state,:page_size,0)!=page_size)
        throw(ArgumentError("BOLD cache settings differ; use overwrite=true to start a fresh download"))
    end
    if isnothing(state)
        query = HTTP.escapeuri("tax:family:"*family)
        info = requester("https://portal.boldsystems.org/api/query?query=$query&extent=full";timeout=120)
        id=get(info,:query_id,nothing)
        isnothing(id) && return NamedTuple[]
        url="https://portal.boldsystems.org/api/documents/"*HTTP.escapeuri(string(id))
        total=0
        for attempt in 1:12
            total=Int(get(requester(url*"?length=0&start=0";timeout=120),:recordsTotal,0))
            total>0 && break
            attempt<12 && sleeper(5)
        end
        # Do not persist an unready query as an empty, completed download.
        if !isnothing(cache_dir) && total>0
            mkpath(cache_dir)
            pages=basename(mktempdir(cache_dir;prefix="pages_"))
            state=(version=1,family=family,page_size=page_size,query_id=string(id),total=total,pages=pages)
            write_json_checkpoint(state_path,state)
        end
    else
        id=state.query_id
        total=Int(state.total)
        url="https://portal.boldsystems.org/api/documents/"*HTTP.escapeuri(string(id))
        show_progress && @info "BOLD download: resuming saved query" family total
    end
    pages_dir=isnothing(state) ? nothing : joinpath(cache_dir,state.pages)
    result=NamedTuple[]
    start=0
    while start<total
        show_progress && @info "Downloading BOLD records" family received=start total
        page_path=isnothing(pages_dir) ? nothing : joinpath(pages_dir,"$start.json")
        cached=!isnothing(page_path) && isfile(page_path)
        page=cached ? JSON3.read(read(page_path,String)) : requester(url*"?length=$page_size&start=$start";timeout)
        records=get(page,:data,[])
        isempty(records) && throw(ErrorException("BOLD returned an empty page before recordsTotal was reached"))
        length(records)<=min(page_size,total-start) || throw(ErrorException("BOLD page exceeds the expected record count"))
        !cached && !isnothing(page_path) && write_json_checkpoint(page_path,page)
        cached && show_progress && @info "BOLD download: reusing saved page" family start records=length(records)
        for r in records
            marker=get(r,:marker_code,nothing)
            if isnothing(marker)
                marker=last(split(string(get(r,:record_id,"")),'.'))
            end
            marker=="COI-5P" || continue
            push!(result,(processid=get(r,:processid,nothing),identification=get(r,:identification,nothing),sequence=get(r,:nuc,nothing)))
        end
        start+=length(records)
        show_progress && @info "BOLD download progress" family received=start total
    end
    result
end

valid_species(s) = !isnothing(s) && !ismissing(s) && count(==(' '),s)==1 &&
    !occursin(r"cf\.|sp\.|/|[0-9]",s) && !isempty(s)

"""Download, filter, translate and align species references for all input families.
`records_provider(family)` may supply local BOLD-like records (`processid`,
`identification`, `sequence`). Outputs follow the original R directory layout.
`show_progress=true` reports families, downloads, filtering and each species
before alignment. Filtering updates use `progress_interval` seconds (default 5).
Complete BOLD pages are cached under `alignment_dir/bold_downloads/<family>`.
With `overwrite=false`, interrupted downloads resume from that cache;
`overwrite=true` starts a fresh query and download.
"""
function align_species_seq!(taxonomy; alignment_dir=pwd(), overwrite=false,
                           timeout=7200, min_length=500, max_length=1000,
                           genetic_code="5", records_provider=nothing,
                           show_progress=true, progress_interval=5.0,
                           threads=:auto, shard_index=1, shard_count=1, kwargs...)
    df=DataFrame(taxonomy)
    require_columns(df,(:family,:genus,:species))
    0<=min_length<=max_length || throw(ArgumentError("Invalid sequence length limits"))
    genetic_table(genetic_code)
    isfinite(progress_interval) && progress_interval>0 || throw(ArgumentError("progress_interval must be finite and positive"))
    thread_budget=threads===:auto ? Threads.nthreads() :
        threads isa Integer && threads>0 ? min(Int(threads),Threads.nthreads()) :
        throw(ArgumentError("threads must be :auto or a positive integer"))
    families=select_shard(collect(skipmissing(df.family));shard_index,shard_count)
    started=time()
    skipped=completed=0
    show_progress && @info "Species reference alignment: started" families=length(families)
    for (family_index,family) in enumerate(families)
        family=safe_name(string(family))
        marker=joinpath(alignment_dir,"raw_sequences",family,".complete")
        if isfile(marker) && !overwrite
            skipped+=1
            show_progress && @info "Species references: reusing completed family" family family_index total_families=length(families)
            continue
        end
        # Invalidate the previous completion before any work that may fail.
        rm(marker;force=true)
        family_started=time()
        show_progress && @info "Species references: loading records" family family_index total_families=length(families)
        records=isnothing(records_provider) ? fetch_bold_records(family;timeout,show_progress,
            cache_dir=joinpath(alignment_dir,"bold_downloads",family),overwrite) : records_provider(family)
        total_records=length(records)
        show_progress && @info "Species references: filtering and translating" family records=total_records
        last_report=time()
        groups=Dict{String,Vector{Sequence}}()
        for (record_index,r) in enumerate(records)
            if show_progress && record_index%1024==0 && time()-last_report>=progress_interval
                @info "Species references: filtering progress" family processed=record_index total_records species=length(groups)
                last_report=time()
            end
            valid_species(r.identification) || continue
            (isnothing(r.processid) || ismissing(r.processid) || isnothing(r.sequence) || ismissing(r.sequence)) && continue
            s=String(r.sequence)
            occursin(r"[^ACGTN+.-]",s) && continue
            # ShortRead::clean drops whole records containing non-ACGT bases.
            # Removing ambiguous bases would silently introduce a frameshift.
            occursin(r"[^ACGT]",s) && continue
            min_length<=length(s)<=max_length || continue
            translation=find_translation(s;genetic_code)
            (isnothing(translation) || isempty(translation.amino_acids)) && continue
            sp=safe_name(replace(String(r.identification)," "=>"_"))
            id=sp*"__"*string(r.processid)
            push!(get!(groups,sp,Sequence[]),Sequence(id,s))
        end
        species=sort(collect(keys(groups)))
        show_progress && @info "Species references: filtering completed" family records=total_records species=length(species)
        # Independent species use outer Julia tasks. Limit each batch to the
        # thread budget so completed alignments do not accumulate in memory.
        for batch_start in 1:thread_budget:length(species)
            batch_indices=batch_start:min(batch_start+thread_budget-1,length(species))
            batch_species=species[batch_indices]
            nt_groups=[unique(r->r.id,groups[sp]) for sp in batch_species]
            aa_groups=[[Sequence(r.id,something(find_translation(r.sequence;genetic_code)).amino_acids)
                for r in nt] for nt in nt_groups]
            for (offset,(sp,nt)) in enumerate(zip(batch_species,nt_groups))
                species_index=batch_start+offset-1
                show_progress && @info "Species references: aligning species" family species=sp species_index total_species=length(species) sequences=length(nt)
            end
            aligned=ProfileAlignments.align_sequence_sets(aa_groups;threads=thread_budget,kwargs...)
            for (sp,nt,aa,result) in zip(batch_species,nt_groups,aa_groups,aligned)
                write_fasta(nt,joinpath(alignment_dir,"raw_sequences_nt",family,"species",sp*".fasta"))
                write_fasta(aa,joinpath(alignment_dir,"raw_sequences",family,"species",sp*".fasta"))
                write_fasta(result.records,joinpath(alignment_dir,"alignments",family,"species",sp*".fasta"))
                write_table(joinpath(alignment_dir,"aa_tables",family,sp*".txt"),aa_table(nt;genetic_code))
            end
        end
        if !isempty(groups)
            mkpath(dirname(marker))
            write(marker,"completed\n")
        else
            @warn "No usable references found" family
        end
        completed+=1
        show_progress && @info "Species references: family completed" family family_index total_families=length(families) species=length(species) elapsed_seconds=round(time()-family_started;digits=1)
    end
    show_progress && @info "Species reference alignment: completed" families=length(families) processed=completed reused=skipped elapsed_seconds=round(time()-started;digits=1)
    nothing
end
align_species_seq(args...;kwargs...)=align_species_seq!(args...;kwargs...)

# Records with insertions open columns that are gaps in most references; a
# query is then flagged as an indel against the whole profile. Remove them
# iteratively, unless every record is flagged.
function prune_gap_causing(records)
    removed=Sequence[]
    while length(records)>1
        bad=causes_internal_gaps(alignment_matrix(records))
        (any(bad) && !all(bad)) || break
        append!(removed,records[bad])
        records=trim_columns(records[.!bad])
    end
    (records=records,removed=removed)
end

function write_gap_causing_report(path, n_start, records, removed)
    rm(path;force=true)
    isempty(removed) && return
    mkpath(dirname(path))
    CSV.write(path,DataFrame(sequence=[r.id for r in removed],n_start=fill(n_start,length(removed)),
        n_retained=fill(length(records),length(removed)));delim='	')
end

function merge_reference_level(level; taxon=nothing, alignment_dir=pwd(), overwrite=false,
                               show_progress=true, shard_index=1, shard_count=1,
                               remove_gap_causing=true, log_dir=joinpath(alignment_dir,"logs"), kwargs...)
    level=rank_name(parse_rank(level))
    base=joinpath(alignment_dir,"alignments")
    families=isnothing(taxon) ? (isdir(base) ? filter(f->isdir(joinpath(base,f)),readdir(base)) : String[]) :
        taxon isa AbstractString ? [taxon] : taxon
    families=select_shard(families;shard_index,shard_count)
    sourcelevel=level=="genus" ? "species" : "genus"
    started=time()
    completed=reused=0
    show_progress && @info "Reference profile alignment: started" level families=length(families)
    for (family_index,family) in enumerate(families)
        family=safe_name(string(family))
        show_progress && @info "Reference profiles: loading family" level family family_index total_families=length(families)
        paths=fasta_files(joinpath(base,family,sourcelevel))
        groups=Dict{String,Vector{String}}()
        for path in paths
            name=level=="genus" ? first(split(basename(path),'_')) : family
            push!(get!(groups,name,String[]),path)
        end
        names=sort(collect(keys(groups)))
        for (taxon_index,name) in enumerate(names)
            dest=joinpath(base,family,level,name*".fasta")
            if isfile(dest) && !overwrite
                reused+=1
                show_progress && @info "Reference profiles: reusing completed taxon" level family taxon=name taxon_index total_taxa=length(names)
                continue
            end
            taxon_started=time()
            show_progress && @info "Reference profiles: loading taxon" level family taxon=name taxon_index total_taxa=length(names) files=length(groups[name])
            blocks=read_reference_records.(groups[name])
            sort!(blocks;by=length,rev=true,alg=MergeSort)
            # Earlier aasv versions could place a species in two genus profiles;
            # its records are kept once, in the largest profile.
            seen=Dict{String,String}()
            repeats=0
            blocks=map(blocks) do block
                kept=filter(block) do r
                    raw=replace(r.sequence,"-"=>"")
                    previous=get(seen,r.id,nothing)
                    isnothing(previous) && (seen[r.id]=raw; return true)
                    previous==raw || throw(ArgumentError("Reference ID $(r.id) has different sequences in the $name profiles"))
                    false
                end
                repeats+=length(block)-length(kept)
                length(kept)<length(block) ? trim_columns(kept) : block
            end
            repeats>0 && @warn "Records repeated across reference profiles were merged once" level family taxon=name repeats
            filter!(!isempty,blocks)
            if isempty(blocks)
                show_progress && @info "Reference profiles: no usable profiles" level family taxon=name
                continue
            end
            merges_done=0
            total_merges=length(blocks)-1
            result=reduce(blocks) do a,b
                show_progress && @info "Reference profiles: merging" level family taxon=name merge_index=merges_done+1 total_merges left_sequences=length(a) right_sequences=length(b)
                merged=align_profiles(a,b;kwargs...)
                merges_done+=1
                show_progress && @info "Reference profiles: merge completed" level family taxon=name merges_done total_merges percent=round(100*merges_done/total_merges;digits=1)
                merged
            end
            if remove_gap_causing
                n_start=length(result)
                result,removed=prune_gap_causing(result)
                isnothing(log_dir) || write_gap_causing_report(joinpath(log_dir,"gap_causing_references",
                    level,family,name*".tsv"),n_start,result,removed)
                show_progress && !isempty(removed) && @info "Reference profiles: removed gap-causing records" level family taxon=name removed=length(removed) retained=length(result)
            end
            write_fasta(result,dest)
            completed+=1
            show_progress && @info "Reference profiles: taxon completed" level family taxon=name taxon_index total_taxa=length(names) sequences=length(result) elapsed_seconds=round(time()-taxon_started;digits=1)
        end
        show_progress && @info "Reference profiles: family completed" level family family_index total_families=length(families)
    end
    show_progress && @info "Reference profile alignment: completed" level processed=completed reused elapsed_seconds=round(time()-started;digits=1)
    nothing
end

"""Merge species profiles into genus FASTA alignments, largest profile first.
`show_progress=true` reports families, taxa and each profile merge.
After merging, records whose insertions open gap columns in most references
are removed iteratively (`remove_gap_causing=true`; R aasv relies on
`check_alignment` and `improve_alignment` at each level), and listed
in `log_dir/gap_causing_references/<level>/<family>/<taxon>.tsv`
(`log_dir=nothing` disables the report). Otherwise a few such records make
every query an indel at that level. Existing profiles are reused unless
`overwrite=true`.
"""
align_genus_seq!(;kwargs...) = merge_reference_level(GenusRank;kwargs...)
align_genus_seq(;kwargs...) = align_genus_seq!(;kwargs...)
"""Merge genus profiles into family FASTA alignments, largest profile first.
`show_progress=true` reports families, taxa and each profile merge.
After merging, records whose insertions open gap columns in most references
are removed iteratively (`remove_gap_causing=true`; R aasv relies on
`check_alignment` and `improve_alignment` at each level), and listed
in `log_dir/gap_causing_references/<level>/<family>/<taxon>.tsv`
(`log_dir=nothing` disables the report). Otherwise a few such records make
every query an indel at that level. Existing profiles are reused unless
`overwrite=true`.
"""
align_family_seq!(;kwargs...) = merge_reference_level(FamilyRank;kwargs...)
align_family_seq(;kwargs...) = align_family_seq!(;kwargs...)
