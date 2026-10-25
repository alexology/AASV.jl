const KYTE_DOOLITTLE = Dict(zip(collect("ACDEFGHIKLMNPQRSTVWY"),
    [1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]))
const GRANTHAM_ORDER = collect("SRLPTAVGIF YCHQNKDEMW" |> s->replace(s," "=>""))
const GRANTHAM_PATH = joinpath(@__DIR__,"..","data","grantham.csv")
include_dependency(GRANTHAM_PATH)
const GRANTHAM_MATRIX = let
    rows = readlines(GRANTHAM_PATH)[2:end]
    reduce(vcat,[permutedims(parse.(Int,split(r,',')[2:end])) for r in rows])
end

"""Typed result of a Kyte-Doolittle calculation."""
struct HydrophobicityMetrics
    global_score::Union{Missing,Float64}
    profile::Vector{Union{Missing,Float64}}
end

@enum StructuralVerdict begin
    StructureUnavailable
    StructurePass
    StructureFail
    InvalidORF
end

"""Typed structural assessment used by the computational core."""
struct StructureAssessment
    global_min::Union{Missing,Float64}
    global_max::Union{Missing,Float64}
    global_score::Union{Missing,Float64}
    global_pass::Union{Missing,Bool}
    local_violation_rate::Union{Missing,Float64}
    verdict::StructuralVerdict
end

"""One amino-acid difference with nucleotide and reference coordinates."""
struct MutationObservation
    aa_position::Int
    query_aa::Char
    query_codon::String
    query_nt::UnitRange{Int}
    aligned_sequences::Int
    unique_sequences::Int
    quality_probability::Union{Nothing,Float64}
    quality_sequences::Union{Nothing,Int}
    reference_aa::Char
    reference_codon::String
    reference_nt::UnitRange{Int}
    reference_id::String
    mutated_codon_positions::Vector{Int}
    hamming_distance::Int
    grantham_distance::Union{Nothing,Float64}
end

@enum ASVStatus ASVValid ASVInvalidORF

"""Typed, in-memory result of an ASV audit."""
struct ASVAssessment
    status::ASVStatus
    structure::StructureAssessment
    mutations::Vector{MutationObservation}
end

verdict_label(v::StructuralVerdict)=v===StructurePass ? "PASS" :
    v===StructureFail ? "FAIL" : v===InvalidORF ? "INVALID_ORF" : "UNAVAILABLE"

function assessment_table(a::StructureAssessment)
    DataFrame(global_min=[a.global_min],global_max=[a.global_max],
        global_score=[a.global_score],global_pass=[a.global_pass],
        local_violation_rate=[a.local_violation_rate],
        structural_verdict=[verdict_label(a.verdict)])
end

function grantham_distance(a::Char,b::Char)
    i,j = findfirst(==(a),GRANTHAM_ORDER),findfirst(==(b),GRANTHAM_ORDER)
    isnothing(i) || isnothing(j) ? missing : GRANTHAM_MATRIX[i,j]
end
get_mut_pos(q,r) = join(findall(i->q[i]!=r[i],1:3),",")

"""Kyte–Doolittle mean and sliding-window profile. Unknown residues yield missing."""
function get_hydro_metrics(aa_seq::AbstractString; hydro_window=7)
    hydro_window>0 || throw(ArgumentError("hydro_window must be positive"))
    values=[get(KYTE_DOOLITTLE,c,missing) for c in aa_seq]
    global_score=isempty(values) ? missing : mean(values)
    profile=[mean(values[i:i+hydro_window-1]) for i in 1:max(length(values)-hydro_window+1,0)]
    HydrophobicityMetrics(global_score,profile)
end
get_hydro_metrics(aa_seq::BioSequences.AASeq;kwargs...)=get_hydro_metrics(String(aa_seq);kwargs...)

"""Compare query hydrophobicity against reference bounds without table allocation."""
function assess_structure(query_aa,global_bounds,lower_env,upper_env;
                          violation_threshold=0.1,hydro_window=7)
    0<=violation_threshold<=1 || throw(ArgumentError("Invalid violation threshold"))
    q=get_hydro_metrics(query_aa;hydro_window)
    length(q.profile)==length(lower_env)==length(upper_env) || throw(ArgumentError("Envelope length mismatch"))
    global_pass=ismissing(q.global_score) ? missing : global_bounds[1]<=q.global_score<=global_bounds[2]
    usable=[i for i in eachindex(q.profile) if !ismissing(q.profile[i]) && !ismissing(lower_env[i]) && !ismissing(upper_env[i])]
    rate=isempty(usable) ? missing : count(i->!(lower_env[i]<=q.profile[i]<=upper_env[i]),usable)/length(usable)
    verdict=ismissing(global_pass) || ismissing(rate) ? StructureUnavailable :
        global_pass && rate<violation_threshold ? StructurePass : StructureFail
    StructureAssessment(Float64(global_bounds[1]),Float64(global_bounds[2]),
        q.global_score,global_pass,rate,verdict)
end

"""Tabular compatibility adapter for `assess_structure`."""
check_asv_structure(args...;kwargs...)=assessment_table(assess_structure(args...;kwargs...))

function empty_aa()
    DataFrame(aa_pos=Int[],illegal_aa=Union{Missing,String}[],dna_triplet=Union{Missing,String}[],
        nt_coords=Union{Missing,String}[],n_sequences_ali=Union{Missing,Int}[],n_unique_seq_ali=Union{Missing,Int}[],
        median_quality_prob=Union{Missing,Float64}[],n_sequences_prob=Union{Missing,Int}[],
        ref_aa=Union{Missing,String}[],dna_triplet_ref=Union{Missing,String}[],nt_coords_ref=Union{Missing,String}[],
        ref_id=Union{Missing,String}[],triplet_mut_pos=Union{Missing,String}[],ham_dist=Union{Missing,Int}[],
        grantham_dist=Union{Missing,Float64}[])
end

function mutation_table(observations::AbstractVector{MutationObservation})
    out=empty_aa()
    for x in observations
        push!(out,(x.aa_position,string(x.query_aa),x.query_codon,
            "$(first(x.query_nt))-$(last(x.query_nt))",x.aligned_sequences,x.unique_sequences,
            something(x.quality_probability,missing),something(x.quality_sequences,missing),
            string(x.reference_aa),x.reference_codon,
            "$(first(x.reference_nt))-$(last(x.reference_nt))",x.reference_id,
            join(x.mutated_codon_positions,','),x.hamming_distance,
            something(x.grantham_distance,missing)))
    end
    out
end

function assessment_tables(result::ASVAssessment)
    aa=mutation_table(result.mutations)
    if result.status===ASVInvalidORF
        push!(aa,(-1,missing,missing,missing,missing,missing,missing,missing,
            missing,missing,missing,missing,missing,missing,missing))
    end
    (hydro=assessment_table(result.structure),aa=aa)
end

function unavailable_assessment(verdict=StructureUnavailable)
    StructureAssessment(missing,missing,missing,missing,missing,verdict)
end
unavailable_hydro(verdict="UNAVAILABLE")=assessment_table(unavailable_assessment(
    verdict=="INVALID_ORF" ? InvalidORF : StructureUnavailable))

function no_orf_result()
    assessment_tables(ASVAssessment(ASVInvalidORF,
        unavailable_assessment(InvalidORF),MutationObservation[]))
end

function quality_probabilities(path,ee)
    qualities=Vector{Vector{Float64}}()
    open(path) do io
        while !eof(io)
            header=readline(io)
            startswith(header,"@") || throw(ArgumentError("Invalid FASTQ header"))
            sequence=readline(io)
            startswith(readline(io),"+") || throw(ArgumentError("Invalid FASTQ separator"))
            quality=readline(io)
            length(sequence)==length(quality) || throw(ArgumentError("FASTQ quality length mismatch"))
            all(c->33<=Int(c)<=126,quality) || throw(ArgumentError("Invalid Phred+33 quality"))
            errors=[10.0^(-(Int(c)-33)/10) for c in quality]
            sum(errors)<=ee && push!(qualities,errors)
        end
    end
    qualities
end

function mapped_qualities(query_dna,query_id,output_dir,vsearch_path,trimmed_folder,asv_samples,ee)
    dir=joinpath(output_dir,"vsearch_map_back")
    mkpath(dir)
    query=joinpath(dir,query_id*".fasta")
    write_fasta([Sequence(query_id,query_dna)],query)
    paths=sort(filter(p->isfile(p) && (endswith(lowercase(p),".fastq") || endswith(lowercase(p),".fq")) &&
        any(s->occursin(s,basename(p)),asv_samples),readdir(trimmed_folder;join=true)))
    qualities=Vector{Vector{Float64}}()
    for (i,path) in enumerate(paths)
        tab=joinpath(dir,"$(query_id)_$(i)_mapping_results.tsv")
        fq=joinpath(dir,"$(query_id)_$(i)_mapping_results.fastq")
        ids=joinpath(dir,"$(query_id)_$(i)_clean_ids.txt")
        run(`$vsearch_path --usearch_global $query --db $path --id 1.0 --strand plus --maxaccepts 0 --maxrejects 0 --userout $tab --userfields target --notrunclabels`)
        isfile(tab) || error("VSEARCH did not produce mapping output")
        labels=unique(filter(!isempty,strip.(readlines(tab))))
        isempty(labels) && continue
        open(ids,"w") do io
            foreach(label->println(io,label),labels)
        end
        run(`$vsearch_path --fastx_getseqs $path --label_words $ids --fastqout $fq --notrunclabels`)
        append!(qualities,quality_probabilities(fq,ee))
    end
    qualities
end

"""Audit a query against an amino-acid reference profile; report novel residues,
codon coordinates, Grantham distances, ORF incompatibilities and hydrophobicity.
Optional VSEARCH mapping uses Phred+33 error probabilities for read-back quality.
"""
function assess_asv(;query_dna,query_id,ref_aa_aligned,ref_nt,
                    output_dir=nothing,tax_lev=SpeciesRank,genetic_code="5",
                    hydro_threshold=0.1,hydro_window=7,vsearch_path=nothing,
                    trimmed_folder=nothing,asv_samples=String[],ee=40,kwargs...)
    hydro_window>0 || throw(ArgumentError("hydro_window must be positive"))
    0<=hydro_threshold<=1 || throw(ArgumentError("Invalid hydro_threshold"))
    query_id=safe_name(string(query_id))
    tax_lev=rank_name(parse_rank(tax_lev))
    query_dna=uppercase(String(query_dna))
    refs=unique(r->r.sequence,ref_aa_aligned)
    isempty(refs) && throw(ArgumentError("Reference profile cannot be empty"))
    translation=find_translation(query_dna;genetic_code)
    (isnothing(translation) || isempty(translation.amino_acids)) &&
        return ASVAssessment(ASVInvalidORF,unavailable_assessment(InvalidORF),MutationObservation[])
    aa=translation.amino_acids
    frame=translation.frame
    options=merge((mode=:query,),NamedTuple(kwargs))
    aln=align_profiles([Sequence(query_id,aa)],refs;options...)
    isnothing(output_dir) || write_fasta(aln,joinpath(output_dir,"query_alignments",query_id*"_"*safe_name(tax_lev)*".fasta"))
    mat=alignment_matrix(aln)
    (inner_gaps(mat)[1] || causes_internal_gaps(mat)[1]) &&
        return ASVAssessment(ASVInvalidORF,unavailable_assessment(InvalidORF),MutationObservation[])
    q=mat[1,:]
    refmat=mat[2:end,:]
    covered=findall(!=('-'),q)
    shared=[j for j in covered if count(==('-'),refmat[:,j])/size(refmat,1)<0.5]
    complete=[i for i in axes(refmat,1) if all(!=('-'),refmat[i,shared])]
    hydro=unavailable_assessment()
    if !isempty(complete) && length(shared)>=hydro_window
        metrics=[get_hydro_metrics(String(refmat[i,shared]);hydro_window) for i in complete]
        if all(x->!ismissing(x.global_score) && !any(ismissing,x.profile),metrics)
            bounds=extrema([x.global_score for x in metrics])
            lower=[quantile([x.profile[j] for x in metrics],0.025) for j in eachindex(first(metrics).profile)]
            upper=[quantile([x.profile[j] for x in metrics],0.975) for j in eachindex(first(metrics).profile)]
            hydro=assess_structure(String(q[shared]),bounds,lower,upper;
                violation_threshold=hydro_threshold,hydro_window)
        end
    end
    nt_lookup=Dict(r.id=>r.sequence for r in ref_nt)
    # Only references at novel residues are translated, not the whole taxon.
    ref_translations=Dict{String,Union{Nothing,Translation}}()
    illegal=[j for j in covered if !(q[j] in refmat[:,j])]
    qualities=Vector{Vector{Float64}}()
    if !isnothing(vsearch_path)
        (isnothing(trimmed_folder) || isnothing(output_dir)) && throw(ArgumentError("VSEARCH requires trimmed_folder and output_dir"))
        isempty(illegal) || (qualities=mapped_qualities(query_dna,query_id,output_dir,vsearch_path,trimmed_folder,asv_samples,ee))
    end
    observations=MutationObservation[]
    for j in illegal
        codon_index=count(!=('-'),q[1:j])
        nt_start=frame+3*(codon_index-1)
        triplet=query_dna[nt_start:nt_start+2]
        codon_probs=[clamp(1-sum(v[nt_start:nt_start+2]),0,1) for v in qualities if length(v)>=nt_start+2]
        prob=isempty(codon_probs) ? nothing : median(codon_probs)
        for (i,r) in enumerate(refs)
            c=refmat[i,j]
            c=='-' && continue
            haskey(nt_lookup,r.id) || throw(ArgumentError("Missing nucleotide reference for $(r.id)"))
            ref_translation=get!(()->find_translation(nt_lookup[r.id];genetic_code),ref_translations,r.id)
            isnothing(ref_translation) && throw(ArgumentError("Reference $(r.id) has no valid translation"))
            k=count(!=('-'),refmat[i,1:j])
            start=ref_translation.frame+3*(k-1)
            dna=nt_lookup[r.id]
            start+2<=length(dna) || throw(ArgumentError("Reference AA/DNA coordinates disagree"))
            ref_triplet=dna[start:start+2]
            translate_codon(ref_triplet,genetic_table(genetic_code))==c || error("Reference AA/DNA translation disagrees for $(r.id)")
            grantham=grantham_distance(q[j],c)
            push!(observations,MutationObservation(j,q[j],triplet,nt_start:nt_start+2,
                length(refs),length(refs),prob,isempty(codon_probs) ? nothing : length(codon_probs),
                c,ref_triplet,start:start+2,r.id,
                findall(k->triplet[k]!=ref_triplet[k],1:3),
                count(k->triplet[k]!=ref_triplet[k],1:3),
                ismissing(grantham) ? nothing : Float64(grantham)))
        end
    end
    ASVAssessment(ASVValid,hydro,observations)
end

"""Tabular compatibility adapter for `assess_asv`, used by report writers."""
calculate_asv_distance(;kwargs...)=assessment_tables(assess_asv(;kwargs...))

# Nucleotide reference files of the species in a taxon.
taxon_nt_paths(nt_dir,level,name)=level=="species" ? filter(isfile,[joinpath(nt_dir,name*".fasta")]) :
    level=="genus" ? filter(p->startswith(basename(p),name*"_"),fasta_files(nt_dir)) : fasta_files(nt_dir)

"""Load a taxon profile and its nucleotide references once per run.
Profile records need their nucleotide reference to report codons. Records
without one, such as other genera in genus files of earlier aasv versions or
altered IDs, are excluded with a warning; gap-only columns are then removed.
"""
function taxon_references(aa_path,nt_paths)
    loaded=read_references(aa_path)
    nt=Sequence[]
    repeats=loaded.repeats
    for path in nt_paths
        x=read_references(path)
        append!(nt,x.records)
        repeats+=x.repeats
    end
    repeats>0 && @warn "Exact repeats of reference records were ignored" path=aa_path repeats
    ids=Set(r.id for r in nt)
    refs=filter(r->r.id in ids,loaded.records)
    if length(refs)<length(loaded.records)
        excluded=[r.id for r in loaded.records if !(r.id in ids)]
        @warn "Profile records without a nucleotide reference in the taxon were excluded" path=aa_path excluded=length(excluded) examples=first(excluded,5)
        isempty(refs) || (refs=trim_columns(refs))
    end
    (refs=refs,nt=nt)
end

function add_metadata(df,row,level)
    out=copy(df)
    insertcols!(out,1,:ASV_id=>row.ASV_id,:tax_lev=>level,:species=>row.species)
    out
end

"""Analyze every available species/genus/family reference for an ASV taxonomy table.
Required columns: family, genus, species, ASV. Returns sorted taxonomy with ASV IDs.
Writes Excel reports and query FASTA alignments; numeric columns are sample abundances.
"""
function asv_functional_structure!(asv_taxonomy;output_dir,alignment_dir,vsearch_path=nothing,
                                  trimmed_folder=nothing,ee=40,genetic_code="5",
                                  hydro_window=7,hydro_threshold=0.1,
                                  threads=:auto,kwargs...)
    df=DataFrame(asv_taxonomy)
    require_columns(df,(:family,:genus,:species,:ASV))
    !isnothing(vsearch_path) && (isnothing(trimmed_folder) || !isdir(trimmed_folder)) && throw(ArgumentError("VSEARCH requires a trimmed_folder directory"))
    ee>=0 || throw(ArgumentError("ee must be nonnegative"))
    thread_budget=threads===:auto ? Threads.nthreads() :
        threads isa Integer && threads>0 ? min(Int(threads),Threads.nthreads()) :
        throw(ArgumentError("threads must be :auto or a positive integer"))
    sort!(df,:species)
    digits=length(string(nrow(df)))
    df.ASV_id=["ASV"*lpad(string(i),digits,'0') for i in 1:nrow(df)]
    mkpath(output_dir)
    write_table(joinpath(output_dir,"asv_taxonomy_revised.xlsx"),df)
    for dir in ("aa_structure_results","query_alignments","vsearch_map_back")
        mkpath(joinpath(output_dir,dir))
    end
    # Remove only reports owned by this pipeline, leaving unrelated files intact.
    resultdir=joinpath(output_dir,"aa_structure_results")
    for file in readdir(resultdir;join=true)
        if isfile(file) && (endswith(file,"_aa_structure.xlsx") || basename(file)=="hydrophobicity.xlsx")
            rm(file)
        end
    end
    samples=[c for c in propertynames(df) if c!=:ASV_id && Base.nonmissingtype(eltype(df[!,c])) <: Number]
    all_hydro=DataFrame(ASV_id=String[],tax_lev=String[],species=Union{Missing,String}[],
        global_min=Union{Missing,Float64}[],global_max=Union{Missing,Float64}[],global_score=Union{Missing,Float64}[],
        global_pass=Union{Missing,Bool}[],local_violation_rate=Union{Missing,Float64}[],structural_verdict=String[])
    # ASVs of a taxon share its references: load, check and warn once per file.
    # A family's references are released after its last ASV.
    references=Dict{String,Dict{String,NamedTuple}}()
    references_lock=ReentrantLock()
    load_references(family,aa_path,nt_paths)=lock(references_lock) do
        get!(()->taxon_references(aa_path,nt_paths),get!(Dict{String,NamedTuple},references,family),aa_path)
    end
    last_row=Dict(safe_name(string(f))=>i for (i,f) in enumerate(df.family) if !ismissing(f))
    function process_asv(row_index)
        row=df[row_index,:]
        ismissing(row.family) && return (hydro=DataFrame(),aa=DataFrame())
        family=safe_name(string(row.family))
        nt_dir=joinpath(alignment_dir,"raw_sequences_nt",family,"species")
        local_hydro=DataFrame()
        all_aa=DataFrame()
        present=string.([s for s in samples if !ismissing(row[s]) && row[s]>0])
        for rank in (SpeciesRank,GenusRank,FamilyRank)
            level=rank_name(rank)
            name=row[rank_symbol(rank)]
            ismissing(name) && continue
            name=safe_name(replace(string(name)," "=>"_"))
            aa_path=joinpath(alignment_dir,"alignments",family,level,name*".fasta")
            nt_paths=taxon_nt_paths(nt_dir,level,name)
            (!isfile(aa_path) || isempty(nt_paths)) && continue
            refs,nt=load_references(family,aa_path,nt_paths)
            isempty(refs) && continue
            result=calculate_asv_distance(;query_dna=row.ASV,query_id=row.ASV_id,ref_aa_aligned=refs,ref_nt=nt,
                output_dir,tax_lev=level,genetic_code,hydro_threshold,hydro_window,vsearch_path,trimmed_folder,
                asv_samples=present,ee,kwargs...)
            append!(local_hydro,add_metadata(result.hydro,row,level);cols=:union)
            append!(all_aa,add_metadata(result.aa,row,level);cols=:union)
        end
        (hydro=local_hydro,aa=all_aa)
    end
    # Bound each batch by the allocated Julia threads so completed per-ASV
    # reports do not accumulate while preserving deterministic output order.
    for batch_start in 1:thread_budget:nrow(df)
        indices=batch_start:min(batch_start+thread_budget-1,nrow(df))
        results = if thread_budget>1 && length(indices)>1
            fetch_all([Threads.@spawn process_asv(i) for i in indices])
        else
            [process_asv(i) for i in indices]
        end
        for (i,result) in zip(indices,results)
            append!(all_hydro,result.hydro;cols=:union)
            isempty(result.aa) || write_table(joinpath(resultdir,df.ASV_id[i]*"_aa_structure.xlsx"),result.aa)
        end
        filter!(((family,_),)->last_row[family]>last(indices),references)
    end
    write_table(joinpath(resultdir,"hydrophobicity.xlsx"),all_hydro)
    df
end
asv_functional_structure(args...;kwargs...)=asv_functional_structure!(args...;kwargs...)

# Wait for every task before rethrowing the first failure (in task order), so no
# task still writes reports after the caller has seen the error.
function fetch_all(tasks)
    failure=nothing
    for t in tasks
        try
            wait(t)
        catch e
            isnothing(failure) && (failure=e)
        end
    end
    isnothing(failure) || throw(failure)
    fetch.(tasks)
end
