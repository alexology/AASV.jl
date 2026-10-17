const AA3 = Dict(zip(collect("ACDEFGHIKLMNPQRSTVWY"),
    ["Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu","Met","Asn","Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr"]))

getvalue(row,key,default=missing) = hasproperty(row,key) ? row[key] : default
finite_values(values) = Float64[x for x in values if !ismissing(x) && isfinite(x)]
mean_or_missing(values) = (v=finite_values(values); isempty(v) ? missing : mean(v))
round_metric(x) = ismissing(x) ? missing : round(x;digits=3)

"""Configuration for the Functionality Index calculation.
Pass this value positionally to `classify_asv` when the same policy is reused.
"""
struct ClassificationConfig
    grantham_max::Float64
    quality_threshold::Float64
    include_position::Bool
    include_grantham::Bool
    include_quality::Bool
    include_hydro::Bool
    include_conservation::Bool
    conservation_k::Float64
    function ClassificationConfig(;grantham_max=215,quality_threshold=0.999,
            include_position=true,include_grantham=true,include_quality=true,
            include_hydro=true,include_conservation=true,conservation_k=20)
        grantham_max>0 && conservation_k>0 ||
            throw(ArgumentError("Normalization constants must be positive"))
        0<=quality_threshold<=1 ||
            throw(ArgumentError("quality_threshold must be in [0,1]"))
        new(Float64(grantham_max),Float64(quality_threshold),Bool(include_position),
            Bool(include_grantham),Bool(include_quality),Bool(include_hydro),
            Bool(include_conservation),Float64(conservation_k))
    end
end

function empty_verdicts()
    DataFrame(ASV_id=String[],tax_lev=String[],FI=Float64[],Class=String[],
        mean_conservation_weight=Union{Missing,Float64}[],mean_conservation_confidence=Union{Missing,Float64}[],
        mean_normalized_grantham=Union{Missing,Float64}[],normalized_hydro=Float64[],mean_site_penalty=Float64[],
        n_pos1=Int[],n_pos2=Int[],codon_penalty=Float64[],total_penalty=Float64[],Evidence=String[],
        n_aa_substitutions=Int[],n_flagged_aa_substitutions=Int[],n_ref_sequences=Union{Missing,Int}[])
end

"""Compute the R package's COI Functionality Index from Excel or TSV reports.
FI = 1 - (0.35 × mean site penalty + 0.30 × codon penalty + 0.35 × hydro).
Hard ORF incompatibilities (`aa_pos=-1`) force FI=0. Low-quality mutations
are excluded; missing quality is accepted. This score measures compatibility
with functional COI evolution and does not establish mitochondrial/NUMT origin.
"""
function classify_asv(output_dir,config::ClassificationConfig)
    (;grantham_max,quality_threshold,include_position,include_grantham,
        include_quality,include_hydro,include_conservation,conservation_k)=config
    dir=joinpath(output_dir,"aa_structure_results")
    files=isdir(dir) ? sort(filter(p->endswith(lowercase(p),".xlsx") || endswith(lowercase(p),".tsv"),readdir(dir;join=true))) : String[]
    hydro_files=filter(p->occursin("hydrophobicity",lowercase(basename(p))),files)
    include_hydro && isempty(hydro_files) && throw(ArgumentError("Hydrophobicity file not found"))
    hydro=isempty(hydro_files) ? DataFrame() : read_table(first(hydro_files))
    isempty(hydro_files) || require_columns(hydro,(:ASV_id,:tax_lev,:local_violation_rate))
    groups=Dict{Tuple{String,String},Vector{NamedTuple}}()
    hydro_values=Dict{Tuple{String,String},Vector{Union{Missing,Float64}}}()
    for row in eachrow(hydro)
        key=(string(row.ASV_id),string(row.tax_lev))
        get!(groups,key,NamedTuple[])
        push!(get!(hydro_values,key,Union{Missing,Float64}[]),row.local_violation_rate)
    end
    for file in filter(p->occursin("_aa_structure",basename(p)),files)
        df=read_table(file)
        :aa_pos in propertynames(df) || continue
        require_columns(df,(:ASV_id,:tax_lev))
        for row in eachrow(df)
            key=(string(row.ASV_id),string(row.tax_lev))
            push!(get!(groups,key,NamedTuple[]),NamedTuple(row))
        end
    end
    out=empty_verdicts()
    for key in sort(collect(keys(groups)))
        rows=groups[key]
        hard=any(r->!ismissing(r.aa_pos) && r.aa_pos==-1,rows)
        reliable=filter(rows) do r
            q=getvalue(r,:median_quality_prob)
            !ismissing(r.aa_pos) && r.aa_pos!=-1 && (!include_quality || ismissing(q) || q>=quality_threshold)
        end
        positions=sort(unique(r.aa_pos for r in reliable))
        sites=NamedTuple[]
        for pos in positions
            site=filter(r->r.aa_pos==pos,reliable)
            aas=[string(r.ref_aa) for r in site if !ismissing(getvalue(r,:ref_aa))]
            counts=[count(==(aa),aas) for aa in unique(aas)]
            probabilities=isempty(aas) ? Float64[] : counts ./ length(aas)
            entropy=-sum(p*log(p) for p in probabilities;init=0.0)/log(20)
            nrefs=finite_values(getvalue(r,:n_sequences_ali) for r in site)
            nref=isempty(nrefs) ? missing : Int(minimum(nrefs))
            confidence=ismissing(nref) ? 1.0 : 1-exp(-nref/conservation_k)
            conservation=include_conservation ? (1-min(entropy,1))*confidence : 1.0
            distances=[getvalue(r,:grantham_dist) for r in site]
            best=argmin([ismissing(x) ? Inf : x for x in distances])
            label=site[best]
            gd=distances[best]
            normalized=include_grantham && !ismissing(gd) ? min(gd/grantham_max,1) : 0.0
            mutation=string(coalesce(getvalue(label,:triplet_mut_pos),""))
            push!(sites,(weight=conservation,confidence=confidence,grantham=normalized,
                penalty=conservation*normalized,pos1=include_position && occursin("1",mutation),
                pos2=include_position && occursin("2",mutation),nref=nref,label=label))
        end
        n=length(sites)
        n1=count(s->s.pos1,sites); n2=count(s->s.pos2,sites)
        codon=n==0 ? 0.0 : (n1+2*n2)/(2*n)
        site_penalty=n==0 ? 0.0 : mean(s.penalty for s in sites)
        hv=mean_or_missing(get(hydro_values,key,Union{Missing,Float64}[]))
        h=include_hydro && !ismissing(hv) ? min(hv,1.0) : 0.0
        penalty=0.35*site_penalty+0.30*codon+0.35*h
        fi=hard ? 0.0 : clamp(1-penalty,0,1)
        class=hard ? "Severe incompatibility with functional COI." : fi>=0.90 ? "Plausible functional sequence" : "Artifact-NUMTs candidate"
        parts=String[]
        if hard
            push!(parts,"query has a premature stop codon, frameshift, or indel relative to the $(key[2])-level reference")
        elseif n==0
            push!(parts,"no reliable novel amino acid substitutions at $(key[2]) level")
        else
            for s in sites
                r=s.label
                from=get(AA3,only(string(getvalue(r,:ref_aa,"X"))),"X")
                to=get(AA3,only(string(getvalue(r,:illegal_aa,"X"))),"X")
                push!(parts,"$from->$to substitution; codon positions $(getvalue(r,:triplet_mut_pos)); Grantham = $(getvalue(r,:grantham_dist)); conservation weight = $(round_metric(s.weight)); confidence = $(round_metric(s.confidence)); site penalty = $(round_metric(s.penalty))")
            end
        end
        push!(parts,"FI = $(round_metric(fi)) ($class)")
        refs=finite_values(s.nref for s in sites)
        push!(out,(key[1],key[2],round_metric(fi),class,
            n==0 ? missing : round_metric(mean(s.weight for s in sites)),
            n==0 ? missing : round_metric(mean(s.confidence for s in sites)),
            n==0 ? missing : round_metric(mean(s.grantham for s in sites)),round_metric(h),round_metric(site_penalty),
            n1,n2,round_metric(codon),round_metric(penalty),join(parts,"; "),n,count(s->s.penalty>0,sites),
            isempty(refs) ? missing : Int(minimum(refs))))
    end
    for row in eachrow(out)
        other=filter(r->r.ASV_id==row.ASV_id && r.tax_lev!=row.tax_lev,eachrow(out))
        isempty(other) || (row.Evidence*="; other levels: "*join(["$(r.tax_lev) FI = $(r.FI)" for r in other],"; "))
    end
    out
end

classify_asv(output_dir;kwargs...)=classify_asv(output_dir,ClassificationConfig(;kwargs...))
