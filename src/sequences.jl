"""The sequence and FASTA API is shared with ProfileAlignments."""
const Sequence = ProfileAlignments.Sequence
const Alignment = ProfileAlignments.Alignment
const ScoringModel = ProfileAlignments.ScoringModel
const read_fasta = ProfileAlignments.read_fasta
const write_fasta = ProfileAlignments.write_fasta

# FASTA entries in file order, without the unique-ID check of `read_fasta`.
function fasta_entries(path::AbstractString)
    entries=Pair{String,String}[]
    reader=FASTX.FASTAReader(open(path);copy=false)
    try
        for rec in reader
            push!(entries,String(strip(FASTX.description(rec)))=>FASTX.sequence(String,rec))
        end
    catch e
        e isa ErrorException || e isa EOFError || rethrow()
        throw(ArgumentError("Invalid FASTA $path: $(sprint(showerror,e))"))
    finally
        close(reader)
    end
    entries
end

"""Read reference FASTA records, dropping exact repeats of a record (same ID and
sequence) and returning `(records, repeats)`. Reference sets built by earlier aasv
versions can repeat whole records; a repeated ID with a different sequence is
ambiguous and raises `ArgumentError`.
"""
function read_references(path::AbstractString)
    records=Sequence[]
    seen=Dict{String,String}()
    repeats=0
    for (id,sequence) in fasta_entries(path)
        previous=get(seen,id,nothing)
        if isnothing(previous)
            seen[id]=sequence
            push!(records,Sequence(id,sequence))
        elseif previous==sequence
            repeats+=1
        else
            throw(ArgumentError("Reference ID $id is repeated with a different sequence in $path"))
        end
    end
    (records=records,repeats=repeats)
end

"""Reference records of one file; exact repeats are ignored with a warning."""
function read_reference_records(path::AbstractString)
    loaded=read_references(path)
    loaded.repeats>0 && @warn "Exact repeats of reference records were ignored" path repeats=loaded.repeats
    loaded.records
end

"""A stop-free forward translation and its one-based reading frame."""
struct Translation
    amino_acids::String
    frame::Int
end

const CODONS = [string(a,b,c) for a in "TCAG" for b in "TCAG" for c in "TCAG"]
# Tables are baked in at precompilation; recompile whenever the data file changes.
const GENETIC_CODES_PATH = joinpath(@__DIR__,"..","data","genetic_codes.tsv")
include_dependency(GENETIC_CODES_PATH)
const GENETIC_TABLES = Dict(begin
    id,aas=split(line,'\t')
    String(id)=>Dict(zip(CODONS,collect(aas)))
end for line in readlines(GENETIC_CODES_PATH)[2:end])
const IUPAC = Dict('A'=>"A", 'C'=>"C", 'G'=>"G", 'T'=>"T", 'N'=>"ACGT",
    'R'=>"AG", 'Y'=>"CT", 'S'=>"CG", 'W'=>"AT", 'K'=>"GT", 'M'=>"AC",
    'B'=>"CGT", 'D'=>"AGT", 'H'=>"ACT", 'V'=>"ACG")

function genetic_table(code)
    if code isa AbstractDict
        table=Dict(String(k)=>only(string(v)) for (k,v) in code)
        all(c->haskey(table,c),CODONS) || throw(ArgumentError("Custom genetic code requires all 64 codons"))
        return table
    end
    gc=string(code)
    haskey(GENETIC_TABLES,gc) || throw(ArgumentError("Unsupported genetic code $gc; supply a complete codon Dict"))
    GENETIC_TABLES[gc]
end

function translate_codon(codon, table)
    haskey(table, codon) && return table[codon]
    bases = [get(IUPAC, c, "") for c in codon]
    any(isempty, bases) && throw(ArgumentError("Invalid DNA codon: $codon"))
    aas = Set(table[string(a,b,c)] for a in bases[1] for b in bases[2] for c in bases[3])
    length(aas) == 1 ? only(aas) : 'X'
end

"""Return the first stop-free forward translation, or `nothing` when none exists."""
function find_translation(dna::AbstractString; genetic_code="5")
    s = uppercase(String(dna))
    all(c -> haskey(IUPAC,c), s) || throw(ArgumentError("DNA contains invalid characters"))
    table = genetic_table(genetic_code)
    for frame in 1:min(3,length(s))
        aa = String([translate_codon(s[i:i+2], table) for i in frame:3:length(s)-2])
        if !occursin('*', aa)
            return Translation(aa,frame)
        end
    end
    nothing
end
find_translation(dna::BioSequences.NucSeq;kwargs...)=find_translation(String(dna);kwargs...)

"""Compatibility adapter returning the sequence or frame, with `missing` on failure.
New code should use `find_translation`, whose return type does not depend on an option.
"""
function best_translation(dna;genetic_code="5",res_type=:sequence)
    kind=res_type isa Symbol ? res_type : Symbol(res_type)
    kind in (:sequence,:position) || throw(ArgumentError("Invalid res_type"))
    result=find_translation(dna;genetic_code)
    isnothing(result) && return missing
    kind===:sequence ? result.amino_acids : result.frame
end

function aa_table(records::AbstractVector{Sequence}; genetic_code="5")
    seq_ids=String[]; positions=Int[]; triplets=String[]; amino_acids=String[]
    table = genetic_table(genetic_code)
    for r in records
        translation=find_translation(r.sequence;genetic_code)
        isnothing(translation) && continue
        for (position,i) in enumerate(translation.frame:3:length(r.sequence)-2)
            triplet = r.sequence[i:i+2]
            push!(seq_ids,r.id); push!(positions,position); push!(triplets,triplet)
            push!(amino_acids,string(translate_codon(triplet,table)))
        end
    end
    DataFrame(seq_id=seq_ids,position=positions,triplet=triplets,aa=amino_acids)
end

"""Dense 64x64 codon Hamming-similarity matrix in `CODONS` order."""
function codon_similarity_matrix(;genetic_code="5")
    genetic_table(genetic_code)
    [count(k->a[k]==b[k],1:3)/3 for a in CODONS,b in CODONS]
end

"""Tabular compatibility view of `codon_similarity_matrix`."""
function aa_similarity_matrix(;genetic_code="5")
    values=codon_similarity_matrix(;genetic_code)
    DataFrame([(triplets_1=CODONS[i],triplets_2=CODONS[j],value=values[i,j])
        for i in eachindex(CODONS) for j in eachindex(CODONS)])
end

fasta_files(dir) = isdir(dir) ? sort(filter(p->endswith(lowercase(p),".fasta"), readdir(dir; join=true))) : String[]
safe_name(s) = (occursin(r"[\\/:*?\"<>|]",s) || s in (".","..")) ?
    throw(ArgumentError("Unsafe taxon or file name: $s")) : String(s)
function write_table(path, df)
    mkpath(dirname(abspath(path)))
    if endswith(lowercase(path), ".xlsx")
        XLSX.writetable(path, df; overwrite=true)
    else
        CSV.write(path, df; delim=endswith(path,".csv") ? ',' : '\t')
    end
    path
end
read_table(path) = endswith(lowercase(path),".xlsx") ?
    DataFrame(XLSX.readtable(path,1; infer_eltypes=true)) :
    CSV.read(path,DataFrame; delim=endswith(path,".csv") ? ',' : '\t')
function require_columns(df, cols)
    all(c->c in propertynames(df), cols) || throw(ArgumentError("Required columns: $(join(cols, ", "))"))
end
