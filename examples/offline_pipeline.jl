using AASV, DataFrames, CSV

# Small synthetic reference set: runs without network or VSEARCH.
reference_dna = "ATGGCTGCTGCTGCTGCTGCTGCTGCT"
taxonomy = DataFrame(family=["DemoFamily"], genus=["Demo"], species=["Demo species"])
records = [(processid="REF1", identification="Demo species", sequence=reference_dna),
           (processid="REF2", identification="Demo second", sequence=reference_dna)]
base = joinpath(@__DIR__, "..", "results", "demo")
references = joinpath(base, "references")
output = joinpath(base, "analysis")

align_species_seq!(taxonomy; alignment_dir=references, min_length=0,
                   overwrite=true, records_provider=family -> records)
problematic = check_alignments!(; alignment_dir=references)
isempty(problematic) || improve_alignment!(problematic; log_dir=joinpath(references,"logs"))
filter_divergent_references!(; alignment_dir=references)
align_genus_seq!(; alignment_dir=references, overwrite=true)
align_family_seq!(; alignment_dir=references, overwrite=true)

asvs = DataFrame(family=fill("DemoFamily", 3), genus=fill("Demo", 3),
    species=fill("Demo species", 3),
    ASV=[reference_dna, "ATGGCTGCTGATGCTGCTGCTGCTGCT", "TAAATAAATAA"],
    sample_1=[100, 10, 1])
asv_functional_structure!(asvs; output_dir=output, alignment_dir=references)
verdicts = classify_asv(output)
CSV.write(joinpath(output, "classification.csv"), verdicts)
show(stdout, verdicts[:, [:ASV_id, :tax_lev, :FI, :Class]]; allrows=true)
println("\nResults: ", abspath(output))
