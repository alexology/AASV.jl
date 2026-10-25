module AASV

using CSV, DataFrames, HTTP, JSON3, SHA, Statistics, XLSX
import BioSequences, FASTX
using PrecompileTools: @setup_workload, @compile_workload
import ProfileAlignments

export Sequence, Alignment, ScoringModel, Translation, read_fasta, write_fasta,
    find_translation, best_translation, aa_table, aa_similarity_matrix,
    codon_similarity_matrix,
    TaxonomicRank, SpeciesRank, GenusRank, FamilyRank,
    HydrophobicityMetrics, StructureAssessment, StructuralVerdict,
    StructureUnavailable, StructurePass, StructureFail, InvalidORF,
    MutationObservation, ASVAssessment, ASVStatus, ASVValid, ASVInvalidORF,
    assess_structure, assess_asv,
    align_species_seq, align_genus_seq, align_family_seq,
    align_species_seq!, align_genus_seq!, align_family_seq!, check_alignment,
    improve_alignment, improve_alignment!, filter_divergent_seq, select_shard,
    check_alignments!, filter_divergent_references!, check_reference_records,
    filter_divergent_seq!, re_align, re_align!, asv_functional_structure,
    asv_functional_structure!, ClassificationConfig, classify_asv,
    align_sequences, align_profiles,
    calculate_asv_distance, get_hydro_metrics, check_asv_structure

include("domain.jl")
include("sequences.jl")
include("alignment.jl")
include("references.jl")
include("reference_checks.jl")
include("structure.jl")
include("classification.jl")
include("precompile.jl")

end
