"""Taxonomic rank used by reference and ASV workflows."""
@enum TaxonomicRank begin
    SpeciesRank
    GenusRank
    FamilyRank
end

rank_symbol(rank::TaxonomicRank)=rank===SpeciesRank ? :species :
    rank===GenusRank ? :genus : :family
rank_name(rank::TaxonomicRank)=String(rank_symbol(rank))

parse_rank(rank::TaxonomicRank)=rank
function parse_rank(rank::Union{Symbol,AbstractString})
    value=Symbol(lowercase(String(rank)))
    value===:species && return SpeciesRank
    value===:genus && return GenusRank
    value===:family && return FamilyRank
    throw(ArgumentError("Invalid taxonomic rank: $rank"))
end
