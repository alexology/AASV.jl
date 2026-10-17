using Pkg
# Bootstrap the sibling dependency, including on Julia 1.10 where
# Project.toml [sources] is not used by the package resolver.
cd(@__DIR__) do
    Pkg.develop(path="../ProfileAlignments.jl")
    Pkg.instantiate()
end
