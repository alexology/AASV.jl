using Pkg
# Bootstrap the sibling dependency in the local Manifest. Pkg >= 1.11 also
# records `develop(path=...)` in Project.toml [sources]; that machine-specific
# path breaks installation elsewhere, so the project file is restored.
cd(@__DIR__) do
    project = read("Project.toml")
    try
        Pkg.develop(path="../ProfileAlignments.jl")
    finally
        write("Project.toml", project)
    end
    Pkg.instantiate()
end
