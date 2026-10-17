using AASV, Test, DataFrames, CSV
import ProfileAlignments

# julia --project=. test/parity.jl "path/to/Rscript.exe" [path/to/aasv/R/repository]
isempty(ARGS) && error("Supply the Rscript executable path")
rscript = ARGS[1]
source = length(ARGS)>1 ? ARGS[2] : joinpath(@__DIR__,"..","reference","aasv")
isfile(joinpath(source,"R","classify_asv.R")) || error("Original R sources not found")

mktempdir() do root
    cases = [(name="default", pos=5, gr=20, mutation="3", q=0.9999, nref=missing, h=0.0),
        (name="position1", pos=5, gr=20, mutation="1", q=0.9999, nref=missing, h=0.0),
        (name="position2", pos=5, gr=20, mutation="2", q=0.9999, nref=missing, h=0.2),
        (name="lowquality", pos=5, gr=20, mutation="3", q=0.8, nref=missing, h=0.0),
        (name="sparse", pos=5, gr=20, mutation="3", q=0.9999, nref=2, h=0.0),
        (name="manyrefs", pos=5, gr=215, mutation="1,2", q=missing, nref=200, h=0.3),
        (name="hard", pos=-1, gr=20, mutation="3", q=0.9999, nref=missing, h=0.0)]
    for c in cases
        dir=joinpath(root,c.name,"aa_structure_results")
        mkpath(dir)
        AASV.write_table(joinpath(dir,"hydrophobicity.xlsx"),DataFrame(ASV_id=["ASV1","ASV2"],
            tax_lev=["species","species"],local_violation_rate=[c.h,0.0],global_pass=[true,true]))
        AASV.write_table(joinpath(dir,"ASV1_aa_structure.xlsx"),DataFrame(ASV_id=["ASV1"],tax_lev=["species"],
            aa_pos=[c.pos],illegal_aa=["F"],ref_aa=["L"],grantham_dist=[c.gr],triplet_mut_pos=[c.mutation],
            median_quality_prob=[c.q],n_sequences_ali=[c.nref]))
    end
    # Variable reference position: conservation uses every reference, while
    # the substitution label uses the lowest Grantham distance.
    dir=joinpath(root,"variable","aa_structure_results")
    mkpath(dir)
    AASV.write_table(joinpath(dir,"hydrophobicity.xlsx"),DataFrame(ASV_id=["ASV1"],tax_lev=["species"],local_violation_rate=[0.1],global_pass=[true]))
    AASV.write_table(joinpath(dir,"ASV1_aa_structure.xlsx"),DataFrame(ASV_id=fill("ASV1",4),tax_lev=fill("species",4),
        aa_pos=fill(5,4),illegal_aa=fill("F",4),ref_aa=["L","L","I","V"],grantham_dist=[22,22,21,50],
        triplet_mut_pos=fill("2",4),median_quality_prob=fill(0.9999,4),n_sequences_ali=fill(40,4)))
    script=joinpath(@__DIR__,"parity.R")
    run(`$rscript $script $source $root`)
    @testset "Numerical parity with original R classification" begin
        for folder in filter(isdir,readdir(root;join=true))
            actual=sort(classify_asv(folder),[:ASV_id,:tax_lev])
            expected=sort(CSV.read(joinpath(folder,"r_verdicts.tsv"),DataFrame;delim='\t'),[:ASV_id,:tax_lev])
            @test names(actual)==names(expected)
            for col in propertynames(actual)
                col==:Evidence && continue
                @test isequal(actual[!,col],expected[!,col])
            end
        end
        translated=CSV.read(joinpath(root,"r_translation.tsv"),DataFrame;delim='\t')
        for row in eachrow(translated)
            @test isequal(best_translation(row.dna),row.aa)
        end
        hydro=CSV.read(joinpath(root,"r_hydro.tsv"),DataFrame;delim='\t')
        for row in eachrow(hydro)
            @test get_hydro_metrics(row.aa;hydro_window=3).global_score≈row.hydro
        end
        codes=CSV.read(joinpath(root,"r_codes.tsv"),DataFrame;delim='\t')
        for row in eachrow(codes)
            table=AASV.genetic_table(row.id)
            @test String([table[codon] for codon in AASV.CODONS])==row.AAs
        end
        blosum=CSV.read(joinpath(root,"r_blosum.tsv"),DataFrame;delim='\t')
        indices=[findfirst(==(string(c)),names(blosum)) for c in ProfileAlignments.PROTEIN_ALPHABET]
        @test Matrix(blosum[indices,indices])==ProfileAlignments.BLOSUM62
    end
end
