using Test, AASV, BioSequences, DataFrames, CSV, XLSX, HTTP, JSON3, Statistics
import ProfileAlignments
include("reference_checks.jl")
include("divergent_distances.jl")

@testset "Project file is portable" begin
    # Pkg >= 1.11 copies the local ProfileAlignments path from the Manifest into
    # Project.toml on add/rm/update; that path breaks installation elsewhere.
    project=read(joinpath(pkgdir(AASV),"Project.toml"),String)
    @test !occursin(r"(?m)^\[sources\]",project)
end

@testset "FASTA and translation" begin
    @test Sequence === ProfileAlignments.Sequence
    mktempdir() do dir
        records=[Sequence("one full header","atgccagca"),Sequence("two","ATGAAATTT")]
        path=write_fasta(records,joinpath(dir,"sequences.fasta");width=4)
        reread=read_fasta(path)
        @test getfield.(reread,:id)==getfield.(records,:id)
        @test getfield.(reread,:sequence)==getfield.(records,:sequence)
        @test_throws ArgumentError write_fasta([records[1],records[1]],path)
        write(path,"ATGC\n")
        @test_throws ArgumentError read_fasta(path)
        @test nrow(aa_table(records))==6
        @test aa_table(records).triplet[1]=="ATG"
    end
    @test best_translation("ATGCCAGCA")=="MPA"
    translation=find_translation("ATGCCAGCA")
    @test translation isa Translation
    @test translation.amino_acids=="MPA"
    @test translation.frame==1
    @test isnothing(find_translation("TAAATAAATAA"))
    @test best_translation("TAACCAGCA";res_type="position")==2
    @test ismissing(best_translation(""))
    @test ismissing(best_translation("TAAATAAATAA"))
    @test best_translation("ATATGAAGAAGG")=="MWSS"
    @test best_translation("ATATGAAGAAGG";genetic_code="2",res_type="position")!=1
    @test best_translation("GCN")=="A"
    @test best_translation("NNN")=="X"
    @test best_translation(dna"ATGCCAGCA")=="MPA"
    @test Sequence("bio-dna",dna"ACGN-").sequence=="ACGN-"
    @test Sequence("bio-aa",aa"ACDX-").sequence=="ACDX-"
    @test get_hydro_metrics(aa"ACDE").global_score==get_hydro_metrics("ACDE").global_score
    @test_throws ArgumentError best_translation("AB!")
    @test_throws ArgumentError best_translation("ATG";genetic_code="99")
    sim=aa_similarity_matrix()
    dense=codon_similarity_matrix()
    @test size(dense)==(64,64)
    @test dense[findfirst(==("AAA"),AASV.CODONS),findfirst(==("AAT"),AASV.CODONS)]≈2/3
    @test nrow(sim)==4096
    @test only(sim.value[(sim.triplets_1.=="AAA") .& (sim.triplets_2.=="AAA")])==1
    @test only(sim.value[(sim.triplets_1.=="AAA") .& (sim.triplets_2.=="AAT")])≈2/3
end

matrix(seqs)=AASV.alignment_matrix([Sequence(string(i),s) for (i,s) in enumerate(seqs)])
@testset "Profiles and alignment cleanup" begin
    @test AASV.inner_gaps(matrix(["--ABC-","A-CDE-","------"]))==[false,true,false]
    @test AASV.causes_internal_gaps(matrix(["ACXEF","AC-EF","AC-EF","AC-EF"]))==[true,false,false,false]
    @test !any(AASV.causes_internal_gaps(matrix(["ACDEF","--DEF","ACDE-"])))
    records=[Sequence("a","ACDEFG"),Sequence("b","CDEF"),Sequence("dup","ACDEFG")]
    aligned=align_sequences(records)
    @test getfield.(aligned,:id)==["a","b","dup"]
    @test replace.(getfield.(aligned,:sequence),"-"=>"")==getfield.(records,:sequence)
    @test !any(AASV.inner_gaps(AASV.alignment_matrix(aligned)))
    @test aligned[1].sequence==aligned[3].sequence
    merged=align_profiles([Sequence("a","AC-EF"),Sequence("b","AC-EF")],[Sequence("c","ACXEF")])
    @test merged[1].sequence==merged[2].sequence
    @test replace(merged[3].sequence,"-"=>"")=="ACXEF"
    mktempdir() do dir
        path=joinpath(dir,"alignments","Family","species","Genus_species.fasta")
        write_fasta([Sequence("Genus_species__a","ACDEFGXHIKLMNPQ"),Sequence("Genus_species__b","ACDEFG-HIKLMNPQ"),
            Sequence("Genus_species__c","ACDEFG-HIKLMNPR"),Sequence("Genus_species__d","ACDEFG-HIKLMNPS")],path)
        @test check_alignment("Family";alignment_dir=dir,tax_lev=SpeciesRank)==[path]
        improve_alignment!([path];print_messages=false,log_dir=joinpath(dir,"logs"))
        @test getfield.(read_fasta(path),:id)==["Genus_species__b","Genus_species__c","Genus_species__d"]
        @test isempty(check_alignment("Family";alignment_dir=dir))
        @test isfile(joinpath(dir,"logs","removed_sequences.tsv"))
        re_align!([path];verbose=false)
        write_fasta([Sequence("a","AAAA"),Sequence("b","AAAA"),Sequence("c","AAAT"),Sequence("d","WWWW")],path)
        filter_divergent_seq!([path];max_distance=0.3,print_messages=false,log_dir=joinpath(dir,"logs"))
        @test length(read_fasta(path))==3
        @test isfile(joinpath(dir,"logs","removed_divergent.tsv"))
        write_fasta([Sequence("a","AAAA"),Sequence("b","WWWW")],path)
        @test_logs (:warn,r"Every sequence") filter_divergent_seq([path];min_seqs=2,print_messages=false)
        @test length(read_fasta(path))==2
        @test_logs (:info,r"started") (:info,r"filter: file") (:info,r"completed") filter_divergent_seq!([path];print_messages=false,show_progress=true)
        @test_logs filter_divergent_seq!([path];print_messages=false)
        @test_throws ArgumentError filter_divergent_seq!([path];progress_interval=0)
        write_fasta([Sequence("seq_$i","ACDE") for i in 1:1000],path)
        before = read(path,String)
        @test_logs match_mode=:any (:info,r"filter: comparisons") filter_divergent_seq!([path];print_messages=false,show_progress=true,progress_interval=eps())
        @test read(path,String)==before
    end
end

@testset "Hydrophobicity, coordinates and quality" begin
    @test AASV.grantham_distance('S','F')==155
    @test AASV.grantham_distance('C','W')==215
    @test AASV.grantham_distance('L','I')==5
    @test ismissing(AASV.grantham_distance('X','I'))
    h=get_hydro_metrics("AAAA";hydro_window=2)
    @test h isa HydrophobicityMetrics
    @test h.global_score≈1.8
    @test h.profile≈[1.8,1.8,1.8]
    @test isempty(get_hydro_metrics("AA";hydro_window=7).profile)
    @test only(check_asv_structure("AAAA",(1.0,2.0),fill(1.0,3),fill(2.0,3);hydro_window=2).structural_verdict)=="PASS"
    typed_structure=assess_structure("AAAA",(1.0,2.0),fill(1.0,3),fill(2.0,3);hydro_window=2)
    @test typed_structure isa StructureAssessment
    @test typed_structure.verdict==AASV.StructurePass
    nt=[Sequence("Genus_species__one","ATGGCTGCTGCTGCTGCTGCTGCTGCT")]
    aa=[Sequence(nt[1].id,best_translation(nt[1].sequence))]
    clean=calculate_asv_distance(;query_dna=nt[1].sequence,query_id="clean",ref_aa_aligned=aa,ref_nt=nt)
    @test isempty(clean.aa)
    @test only(clean.hydro.local_violation_rate)==0
    query="ATGGCTGCTGATGCTGCTGCTGCTGCT"
    result=calculate_asv_distance(;query_dna=query,query_id="mut",ref_aa_aligned=aa,ref_nt=nt)
    typed=assess_asv(;query_dna=query,query_id="mut",ref_aa_aligned=aa,ref_nt=nt)
    @test typed isa ASVAssessment
    @test typed.status==AASV.ASVValid
    @test only(typed.mutations) isa MutationObservation
    @test only(typed.mutations).query_nt==10:12
    @test result.aa.aa_pos==[4]
    @test result.aa.dna_triplet==["GAT"]
    @test result.aa.dna_triplet_ref==["GCT"]
    @test result.aa.nt_coords==["10-12"]
    @test result.aa.triplet_mut_pos==["2"]
    @test result.aa.grantham_dist==[126]
    hard=calculate_asv_distance(;query_dna="TAAATAAATAA",query_id="bad",ref_aa_aligned=aa,ref_nt=nt)
    @test only(hard.aa.aa_pos)==-1
    mktempdir() do dir
        path=joinpath(dir,"test.fastq")
        write(path,"@a\nATG\n+\nIII\n@b\nATG\n+\n!!!\n")
        q=AASV.quality_probabilities(path,1)
        @test length(q)==1
        @test sum(q[1])≈0.0003
    end
end

function classification_fixture(dir; pos=5, gr=20, mutation="3", quality=0.9999, nref=missing, hydro=0.0)
    resultdir=joinpath(dir,"aa_structure_results")
    mkpath(resultdir)
    AASV.write_table(joinpath(resultdir,"hydrophobicity.xlsx"),DataFrame(ASV_id=["ASV1","ASV2"],tax_lev=["species","species"],local_violation_rate=[hydro,0.0]))
    AASV.write_table(joinpath(resultdir,"ASV1_aa_structure.xlsx"),DataFrame(ASV_id=["ASV1"],tax_lev=["species"],aa_pos=[pos],
        illegal_aa=["F"],ref_aa=["L"],grantham_dist=[gr],triplet_mut_pos=[mutation],median_quality_prob=[quality],n_sequences_ali=[nref]))
end

@testset "Functionality Index" begin
    mktempdir() do dir
        @test_throws ArgumentError classify_asv(dir)
        @test isempty(classify_asv(dir;include_hydro=false))
        classification_fixture(dir)
        v=classify_asv(dir)
        @test classify_asv(dir,ClassificationConfig()).FI==v.FI
        @test v.FI≈[0.967,1.0]
        @test v.n_aa_substitutions==[1,0]
        @test ismissing(v.mean_conservation_weight[2])
        classification_fixture(dir;mutation="1")
        @test classify_asv(dir).FI[1]==0.817
        classification_fixture(dir;mutation="2",hydro=0.2)
        v=classify_asv(dir)
        @test v.FI[1]≈round(1-0.35*20/215-0.3-0.35*0.2;digits=3)
        @test v.codon_penalty[1]==1
        classification_fixture(dir;quality=0.8)
        @test classify_asv(dir).FI[1]==1
        @test classify_asv(dir;include_quality=false).FI[1]==0.967
        classification_fixture(dir;nref=2)
        @test classify_asv(dir).mean_conservation_confidence[1]==round(1-exp(-2/20);digits=3)
        @test classify_asv(dir;include_conservation=false).FI[1]==0.967
        classification_fixture(dir;pos=-1)
        @test classify_asv(dir).FI[1]==0
        @test classify_asv(dir).Class[1]=="Severe incompatibility with functional COI."
        classification_fixture(dir)
        @test classify_asv(dir;include_grantham=false).FI[1]==1
    end
end

const DEMO_DNA="ATGGCTGCTGCTGCTGCTGCTGCTGCT"
@testset "Species completion marker after failed overwrite" begin
    mktempdir() do dir
        taxonomy=DataFrame(family=["Family"],genus=["Genus"],species=["Genus species"])
        records=[(processid="one",identification="Genus species",sequence=DEMO_DNA)]
        marker=joinpath(dir,"raw_sequences","Family",".complete")
        mkpath(dirname(marker))
        write(marker,"completed\n")
        interrupted=family->begin
            @test !isfile(marker)
            throw(InterruptException())
        end
        @test_throws InterruptException align_species_seq!(taxonomy;alignment_dir=dir,
            overwrite=true,records_provider=interrupted,show_progress=false)
        @test !isfile(marker)
        calls=Ref(0)
        provider=family->begin
            calls[]+=1
            records
        end
        align_species_seq!(taxonomy;alignment_dir=dir,min_length=0,
            records_provider=provider,show_progress=false)
        @test calls[]==1
        @test isfile(marker)
        align_species_seq!(taxonomy;alignment_dir=dir,records_provider=provider,show_progress=false)
        @test calls[]==1
    end
end

@testset "Offline full pipeline and reports" begin
    mktempdir() do dir
        taxonomy=DataFrame(family=["Family"],genus=["Genus"],species=["Genus species"])
        records=[(processid="one",identification="Genus species",sequence=DEMO_DNA),
            (processid="two",identification="Genus other",sequence=DEMO_DNA),
            (processid="three",identification="Genus sp.",sequence=DEMO_DNA),
            (processid="ambiguous",identification="Genus species",sequence="ATGNGCTGCT"),
            (processid="gap",identification="Genus species",sequence="ATG-GCTGCT")]
        provider=family->records
        @test_logs (:info,r"alignment: started") (:info,r"loading records") (:info,r"filtering and translating") (:info,r"filtering completed") (:info,r"aligning species") (:info,r"aligning species") (:info,r"family completed") (:info,r"alignment: completed") align_species_seq(taxonomy;alignment_dir=dir,min_length=0,records_provider=provider)
        no_download=family->error("Completed families must not be downloaded again")
        @test_logs (:info,r"alignment: started") (:info,r"reusing completed family") (:info,r"alignment: completed") align_species_seq!(taxonomy;alignment_dir=dir,records_provider=no_download)
        @test_logs align_species_seq!(taxonomy;alignment_dir=dir,records_provider=no_download,show_progress=false)
        @test_throws ArgumentError align_species_seq!(taxonomy;alignment_dir=dir,progress_interval=0)
        @test length(AASV.fasta_files(joinpath(dir,"alignments","Family","species")))==2
        @test length(read_fasta(joinpath(dir,"raw_sequences_nt","Family","species","Genus_species.fasta")))==1
        @test_logs (:info,r"alignment: started") (:info,r"loading family") (:info,r"loading taxon") (:info,r"profiles: merging") (:info,r"merge completed") (:info,r"taxon completed") (:info,r"family completed") (:info,r"alignment: completed") align_genus_seq(;alignment_dir=dir)
        @test_logs (:info,r"alignment: started") (:info,r"loading family") (:info,r"loading taxon") (:info,r"taxon completed") (:info,r"family completed") (:info,r"alignment: completed") align_family_seq(;alignment_dir=dir)
        @test_logs (:info,r"alignment: started") (:info,r"loading family") (:info,r"reusing completed taxon") (:info,r"family completed") (:info,r"alignment: completed") align_genus_seq!(;alignment_dir=dir)
        @test_logs align_genus_seq!(;alignment_dir=dir,show_progress=false)
        @test_logs align_family_seq!(;alignment_dir=dir,show_progress=false)
        @test length(read_fasta(joinpath(dir,"alignments","Family","genus","Genus.fasta")))==2
        @test isfile(joinpath(dir,"aa_tables","Family","Genus_species.txt"))
        input=DataFrame(family=fill("Family",3),genus=fill("Genus",3),species=fill("Genus species",3),
            ASV=[DEMO_DNA,"ATGGCTGCTGATGCTGCTGCTGCTGCT","TAAATAAATAA"],sample1=[10,5,1])
        output=joinpath(dir,"results")
        revised=asv_functional_structure(input;output_dir=output,alignment_dir=dir,threads=2)
        @test revised.ASV_id==["ASV1","ASV2","ASV3"]
        v=classify_asv(output)
        @test nrow(v)==9
        # Each row lists the other levels of the same ASV in table order.
        for row in eachrow(v)
            others=[r for r in eachrow(v) if r.ASV_id==row.ASV_id && r.tax_lev!=row.tax_lev]
            @test endswith(row.Evidence,"; other levels: "*
                join(["$(r.tax_lev) FI = $(r.FI)" for r in others],"; "))
        end
        @test all(==(1),v.FI[v.ASV_id.=="ASV1"])
        @test all(==(0),v.FI[v.ASV_id.=="ASV3"])
        @test all(<(1),v.FI[v.ASV_id.=="ASV2"])
        @test isfile(joinpath(output,"query_alignments","ASV2_species.fasta"))
        # Repeated runs remove stale mutation reports.
        asv_functional_structure(input[1:1,:];output_dir=output,alignment_dir=dir)
        @test nrow(classify_asv(output))==3
        @test !isfile(joinpath(output,"aa_structure_results","ASV2_aa_structure.xlsx"))
        @test_throws ArgumentError asv_functional_structure(input;output_dir=output,
            alignment_dir=dir,threads=0)
    end
end

@testset "Family merge progress" begin
    mktempdir() do dir
        for (genus,sequence) in [("Alpha","ACDEFG"),("Beta","ACDEFG"),("Gamma","ACDEYG")]
            write_fasta([Sequence(genus,sequence)],joinpath(dir,"alignments","Family","genus",genus*".fasta"))
        end
        @test_logs (:info,r"alignment: started") (:info,r"loading family") (:info,r"loading taxon") (:info,r"profiles: merging") (:info,r"merge completed") (:info,r"profiles: merging") (:info,r"merge completed") (:info,r"taxon completed") (:info,r"family completed") (:info,r"alignment: completed") align_family_seq!(;alignment_dir=dir)
        path=joinpath(dir,"alignments","Family","family","Family.fasta")
        before=read(path,String)
        @test_logs align_family_seq!(;alignment_dir=dir,overwrite=true,show_progress=false)
        @test read(path,String)==before
        @test length(read_fasta(path))==3
    end
end

@testset "BOLD persistent download recovery" begin
    mktempdir() do dir
        urls=String[]
        interrupted=Ref(true)
        generation=Ref("original")
        requester=function(url;timeout=0)
            push!(urls,url)
            if occursin("/query?",url)
                return (query_id=generation[],)
            elseif occursin("length=0",url)
                return (recordsTotal=3,)
            end
            offset=parse(Int,last(split(url,"start=")))
            interrupted[] && offset==1 && throw(InterruptException())
            # Non-COI records still count toward the next page offset.
            marker=offset==1 ? "other" : "COI-5P"
            (data=[(processid="$(generation[])_$offset",identification="Genus species",
                    marker_code=marker,nuc="ATG")],)
        end
        fetch(;kwargs...)=AASV.fetch_bold_records("Family";cache_dir=dir,page_size=1,
            requester,sleeper=_->nothing,kwargs...)
        @test_throws InterruptException fetch()
        state=JSON3.read(read(joinpath(dir,"download.json"),String))
        pages=joinpath(dir,state.pages)
        @test isfile(joinpath(pages,"0.json"))
        @test !isfile(joinpath(pages,"1.json"))
        # A temporary file left by a terminated process is not a checkpoint.
        write(joinpath(pages,"unfinished.tmp"),"{\"data\":[")
        empty!(urls)
        interrupted[]=false
        records=fetch()
        @test getfield.(records,:processid)==["original_0","original_2"]
        @test length(urls)==2
        @test endswith(urls[1],"start=1")
        @test endswith(urls[2],"start=2")
        @test all(url->occursin("/documents/original?",url),urls)
        empty!(urls)
        @test fetch()==records
        @test isempty(urls) # Fully cached downloads work without network access.
        @test_throws ArgumentError AASV.fetch_bold_records("Different family";
            cache_dir=dir,page_size=1,requester)
        @test_throws ArgumentError AASV.fetch_bold_records("Family";
            cache_dir=dir,page_size=2,requester)
        generation[]="fresh"
        fresh=fetch(overwrite=true)
        @test getfield.(fresh,:processid)==["fresh_0","fresh_2"]
        @test length(urls)==5 # New query, count, and all three pages.
        @test JSON3.read(read(joinpath(dir,"download.json"),String)).pages!=state.pages
        failed=(url;timeout=0)->throw(InterruptException())
        @test_throws InterruptException AASV.fetch_bold_records("Family";
            cache_dir=dir,page_size=1,requester=failed,overwrite=true)
        @test !isfile(joinpath(dir,"download.json"))
    end
    mktempdir() do dir
        requester=(url;timeout=0)->occursin("/query?",url) ?
            (query_id="empty",) : occursin("length=0",url) ? (recordsTotal=1,) : (data=[],)
        @test_throws ErrorException AASV.fetch_bold_records("Family";cache_dir=dir,requester)
        state=JSON3.read(read(joinpath(dir,"download.json"),String))
        @test isempty(readdir(joinpath(dir,state.pages)))
    end
end

@testset "BOLD pagination and retry" begin
    urls=String[]
    requester=function(url;timeout=0)
        push!(urls,url)
        if occursin("/query?",url)
            JSON3.read("{\"query_id\":\"demo\"}")
        elseif occursin("length=0",url)
            JSON3.read("{\"recordsTotal\":2}")
        else
            JSON3.read("{\"data\":[{\"processid\":\"x\",\"identification\":\"Genus species\",\"marker_code\":\"COI-5P\",\"nuc\":\"ATG\"}]}")
        end
    end
    records=AASV.fetch_bold_records("Test family";page_size=1,requester,sleeper=_->nothing)
    @test length(records)==2
    @test endswith(last(urls),"start=1")
    @test occursin("Test%20family",first(urls))
    @test_logs (:info,r"Downloading BOLD") (:info,r"BOLD download progress") (:info,r"Downloading BOLD") (:info,r"BOLD download progress") AASV.fetch_bold_records("Test family";page_size=1,requester,sleeper=_->nothing,show_progress=true)
    statuses=[503,503,200,503,200]
    pauses=Float64[]
    transport=(url;kwargs...)->HTTP.Response(popfirst!(statuses),"{}")
    perform=AASV.make_retry_performer(;initial_wait=1,max_wait=10,max_attempts=5,transport,sleeper=x->push!(pauses,x))
    perform("test"); perform("test")
    @test pauses==[1,2,4]
    bad=AASV.make_retry_performer(;transport=(url;kwargs...)->HTTP.Response(404),sleeper=_->error("must not retry"))
    @test_throws ArgumentError bad("test")
end
