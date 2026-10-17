@testset "Persistent family alignment checks" begin
    mktempdir() do dir
        path=joinpath(dir,"alignments","Family","species","Genus_species.fasta")
        write_fasta([Sequence("a","AC-DE"),Sequence("b","ACXDE")],path)
        checkpoint=joinpath(dir,"checkpoints","check_alignment","alignments","species","Family.json")
        @test check_alignments!(;alignment_dir=dir,show_progress=false)==[path]
        @test isfile(checkpoint)
        decorate=function(path)
            state=JSON3.read(read(path,String))
            payload=Dict{String,Any}(String(k)=>v for (k,v) in pairs(state))
            payload["sentinel"]="reused"
            AASV.write_json_checkpoint(path,payload)
        end
        decorate(checkpoint)
        @test check_alignments!(;alignment_dir=dir,show_progress=false)==[path]
        @test haskey(JSON3.read(read(checkpoint,String)),:sentinel)
        @test check_alignments!(;alignment_dir=dir,resume=false,show_progress=false)==[path]
        @test !haskey(JSON3.read(read(checkpoint,String)),:sentinel)
        decorate(checkpoint)
        # Same file size, different content: SHA-256 must invalidate the result.
        write_fasta([Sequence("a","ACXDE"),Sequence("b","ACXDE")],path)
        @test isempty(check_alignments!(;alignment_dir=dir,show_progress=false))
        @test !haskey(JSON3.read(read(checkpoint,String)),:sentinel)
        added=joinpath(dirname(path),"Other.fasta")
        write_fasta([Sequence("a","A-CD"),Sequence("b","AXCD")],added)
        @test check_alignments!(;alignment_dir=dir,show_progress=false)==[added]
        rm(added)
        @test isempty(check_alignments!(;alignment_dir=dir,show_progress=false))
        write_fasta([Sequence("a","A-CD"),Sequence("b","AXCD")],
            joinpath(dir,"alignments","OtherFamily","species","Other.fasta"))
        @test isempty(check_alignments!(;alignment_dir=dir,families=["Family"],show_progress=false))
        @test length(check_alignments!(;alignment_dir=dir,show_progress=false))==1
        @test isempty(check_alignments!(;alignment_dir=dir,tax_lev=GenusRank,show_progress=false))
        @test_throws ArgumentError check_alignments!(;alignment_dir=dir,type="invalid")
        # A corrupt checkpoint must not prevent recomputing the analysis.
        write(checkpoint,"{")
        @test_logs (:warn,r"Unreadable reference checkpoint") check_alignments!(;
            alignment_dir=dir,families="Family",show_progress=false)
    end
end

@testset "Deterministic cluster shards" begin
    families=["Delta","Alpha","Gamma","Beta","Alpha"]
    shards=[select_shard(families;shard_index=i,shard_count=3) for i in 1:3]
    @test shards == [["Alpha","Gamma"],["Beta"],["Delta"]]
    @test sort(vcat(shards...)) == sort(unique(families))
    @test isempty(intersect(shards[1],shards[2]))
    @test_throws ArgumentError select_shard(families;shard_index=0,shard_count=3)
    @test_throws ArgumentError select_shard(families;shard_index=4,shard_count=3)
    @test_throws ArgumentError select_shard(families;shard_count=0)
end

@testset "Persistent family divergence filter" begin
    mktempdir() do dir
        path=joinpath(dir,"alignments","Family","species","Genus_species.fasta")
        original=[Sequence("a","AAAA"),Sequence("b","AAAA"),Sequence("c","AAAT"),Sequence("d","WWWW")]
        write_fasta(original,path)
        cache=joinpath(dir,"checkpoints","filter_divergent","species","Family")
        checkpoint=joinpath(cache,"complete.json")
        report=joinpath(dir,"logs","divergent_references","species","Family.tsv")
        run_filter(;kwargs...)=filter_divergent_references!(;alignment_dir=dir,
            max_distance=0.3,min_seqs=2,show_progress=false,kwargs...)
        @test run_filter()==[path]
        @test length(read_fasta(path))==3
        state=JSON3.read(read(checkpoint,String))
        @test String(only(state.removed).sequence)=="d"
        @test Dict(String(k)=>String(v) for (k,v) in pairs(state.fingerprint))==AASV.reference_fingerprint([path])
        before_report=read(report,String)
        payload=Dict{String,Any}(String(k)=>v for (k,v) in pairs(state))
        payload["sentinel"]="reused"
        AASV.write_json_checkpoint(checkpoint,payload)
        @test run_filter()==[path]
        @test haskey(JSON3.read(read(checkpoint,String)),:sentinel)
        @test read(report,String)==before_report
        rm(report)
        run_filter()
        @test read(report,String)==before_report
        # A source replacement invalidates the family and is filtered again.
        write_fasta(original,path)
        run_filter()
        @test length(read_fasta(path))==3
        @test !haskey(JSON3.read(read(checkpoint,String)),:sentinel)
        @test nrow(CSV.read(report,DataFrame;delim='\t'))==1
        # A changed threshold must repeat filtering even with identical inputs.
        filter_divergent_references!(;alignment_dir=dir,max_distance=0.1,min_seqs=2,show_progress=false)
        @test length(read_fasta(path))==2
        @test JSON3.read(read(checkpoint,String)).parameters.max_distance==0.1
        @test only(CSV.read(report,DataFrame;delim='\t').sequence)=="c"
        @test isempty(filter_divergent_references!(;alignment_dir=dir,tax_lev=GenusRank,show_progress=false))
        @test_throws ArgumentError filter_divergent_references!(;alignment_dir=dir,max_distance=2)
        @test_throws ArgumentError filter_divergent_references!(;alignment_dir=dir,progress_interval=0)
    end
end

@testset "Interrupted filter installation recovery" begin
    mktempdir() do dir
        directory=joinpath(dir,"alignments","Family","species")
        paths=[joinpath(directory,"One.fasta"),joinpath(directory,"Two.fasta")]
        original=[Sequence("a","AAAA"),Sequence("b","AAAA"),Sequence("c","AAAT"),Sequence("d","WWWW")]
        for path in paths
            write_fasta(original,path)
        end
        before=AASV.reference_fingerprint(paths)
        cache=joinpath(dir,"checkpoints","filter_divergent","species","Family")
        mkpath(cache)
        stage=mktempdir(cache;prefix="prepared_",cleanup=false)
        for path in paths
            write_fasta(original[1:3],joinpath(stage,basename(path)))
        end
        after=AASV.reference_fingerprint(AASV.fasta_files(stage))
        rows=[(file=basename(path),taxon=splitext(basename(path))[1],sequence="d",
            distance=1.0,n_start=4,n_retained=3) for path in paths]
        pending=(version=AASV.DIVERGENT_FILTER_VERSION,
            parameters=(rank="species",max_distance=0.3,min_seqs=4),
            before=before,fingerprint=after,removed=rows,stage=basename(stage))
        AASV.write_json_checkpoint(joinpath(cache,"pending.json"),pending)
        # Simulate termination after only the first output was installed.
        cp(joinpath(stage,basename(first(paths))),first(paths);force=true)
        # Windows may also be interrupted between removing and renaming a file.
        rm(last(paths))
        @test filter_divergent_references!(;alignment_dir=dir,max_distance=0.3,
            show_progress=false)==paths
        @test all(length(read_fasta(path))==3 for path in paths)
        @test !isfile(joinpath(cache,"pending.json"))
        @test !isdir(stage)
        state=JSON3.read(read(joinpath(cache,"complete.json"),String))
        @test length(state.removed)==2
        report=joinpath(dir,"logs","divergent_references","species","Family.tsv")
        @test nrow(CSV.read(report,DataFrame;delim='\t'))==2
    end
end
