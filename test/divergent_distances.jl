@testset "Optimized divergence distances match exhaustive character comparisons" begin
    oracle=function(records)
        mat=AASV.alignment_matrix(records)
        nearest=fill(Inf,length(records))
        for i in 1:length(records)-1, j in i+1:length(records)
            distance=AASV.aligned_distance(mat[i,:],mat[j,:])
            nearest[i]=min(nearest[i],distance)
            nearest[j]=min(nearest[j],distance)
        end
        nearest
    end
    fixtures=[String[],[""],["",""],["----","----"],
        ["AAAA","AAAA","AAAT","WWWW"],
        ["AA--","--AA","----","AA--"],
        ["A--A","----","A--A","--AA"],
        ["--AA--","--AA--","AA----","----AA"],
        ["A-BC","AXBC","A-BC","----"]]
    for sequences in fixtures
        records=[Sequence(string(i),sequence) for (i,sequence) in enumerate(sequences)]
        @test AASV.nearest_aligned_distances(records)==oracle(records)
    end
    # Fixed local generator makes failures reproducible without changing global RNG state.
    state=Ref(UInt32(187))
    nextindex=function(n)
        state[]=state[]*UInt32(1664525)+UInt32(1013904223)
        Int(state[]%UInt32(n))+1
    end
    alphabet=collect("ACDEFGHIKLMNPQRSTVWY---")
    for trial in 1:100
        width=nextindex(70)
        unique_count=nextindex(20)
        pool=[String([alphabet[nextindex(length(alphabet))] for _ in 1:width]) for _ in 1:unique_count]
        records=[Sequence(string(i),pool[nextindex(unique_count)]) for i in 1:nextindex(45)]
        @test AASV.nearest_aligned_distances(records)==oracle(records)
    end
    @test_throws ArgumentError AASV.nearest_aligned_distances([Sequence("a","AA"),Sequence("b","A")])
end
