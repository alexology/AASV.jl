using AASV

# The former implementation is retained here as a correctness/performance baseline.
function original_nearest(records)
    mat=AASV.alignment_matrix(records)
    nearest=fill(Inf,length(records))
    for i in 1:length(records)-1, j in i+1:length(records)
        d=AASV.aligned_distance(mat[i,:],mat[j,:])
        nearest[i]=min(nearest[i],d)
        nearest[j]=min(nearest[j],d)
    end
    nearest
end

function benchmark_case(n,unique_count,width)
    state=Ref(UInt32(73))
    alphabet=collect("ACDEFGHIKLMNPQRSTVWY---")
    residue=function()
        state[]=state[]*UInt32(1664525)+UInt32(1013904223)
        alphabet[Int(state[]%UInt32(length(alphabet)))+1]
    end
    pool=[String([residue() for _ in 1:width]) for _ in 1:unique_count]
    records=[Sequence(string(i),pool[mod1(i,unique_count)]) for i in 1:n]
    @assert original_nearest(records)==AASV.nearest_aligned_distances(records)
    GC.gc()
    old=@timed original_nearest(records)
    GC.gc()
    new=@timed AASV.nearest_aligned_distances(records)
    @assert old.value==new.value
    println((sequences=n,unique_sequences=unique_count,columns=width,
        old_seconds=old.time,new_seconds=new.time,speedup=old.time/new.time,
        old_bytes=old.bytes,new_bytes=new.bytes))
end

# Warm up both paths before measuring.
sample=[Sequence("a","AAAA"),Sequence("b","AAAT")]
original_nearest(sample)
AASV.nearest_aligned_distances(sample)
benchmark_case(400,400,220)
benchmark_case(1000,50,220)
