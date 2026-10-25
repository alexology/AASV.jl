# Compile the ASV analysis and its reports on a small offline reference set,
# so the first call in a session does not wait for compilation.
@setup_workload begin
    dna="ATGGCTGCTGCTGCTGCTGCTGCTGCT"
    @compile_workload begin
        mktempdir() do dir
            nt=[Sequence("Genus_a__1",dna),Sequence("Genus_a__2",dna)]
            aa=[Sequence(r.id,"MAAAAAAAA") for r in nt]
            write_fasta(nt,joinpath(mkpath(joinpath(dir,"raw_sequences_nt","Family","species")),"Genus_a.fasta"))
            for (level,name) in (("species","Genus_a"),("genus","Genus"),("family","Family"))
                write_fasta(aa,joinpath(mkpath(joinpath(dir,"alignments","Family",level)),name*".fasta"))
            end
            input=DataFrame(family=["Family","Family"],genus=["Genus","Genus"],species=["Genus a","Genus a"],
                ASV=[dna,"ATGGCTGCTGATGCTGCTGCTGCTGCT"],sample1=[3,2])
            output=joinpath(dir,"results")
            asv_functional_structure!(input;output_dir=output,alignment_dir=dir,threads=1)
            classify_asv(output)
            check_reference_records(;alignment_dir=dir)
        end
    end
end
