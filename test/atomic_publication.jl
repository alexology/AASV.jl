@testset "Atomic rename never copies after a failure" begin
    mktempdir() do dir
        source, destination = joinpath(dir, "source"), joinpath(dir, "destination")
        write(source, "new")
        if Sys.iswindows()
            handle = ccall((:CreateFileW, "kernel32"), stdcall, Ptr{Cvoid},
                (Cwstring, UInt32, UInt32, Ptr{Cvoid}, UInt32, UInt32, Ptr{Cvoid}),
                source, 0x80000000, 0x00000001, C_NULL, 0x00000003, 0x00000000, C_NULL)
            @test handle != Ptr{Cvoid}(typemax(UInt))
            try
                @test_throws Base.IOError AASV.rename_atomic(source, destination)
                @test !isfile(destination)
                @test read(source, String) == "new"
            finally
                ccall((:CloseHandle, "kernel32"), stdcall, Cint, (Ptr{Cvoid},), handle)
            end
        end
        write(destination, "old")
        AASV.rename_atomic(source, destination)
        @test read(destination, String) == "new"
        @test !isfile(source)
    end
end
