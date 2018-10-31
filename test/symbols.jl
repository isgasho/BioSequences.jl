# NOTE: Most tests related to biological symbols are located in BioSymbols.jl.
@testset "Symbols" begin
    @testset "DNA" begin
        @test DNA_A === BioSymbols.DNA_A
        @test ACGT  === BioSymbols.ACGT
        @test ACGTN === BioSymbols.ACGTN
        @test typeof(DNA_A) === BioSymbols.DNA
    end

    @testset "RNA" begin
        @test RNA_A === BioSymbols.RNA_A
        @test ACGU  === BioSymbols.ACGU
        @test ACGUN === BioSymbols.ACGUN
        @test typeof(RNA_A) === BioSymbols.RNA
    end

    @testset "AminoAcid" begin
        @test AA_A === BioSymbols.AA_A
        @test typeof(AA_A) === BioSymbols.AminoAcid
    end

    @testset "Predicate functions" begin
        @test iscompatible(DNA_A, DNA_N)
        @test isambiguous(DNA_N)
        @test iscertain(DNA_A)
        @test isgap(DNA_Gap)
        @test ispurine(DNA_A)
        @test ispyrimidine(DNA_C)
        @test isGC(DNA_G)
    end

    @testset "Misc. functions" begin
        @test length(BioSymbols.alphabet(DNA)) == 16
        @test BioSymbols.gap(DNA) === DNA_Gap
        @test BioSymbols.complement(DNA_A) === DNA_T
    end
end
