@testset "EachSkipmerIterator" begin
    seq = GeneralSequence{DNAAlphabet{2}}("GAGGGGGATGCCCCTCTTTGAGCCCAAGG")
    ST = Skipmer{DNAAlphabet{2},2,3,14}
    ans = [ST("CAAGGGGCTCCCTC"), ST("AGGGATCCCTTTGA"), ST("CTAAGAGGCACCCC"),
           ST("GCCAAGGGGCTCCC"), ST("GGATCCCTTTGACC"), ST("GGCTAAGAGGCACC"),
           ST("GAGCCCCTTGGCCA"), ST("ATCCCTTTGACCAA"), ST("CTGGCTAAGAGGCA"),
           ST("CCTGGCCAAGGGGC")]
           
    @test eltype(BioSequences.EachSkipmerIterator{ST,UInt64,typeof(seq)}) == ST
    @test Base.IteratorSize() == Base.HasLength()
    @test Base.IteratorEltype(BioSequences.EachSkipmerIterator{ST,UInt64,typeof(seq)}) == Base.HasEltype()
    @test BioSequences.kmersize(BioSequences.EachSkipmerIterator(ST, seq)) == 14
    @test BioSequences.firstoffset(BioSequences.EachSkipmerIterator(ST, seq)) == 26
    @test collect(BioSequences.EachSkipmerIterator(ST, seq)) == ans
    @test_throws ArgumentError BioSequences.EachSkipmerIterator(Skipmer{DNAAlphabet{2},1,3,14}, seq)
end