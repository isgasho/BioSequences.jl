@testset "CanonicalSkipmers" begin
    seq = GeneralSequence{DNAAlphabet{2}}("GAGGGGGATGCCCCTCTTTGAGCCCAAGG")
    ST = Skipmer{DNAAlphabet{2},2,3,14}
    ans = [ST("CAAGGGGCTCCCTC"), ST("AGGGATCCCTTTGA"), ST("CTAAGAGGCACCCC"),
           ST("GCCAAGGGGCTCCC"), ST("GGATCCCTTTGACC"), ST("GGCTAAGAGGCACC"),
           ST("GAGCCCCTTGGCCA"), ST("ATCCCTTTGACCAA"), ST("CTGGCTAAGAGGCA"),
           ST("CCTGGCCAAGGGGC")]
           
    @test eltype(BioSequences.CanonicalSkipmers{ST,UInt64,typeof(seq)}) == ST
    @test Base.IteratorSize(BioSequences.CanonicalSkipmers{ST,UInt64,typeof(seq)}) == Base.HasLength()
    @test Base.IteratorEltype(BioSequences.CanonicalSkipmers{ST,UInt64,typeof(seq)}) == Base.HasEltype()
    @test BioSequences.kmersize(BioSequences.CanonicalSkipmers(ST, seq)) == 14
    @test BioSequences.kmer_mask(BioSequences.CanonicalSkipmers(ST, seq)) == 0x000000000fffffff
    @test BioSequences.firstoffset(BioSequences.CanonicalSkipmers(ST, seq)) == 26
    @test collect(BioSequences.CanonicalSkipmers(ST, seq)) == ans
    @test_throws ArgumentError BioSequences.CanonicalSkipmers(Skipmer{DNAAlphabet{2},1,3,14}, seq)
end