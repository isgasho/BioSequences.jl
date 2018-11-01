@testset "EachSkipmerIterator" begin
    seq = GeneralSequence{DNAAlphabet{2}}("GAGGGGGATGCCCCTCTTTGAGCCCAAGG")
    ST = Skipmer{DNAAlphabet{2},2,3,14}
    ans = [ST("CAAGGGGCTCCCTC"), ST("AGGGATCCCTTTGA"), ST("CTAAGAGGCACCCC"),
           ST("GCCAAGGGGCTCCC"), ST("GGATCCCTTTGACC"), ST("GGCTAAGAGGCACC"),
           ST("GAGCCCCTTGGCCA"), ST("ATCCCTTTGACCAA"), ST("CTGGCTAAGAGGCA"),
           ST("CCTGGCCAAGGGGC")]
    
    @test collect(BioSequences.EachSkipmerIterator(ST, seq)) == ans
end