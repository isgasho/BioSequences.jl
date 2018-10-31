@testset "Alphabets" begin

    @testset "Indexing" begin
        function test_generic_alphabet_indexing(T::Type)
            @test Base.IteratorSize(T) == Base.HasLength()
            @test Base.IteratorSize(T) == Base.HasLength()
            @test Base.IteratorEltype(T) == Base.HasEltype()
            @test Base.IteratorEltype(T) == Base.HasEltype()
            @test firstindex(T()) == 1
            @test firstindex(T()) == 1
        end
        
        function test_alphabet_lastindex(T::Type, val)
            @test lastindex(T()) == length(T()) == val
        end
        
        function test_alphabet_indexing(T::Type, symbols)
            for i in eachindex(T())
                @test T()[i] == symbols[i]
            end
        end
        
        function test_alphabet_iteration(T::Type, symbols)
            @test all(collect(T()) .== symbols)
        end
            
    
        @testset "DNA" begin
            dna4sym = (DNA_Gap, DNA_A, DNA_C, DNA_M,
                       DNA_G, DNA_R, DNA_S, DNA_V,
                       DNA_T, DNA_W, DNA_Y, DNA_H,
                       DNA_K, DNA_D, DNA_B, DNA_N)
                       
            dna2sym = (DNA_A, DNA_C, DNA_G, DNA_T)
                       
            test_generic_alphabet_indexing(DNAAlphabet{2})
            test_generic_alphabet_indexing(DNAAlphabet{4})
            
            test_alphabet_lastindex(DNAAlphabet{2}, 4)
            test_alphabet_lastindex(DNAAlphabet{4}, 16)
            
            @test eltype(DNAAlphabet{2}) == DNA
            @test eltype(DNAAlphabet{4}) == DNA
            
            test_alphabet_indexing(DNAAlphabet{2}, dna2sym)
            test_alphabet_indexing(DNAAlphabet{4}, dna4sym)
            
            test_alphabet_iteration(DNAAlphabet{4}, dna4sym)
            test_alphabet_iteration(DNAAlphabet{2}, dna2sym)
        end

        @testset "RNA" begin
            rna4sym = (RNA_Gap, RNA_A, RNA_C, RNA_M,
                       RNA_G, RNA_R, RNA_S, RNA_V,
                       RNA_U, RNA_W, RNA_Y, RNA_H,
                       RNA_K, RNA_D, RNA_B, RNA_N)
                   
            rna2sym = (RNA_A, RNA_C, RNA_G, RNA_U)
        
            test_generic_alphabet_indexing(RNAAlphabet{2})
            test_generic_alphabet_indexing(RNAAlphabet{4})
            
            test_alphabet_lastindex(RNAAlphabet{2}, 4)
            test_alphabet_lastindex(RNAAlphabet{4}, 16)
            
            @test eltype(RNAAlphabet{2}) == RNA
            @test eltype(RNAAlphabet{4}) == RNA
            
            test_alphabet_indexing(RNAAlphabet{2}, rna2sym)
            test_alphabet_indexing(RNAAlphabet{4}, rna4sym)
            
            test_alphabet_iteration(RNAAlphabet{4}, rna4sym)
            test_alphabet_iteration(RNAAlphabet{2}, rna2sym)
        end
        @testset "AminoAcid" begin
            aasym = (AA_A, AA_R, AA_N, AA_D,
                     AA_C, AA_Q, AA_E, AA_G,
                     AA_H, AA_I, AA_L, AA_K,
                     AA_M, AA_F, AA_P, AA_S,
                     AA_T, AA_W, AA_Y, AA_V,
                     AA_O, AA_U, AA_B, AA_J,
                     AA_Z, AA_X, AA_Term, AA_Gap)
            
            test_generic_alphabet_indexing(AminoAcidAlphabet)
            test_alphabet_lastindex(AminoAcidAlphabet, 28)
            @test eltype(AminoAcidAlphabet) == AminoAcid
            test_alphabet_indexing(AminoAcidAlphabet, aasym)
            test_alphabet_iteration(AminoAcidAlphabet, aasym)
        end
        @testset "Char" begin
            charsym = collect('\0':'\U10ffff')
            test_generic_alphabet_indexing(CharAlphabet)
            test_alphabet_lastindex(CharAlphabet, 1114112)
            @test eltype(CharAlphabet) == Char
            test_alphabet_indexing(CharAlphabet, charsym)
            test_alphabet_iteration(CharAlphabet, charsym)
        end
    end
    
    @testset "Promotion" begin
        @testset "DNA" begin
            @test Base.promote_rule(DNAAlphabet{2}, DNAAlphabet{4}) == DNAAlphabet{4}
            @test Base.promote_rule(DNAAlphabet{4}, DNAAlphabet{2}) == DNAAlphabet{4}
            @test Base.promote_rule(DNAAlphabet{4}, DNAAlphabet{4}) == DNAAlphabet{4}
            @test Base.promote_rule(DNAAlphabet{2}, DNAAlphabet{2}) == DNAAlphabet{2}
        end
        @testset "RNA" begin
            @test Base.promote_rule(RNAAlphabet{2}, RNAAlphabet{4}) == RNAAlphabet{4}
            @test Base.promote_rule(RNAAlphabet{4}, RNAAlphabet{2}) == RNAAlphabet{4}
            @test Base.promote_rule(RNAAlphabet{4}, RNAAlphabet{4}) == RNAAlphabet{4}
            @test Base.promote_rule(RNAAlphabet{2}, RNAAlphabet{2}) == RNAAlphabet{2}
        end
    end
    
    @testset "Minimal Alphabet" begin
        @test BioSequences.minimal_alphabet(DNAAlphabet{4}) == DNAAlphabet{2}
        @test BioSequences.minimal_alphabet(DNAAlphabet{2}) == DNAAlphabet{2}
        @test BioSequences.minimal_alphabet(RNAAlphabet{4}) == RNAAlphabet{2}
        @test BioSequences.minimal_alphabet(RNAAlphabet{2}) == RNAAlphabet{2}
    end
    
    @testset "BitsPerSymbol" begin
        BitsPerSymbol = BioSequences.BitsPerSymbol
        @testset "DNA" begin
            @test BitsPerSymbol(DNAAlphabet{2}()) == BitsPerSymbol{2}()
            @test BitsPerSymbol(DNAAlphabet{4}()) == BitsPerSymbol{4}()
        end
        @testset "RNA" begin
            @test BitsPerSymbol(RNAAlphabet{2}()) == BitsPerSymbol{2}()
            @test BitsPerSymbol(RNAAlphabet{4}()) == BitsPerSymbol{4}()
        end
        @testset "AminoAcid" begin
            @test BitsPerSymbol(AminoAcidAlphabet()) == BitsPerSymbol{8}()
        end
        @testset "Char" begin
            @test BitsPerSymbol(CharAlphabet()) == BitsPerSymbol{32}()
        end
        @testset "Void" begin
            @test BitsPerSymbol(BioSequences.VoidAlphabet()) == BitsPerSymbol{0}()
        end
    end
    
    @testset "Encoder" begin
        encode = BioSequences.encode
        EncodeError = BioSequences.EncodeError

        @testset "DNA" begin
            # 2 bits
            @test encode(DNAAlphabet{2}(), DNA_A) === 0x00
            @test encode(DNAAlphabet{2}(), DNA_C) === 0x01
            @test encode(DNAAlphabet{2}(), DNA_G) === 0x02
            @test encode(DNAAlphabet{2}(), DNA_T) === 0x03
            @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_M)
            @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_N)
            @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_Gap)

            # 4 bits
            for nt in BioSymbols.alphabet(DNA)
                @test encode(DNAAlphabet{4}(), nt) === reinterpret(UInt8, nt)
            end
            @test_throws EncodeError encode(DNAAlphabet{4}(), reinterpret(DNA, 0b10000))
        end

        @testset "RNA" begin
            # 2 bits
            @test encode(RNAAlphabet{2}(), RNA_A) === 0x00
            @test encode(RNAAlphabet{2}(), RNA_C) === 0x01
            @test encode(RNAAlphabet{2}(), RNA_G) === 0x02
            @test encode(RNAAlphabet{2}(), RNA_U) === 0x03
            @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_M)
            @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_N)
            @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_Gap)

            # 4 bits
            for nt in BioSymbols.alphabet(RNA)
                @test encode(RNAAlphabet{4}(), nt) === reinterpret(UInt8, nt)
            end
            @test_throws EncodeError encode(RNAAlphabet{4}(), reinterpret(RNA, 0b10000))
        end

        @testset "AminoAcid" begin
            @test encode(AminoAcidAlphabet(), AA_A) === 0x00
            for aa in BioSymbols.alphabet(AminoAcid)
                @test encode(AminoAcidAlphabet(), aa) === convert(UInt8, aa)
            end
            @test_throws BioSequences.EncodeError encode(AminoAcidAlphabet(), BioSymbols.AA_INVALID)
        end
    end
    
    
    @testset "Decoder" begin
        decode = BioSequences.decode
        DecodeError = BioSequences.DecodeError

        @testset "DNA" begin
            # 2 bits
            @test decode(DNAAlphabet{2}(), 0x00) === DNA_A
            @test decode(DNAAlphabet{2}(), 0x01) === DNA_C
            @test decode(DNAAlphabet{2}(), 0x02) === DNA_G
            @test decode(DNAAlphabet{2}(), 0x03) === DNA_T
            @test_throws DecodeError decode(DNAAlphabet{2}(), 0x04)
            @test_throws DecodeError decode(DNAAlphabet{2}(), 0x0e)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(DNAAlphabet{4}(), x) === reinterpret(DNA, x)
            end
            @test_throws DecodeError decode(DNAAlphabet{4}(), 0b10000)
        end

        @testset "RNA" begin
            # 2 bits
            @test decode(RNAAlphabet{2}(), 0x00) === RNA_A
            @test decode(RNAAlphabet{2}(), 0x01) === RNA_C
            @test decode(RNAAlphabet{2}(), 0x02) === RNA_G
            @test decode(RNAAlphabet{2}(), 0x03) === RNA_U
            @test_throws DecodeError decode(RNAAlphabet{2}(), 0x04)
            @test_throws DecodeError decode(RNAAlphabet{2}(), 0x0e)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(RNAAlphabet{4}(), x) === reinterpret(RNA, x)
            end
            @test_throws DecodeError decode(RNAAlphabet{4}(), 0b10000)
        end

        @testset "AminoAcid" begin
            @test decode(AminoAcidAlphabet(), 0x00) === AA_A
            for x in 0x00:0x1b
                @test decode(AminoAcidAlphabet(), x) === convert(AminoAcid, x)
            end
            @test_throws BioSequences.DecodeError decode(AminoAcidAlphabet(), 0x1c)
        end
    end

end