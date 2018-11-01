module TestBioSequences

using Test
using Random
using LinearAlgebra: normalize
import BioSymbols
using BioSequences
using IntervalTrees: IntervalValue
using StatsBase
using YAML
using BioCore

const codons = [
    "AAA", "AAC", "AAG", "AAU",
    "ACA", "ACC", "ACG", "ACU",
    "AGA", "AGC", "AGG", "AGU",
    "AUA", "AUC", "AUG", "AUU",
    "CAA", "CAC", "CAG", "CAU",
    "CCA", "CCC", "CCG", "CCU",
    "CGA", "CGC", "CGG", "CGU",
    "CUA", "CUC", "CUG", "CUU",
    "GAA", "GAC", "GAG", "GAU",
    "GCA", "GCC", "GCG", "GCU",
    "GGA", "GGC", "GGG", "GGU",
    "GUA", "GUC", "GUG", "GUU",
    "UAA", "UAC", "UAG", "UAU",
    "UCA", "UCC", "UCG", "UCU",
    "UGA", "UGC", "UGG", "UGU",
    "UUA", "UUC", "UUG", "UUU",
    # translatable ambiguities in the standard code
    "CUN", "CCN", "CGN", "ACN",
    "GUN", "GCN", "GGN", "UCN"
]

function random_translatable_rna(n)
    probs = fill(1.0 / length(codons), length(codons))
    cumprobs = cumsum(probs)
    r = rand()
    x = Vector{String}(undef, n)
    for i in 1:n
        x[i] = codons[searchsorted(cumprobs, rand()).start]
    end

    return string(x...)
end

function get_bio_fmt_specimens()
    path = joinpath(dirname(@__FILE__), "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
end

# The generation of random test cases...

function random_array(n::Integer, elements, probs)
    cumprobs = cumsum(probs)
    x = Vector{eltype(elements)}(undef, n)
    for i in 1:n
        x[i] = elements[searchsorted(cumprobs, rand()).start]
    end
    return x
end

# Return a random DNA/RNA sequence of the given length.
function random_seq(n::Integer, nts, probs)
    cumprobs = cumsum(probs)
    x = Vector{Char}(undef, n)
    for i in 1:n
        x[i] = nts[searchsorted(cumprobs, rand()).start]
    end
    return String(x)
end

function random_seq(::Type{A}, n::Integer) where {A<:Alphabet}
    # TODO: Resolve the use of symbols(A()).
    nts = symbols(A())
    probs = Vector{Float64}(undef, length(nts))
    fill!(probs, 1 / length(nts))
    return GeneralSequence{A}(random_seq(n, nts, probs))
end

function random_dna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'T', 'N'], probs)
end

function random_rna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'U', 'N'], probs)
end

function random_aa(len)
    return random_seq(len,
        ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X' ],
        push!(fill(0.049, 20), 0.02))
end

function intempdir(fn::Function, parent=tempdir())
    dirname = mktempdir(parent)
    try
        cd(fn, dirname)
    finally
        rm(dirname, recursive=true)
    end
end

function random_dna_kmer(len)
    return random_dna(len, [0.25, 0.25, 0.25, 0.25])
end

function random_rna_kmer(len)
    return random_rna(len, [0.25, 0.25, 0.25, 0.25])
end

function random_dna_kmer_nucleotides(len)
    return random_array(len, [DNA_A, DNA_C, DNA_G, DNA_T],
                        [0.25, 0.25, 0.25, 0.25])
end

function random_rna_kmer_nucleotides(len)
    return random_array(len, [RNA_A, RNA_C, RNA_G, RNA_U],
                        [0.25, 0.25, 0.25, 0.25])
end

function dna_complement(seq::AbstractString)
    seqc = Vector{Char}(undef, length(seq))
    complementer = Dict(zip("-ACGTSWYRKMDVHBN", "-TGCASWRYMKHBDVN"))
    for (i, c) in enumerate(seq)
        seqc[i] = complementer[c]
    end
    return String(seqc)
end

function rna_complement(seq::AbstractString)
    seqc = Vector{Char}(undef, length(seq))
    complementer = Dict(zip("-ACGUSWYRKMDVHBN", "-UGCASWRYMKHBDVN"))
    for (i, c) in enumerate(seq)
        seqc[i] = complementer[c]
    end
    return String(seqc)
end

function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
end

include("symbols.jl")
include("alphabets.jl")

@test gap(Char) == '-'

@testset "BioSequences" begin
    @test BioSequence() == DNASequence()
    @test BioSequences.encoded_data_eltype(DNASequence()) == UInt64
    @test Alphabet(DNASequence()) == DNAAlphabet{4}()
    @test BioSequences.BitsPerSymbol(DNASequence()) == BioSequences.BitsPerSymbol{4}()
    @test BioSequences.BitsPerSymbol(RNASequence()) == BioSequences.BitsPerSymbol{4}()
    @test BioSequences.bits_per_symbol(DNASequence()) == 4
    @test BioSequences.bits_per_symbol(RNASequence()) == 4
    a = dna"A-CG-G"; b = rna"A-CG-G"; c = aa"AK-MV-";
    @test BioSequences.symbols_per_data_element(a) == 16
    @test BioSequences.firstbitindex(a) == BioSequences.bitindex(a, 1)
    @test BioSequences.lastbitindex(a) == BioSequences.bitindex(a, lastindex(a))
    @test BioSequences.bindata_mask(a) == 0x000000000000000f
    @test eltype(a) == DNA
    @test eltype(b) == RNA
    @test eltype(typeof(a)) == DNA
    @test eltype(typeof(b)) == RNA
    @test size(a) == (length(a),)
    @test size(b) == (length(b),)
    @test ungap(a) == dna"ACGG"
    @test ungap(b) == rna"ACGG"
    @test ungap(c) == aa"AKMV"
    @test ungap!(a) === a && a == dna"ACGG"
    @test ungap!(b) === b && b == rna"ACGG"
    @test ungap!(c) === c && c == aa"AKMV"
end

@testset "GeneralSequences" begin
    include("mutablesequences/conversion.jl")
    include("mutablesequences/basics.jl")
    include("mutablesequences/hashing.jl")
    include("mutablesequences/iteration.jl")
    include("mutablesequences/subseq.jl")
    include("mutablesequences/mutability.jl")
    include("mutablesequences/print.jl")
    include("mutablesequences/transformations.jl")
    include("mutablesequences/mutability.jl")
    include("mutablesequences/predicates.jl")
    include("mutablesequences/find.jl")
    include("mutablesequences/counting.jl")
    include("mutablesequences/gc_content.jl")
    include("mutablesequences/ambiguous.jl")
    include("mutablesequences/shuffle.jl")
end

@testset "ReferenceSequences" begin
    include("refseq/conversion.jl")
    include("refseq/basics.jl")
    include("refseq/long.jl")
    include("refseq/print.jl")
    include("refseq/random.jl")
end

include("composition.jl")

@testset "Skipmers" begin
    include("kmers/conversion.jl")
    include("kmers/comparisons.jl")
    include("kmers/length.jl")
    include("kmers/arithmetic.jl")
    include("kmers/access.jl")
    include("kmers/random.jl")
    include("kmers/find.jl")
    include("kmers/print.jl")
    include("kmers/transformations.jl")
    include("kmers/mismatches.jl")
    include("kmers/eachkmer.jl")
    include("kmers/debruijn_neighbors.jl")
    include("kmers/shuffle.jl")
end

@testset "Search" begin
    include("search/exact.jl")
    include("search/approximate.jl")
    include("search/regex.jl")
    include("search/pwm.jl")
end

include("translation.jl")

include("demultiplexer.jl")

@testset "Reading and Writing" begin
    include("io/FASTA.jl")
    include("io/FASTQ.jl")
    include("io/twobit.jl")
    include("io/abif.jl")
end

end
