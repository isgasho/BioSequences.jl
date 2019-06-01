# K-mer
# =====
#
# Compact k-mer sequence type.
#
# A Kmer is a sequence ≤ 32nt, without any 'N's, packed in a single 64 bit
# value.  While BioSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Representation
# --------------
#
# Four kinds of nucleotides are encoded as follows:
#
#   nucleotide | binary
#   ---------- | ------
#       A      |   00
#       C      |   01
#       G      |   10
#     T / U    |   11
#
# NucleicAcids are filled from MSBs to LSBs and right-aligned so that all k-mers
# are lexicographically ordered. For example, the memory layout of "TACG" is:
#   64-bit: 0b 00 00 … 00 11 00 01 10
#    4-mer:                T  A  C  G

primitive type Kmer{T<:NucleicAcid, K} <: Sequence 64 end

const DNAKmer{K} = Kmer{DNA, K}
const RNAKmer{K} = Kmer{RNA, K}
const DNACodon = DNAKmer{3}
const RNACodon = RNAKmer{3}

#take the subsequence of kmer generate a new type using the index ind
#
# sub_seq(Kmer{DNA,8},3) generates a kmer of type Kmer{DNA,6} which contains the same nucleic acid sequence with
# the initial kmer starting from index 3
#
#
#

"""
Special way of subsequencing kmers since we define each kmer type separately
"""
sub_seq(kmer::Kmer{T,K},ind::Int64) where{T,K} = Kmer{T,K-ind+1}(String(kmer)[ind:end])
sub_seq2(kmer::Kmer{T,K},ind::Int64) where{T,K} = Kmer{T,K-ind+1}(String(kmer)[ind+1:end])



function Kmer(nts::T...) where {T<:NucleicAcid}
    return make_kmer(nts)
end

function DNACodon(x::DNA, y::DNA, z::DNA)
    return make_kmer((x, y, z))
end

function RNACodon(x::RNA, y::RNA, z::RNA)
    return make_kmer((x, y, z))
end


# Conversion
# ----------

function Kmer{T, K}(x::UInt64) where {T, K}
    checkkmer(Kmer{T,K})
    mask = ~UInt64(0) >> (64 - 2K)
    return reinterpret(Kmer{T, K}, x & mask)
end

UInt64(x::Kmer) = reinterpret(UInt64, x)
Base.convert(::Type{UInt64}, x::Kmer) = reinterpret(UInt64, x)

DNAKmer{K}(x::RNAKmer{K}) where {K} = reinterpret(DNAKmer{K}, x)
RNAKmer{K}(x::DNAKmer{K}) where {K} = reinterpret(RNAKmer{K}, x)

function Kmer{T,K}(seq::AbstractString) where {T,K}
    return make_kmer(Kmer{T,K}, seq)
end

function Kmer{T,K}(seq::BioSequence{A}) where {T,K,A<:DNAAlphabet}
    return make_kmer(Kmer{DNA,K}, seq)
end

function Kmer{T,K}(seq::BioSequence{A}) where {T,K,A<:RNAAlphabet}
    return make_kmer(Kmer{RNA,K}, seq)
end

Kmer{T,K}(x::Kmer{T,K}) where {T,K} = x
Kmer{T}(seq::AbstractString) where {T} = Kmer{T,length(seq)}(seq)
Kmer(seq::BioSequence{A}) where {A<:DNAAlphabet} = Kmer{DNA,length(seq)}(seq)
Kmer(seq::BioSequence{A}) where {A<:RNAAlphabet} = Kmer{RNA,length(seq)}(seq)
DNAKmer(seq::BioSequence{A}) where {A<:DNAAlphabet} = DNAKmer{length(seq)}(seq)
RNAKmer(seq::BioSequence{A}) where {A<:RNAAlphabet} = RNAKmer{length(seq)}(seq)

# create a kmer from a sequence whose elements are convertible to a nucleotide
function make_kmer(::Type{Kmer{T,K}}, seq) where {T,K}
    seqlen = length(seq)
    if seqlen > 32
        throw(ArgumentError("cannot create a k-mer loger than 32nt"))
    elseif seqlen != K
        throw(ArgumentError("cannot create a $(K)-mer from a sequence of length $(seqlen)"))
    end

    x = UInt64(0)
    for c in seq
        nt = convert(T, c)
        if isambiguous(nt)
            throw(ArgumentError("cannot create a k-mer with ambiguous nucleotides"))
        elseif isgap(nt)
            throw(ArgumentError("cannot create a k-mer with gaps"))
        end
        x = (x << 2) | UInt64(trailing_zeros(nt))
    end

    return Kmer{T,K}(x)
end

make_kmer(seq::NTuple{K,T}) where {K,T} = make_kmer(Kmer{T,K}, seq)

BioSequence(x::DNAKmer{K}) where {K} = DNASequence(x)
BioSequence(x::RNAKmer{K}) where {K} = RNASequence(x)
BioSequence{A}(x::DNAKmer{K}) where {A<:DNAAlphabet,K} = BioSequence{A}([nt for nt in x])
BioSequence{A}(x::RNAKmer{K}) where {A<:RNAAlphabet,K} = BioSequence{A}([nt for nt in x])
Base.convert(::Type{S}, seq::Kmer) where {S<:AbstractString} = S([Char(x) for x in seq])
Base.String(seq::Kmer) = convert(String, seq)


# Basic Functions
# ---------------

BioSymbols.alphabet(::Type{DNAKmer{k}}) where {k} = (DNA_A, DNA_C, DNA_G, DNA_T)
BioSymbols.alphabet(::Type{RNAKmer{k}}) where {k} = (RNA_A, RNA_C, RNA_G, RNA_U)

Base.hash(x::Kmer, h::UInt) = hash(UInt64(x), h)

kmersize(::Type{Kmer{T,k}}) where {T,k} = k
kmersize(kmer::Kmer) = kmersize(typeof(kmer))
Base.length(x::Kmer{T, K}) where {T,K} = kmersize(x)
Base.eltype(::Type{Kmer{T,k}}) where {T,k} = T

@inline function inbounds_getindex(x::Kmer{T,K}, i::Integer) where {T,K}
    return reinterpret(T, 0x01 << ((UInt64(x) >> 2(K - i)) & 0b11))
end

Base.summary(x::DNAKmer{k}) where {k} = string("DNA ", k, "-mer")
Base.summary(x::RNAKmer{k}) where {k} = string("RNA ", k, "-mer")

Base.:-(x::Kmer{T,K}, y::Integer) where {T,K} = Kmer{T,K}(UInt64(x) - y % UInt64)
Base.:+(x::Kmer{T,K}, y::Integer) where {T,K} = Kmer{T,K}(UInt64(x) + y % UInt64)
Base.:+(x::Integer, y::Kmer{T,K}) where {T,K} = y + x
Base.:(==)(x::Kmer{T,k}, y::Kmer{T,k}) where {T,k} = UInt64(x) == UInt64(y)
Base.isless(x::Kmer{T,K}, y::Kmer{T,K}) where {T,K} = isless(UInt64(x), UInt64(y))

function Base.typemin(::Type{Kmer{T,K}}) where {T,K}
    checkkmer(Kmer{T,K})
    return reinterpret(Kmer{T,K}, UInt64(0))
end

function Base.typemax(::Type{Kmer{T,K}}) where {T,K}
    checkkmer(Kmer{T,K})
    return reinterpret(Kmer{T,K}, ~UInt64(0) >> (64 - 2K))
end

@inline function checkkmer(::Type{Kmer{T,K}}) where {T,K}
    if !(1 ≤ K ≤ 32)
        throw(ArgumentError("the length K must be within 1..32"))
    end
end


# Other functions
# ---------------

"""
    complement(kmer::Kmer)

Return the complement of `kmer`.
"""
BioSymbols.complement(x::Kmer{T,k}) where {T,k} = Kmer{T,k}(~UInt64(x))

"""
    reverse(kmer::Kmer)

Return the reverse of `kmer`.
"""
Base.reverse(x::Kmer{T,k}) where {T,k} = Kmer{T,k}(nucrev2(UInt64(x)) >> (64 - 2k))

"""
    reverse_complement(kmer::Kmer)

Return the reverse complement of `kmer`
"""
reverse_complement(x::Kmer) = complement(reverse(x))

"""
    mismatches(a::Kmer, b::Kmer)

Return the number of mismatches between `a` and `b`.
"""
function mismatches(a::Kmer{T,k}, b::Kmer{T,k}) where {T,k}
    return count_nonzero_bitpairs(UInt64(a) ⊻ UInt64(b))
end

"""
    canonical(kmer::Kmer)

Return the canonical k-mer of `x`.

A canonical k-mer is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting k-mers in data that is not strand specific,
and thus observing k-mer is equivalent to observing its reverse complement.
"""
function canonical(x::Kmer)
    y = reverse_complement(x)
    return x < y ? x : y
end

function Base.rand(::Type{Kmer{T,k}}) where {T,k}
    return Kmer{T,k}(rand(UInt64))
end

function Base.rand(::Type{Kmer{T,k}}, size::Integer) where {T,k}
    return [rand(Kmer{T,k}) for _ in 1:size]
end


# K-mer neighbor
# --------------

# neighbors on a de Bruijn graph
struct KmerNeighborIterator{T, K}
    x::Kmer{T, K}
end

"""
    neighbors(kmer::Kmer)

Return an iterator through k-mers neighboring `kmer` on a de Bruijn graph.
"""
neighbors(x::Kmer{T,K}) where {T,K} = KmerNeighborIterator{T,K}(x)

Base.length(::KmerNeighborIterator) = 4
Base.eltype(::Type{KmerNeighborIterator{T,k}}) where {T,k} = Kmer{T,k}

function Base.iterate(it::KmerNeighborIterator{T, K}, i::UInt64=UInt64(0)) where {T,K}
    if i == 4
        return nothing
    else
        return Kmer{T,K}((UInt64(it.x) << 2) | i), i + 1
    end
end


# Counters
# --------

function gc_content(kmer::Kmer{T,k}) where {T,k}
    if k == 0
        return 0.0
    else
        return (count_g(kmer) + count_c(kmer)) / k
    end
end

function count_a(kmer::Kmer{T,k}) where {T,k}
    return count_a(reinterpret(UInt64, kmer)) - (32 - k)
end

function count_c(kmer::Kmer{T,k}) where {T,k}
    return count_c(reinterpret(UInt64, kmer))
end

function count_g(kmer::Kmer{T,k}) where {T,k}
    return count_g(reinterpret(UInt64, kmer))
end

function count_t(kmer::Kmer{T,k}) where {T,k}
    return count_t(reinterpret(UInt64, kmer))
end

# Count A, C, T/U, G respectively in a kmer stored in a UInt64
function count_a(x::UInt64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & 0x5555555555555555)
end
count_c(x::UInt64) = count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
count_g(x::UInt64) = count_ones(((x >>> 1) & (~x)) & 0x5555555555555555)
count_t(x::UInt64) = count_ones((x    & (x >>> 1)) & 0x5555555555555555)


# Shuffle
# -------

function Random.shuffle(kmer::Kmer{T,k}) where {T,k}
    # Fisher-Yates shuffle
    for i in 1:k-1
        j = rand(i:k)
        kmer = swap(kmer, i, j)
    end
    return kmer
end


# Swap two nucleotides at `i` and `j`.
function swap(kmer::Kmer{T,k}, i, j) where {T,k}
    i = 2k - 2i
    j = 2k - 2j
    b = convert(UInt64, kmer)
    x = ((b >> i) ⊻ (b >> j)) & UInt64(0x03)
    return Kmer{T,k}(b ⊻ ((x << i) | (x << j)))
end


# random k-mer generators
# -----------------------

"""
    generate_kmer(::Type{Kmer{T,K}},len::Int64)where{T,K}

Return a random kmer of  length len (len = K)

generate_kmer generates a random string of length k using the corresponding alphabet for the NucleicAcid type
and converts it into a Kmer variable
"""
function generate_kmer(::Type{Kmer{T,K}},len::Int64) where{T,K}
    if T!=DNA && T!=RNA
        throw(ArgumentError("Cannot generate a kmer of  type $(T)"))
    end
    if K!=len
        throw(ArgumentError("cannot create a $(K)-mer of length $(len)"))
    end
    alp = alphabet(A)
    str_alp = ""
    for a in alp
        str_alp = str_alp*string(convert(Char,a))
    end

    Kmer{T,K}(randstring(str_alp,len))
end

"""
    generate_random_kmers(::Type{X},len::Int64,size::Int64) where{X<:NucleicAcid}

Returns size k-mers of length len of type X

generate_random_kmers uses the generate_kmer function to generate a vector of k-mers in a single line.
This is useful to gain time during creating toy datasets for experimentations
Also acts in a similar spirit to monitors frequently used for hardwares to ensure the function works properly.
"""
function generate_random_kmers(::Type{X},len::Int64,size::Int64) where{X<:NucleicAcid}
    seqs = Vector{Kmer{X,len}}()
    for i in 1:size
        seq = generate_kmer(Kmer{X,len},len)
        push!(seqs,seq)
    end
    seqs
end


# String literal
# --------------

macro kmer_str(seq)
    return DNAKmer(remove_newlines(seq))
end
