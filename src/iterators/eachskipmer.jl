
struct EachSkipmerIterator{SK <: Skipmer, UT <: Unsigned, SQ <: BioSequence}
    seq::SQ
    cycle_pos::Vector{UInt8}
    last_unknown::Vector{Int64}
    fkmer::Vector{UT}
    rkmer::Vector{UT}
end

@inline Base.IteratorSize(::Type{T}) where T <: EachSkipmerIterator = Base.HasLength()
@inline Base.IteratorEltype(::Type{T}) where T <: EachSkipmerIterator = Base.HasEltype()
@inline Base.eltype(::Type{EachSkipmerIterator{SK, UT, SQ}}) where {SK <: Skipmer, UT <: Unsigned, SQ <: BioSequence} = SK
@inline Base.length(it::EachSkipmerIterator) = length(it.seq) - span(eltype(it)) + 1

@inline kmersize(x::EachSkipmerIterator) = kmersize(eltype(x))

@inline firstoffset(x::EachSkipmerIterator) = (kmersize(x) - 1) * 2

@inline function kmer_mask(x::EachSkipmerIterator{SK,UT}) where {SK <: Skipmer, UT <: Unsigned}
    return (UT(1) << (kmersize(SK) * 2)) - 1
end

function EachSkipmerIterator(::Type{SK}, seq::SQ) where {SK <: Skipmer, SQ <: BioSequence}
    checkskipmer(SK)
    if span(SK) > length(seq)
        throw(ArgumentError(string("The span of ", SK, " (", span(SK), ") is greater than the input sequence length (", length(seq), ").")))
    end
    last_unknown = Vector{Int64}(undef, cycle_len(SK))
    fkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    rkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    cycle_pos = Vector{UInt8}(undef, cycle_len(SK))
    return EachSkipmerIterator{SK, encoded_data_eltype(SK), SQ}(seq, cycle_pos, last_unknown, fkmer, rkmer)
end

@inline function init_iterator!(it::EachSkipmerIterator)
    N = cycle_len(eltype(it))
    @inbounds for i in 1:N
        it.cycle_pos[i] = N - i
        it.last_unknown[i] = -1
        it.fkmer[i] = 0
        it.rkmer[i] = 0
    end
end

@inline function _consider_position!(it::EachSkipmerIterator, pos)
    N = cycle_len(eltype(it))
    M = bases_per_cycle(eltype(it))
    for ni in 1:N
        it.cycle_pos[ni] += 1
        if it.cycle_pos[ni] == N
            it.cycle_pos[ni] = 0
        end
        if it.cycle_pos[ni] < M
            fbits = extract_encoded_symbol(bitindex(it.seq, pos), encoded_data(it.seq))
            rbits = ~fbits & 0x03
            it.fkmer[ni] = ((it.fkmer[ni] << 2) | fbits) & kmer_mask(it)
            it.rkmer[ni] = (it.rkmer[ni] >> 2) | (UInt64(rbits) << firstoffset(it))
        end
    end
end

function Base.iterate(it::EachSkipmerIterator{SK, UT, SQ}) where
        {SK, UT, A <: NucleicAcidAlphabet{2}, SQ <: BioSequence{A}}
    N = cycle_len(eltype(it))
    M = bases_per_cycle(eltype(it))
    S = span(eltype(it))
    init_iterator!(it)
    for pos in 1:S
        _consider_position!(it, pos)
    end
    fkmer = first(it.fkmer)
    rkmer = first(it.rkmer)
    outkmer = ifelse(fkmer < rkmer, fkmer, rkmer)
    return reinterpret(eltype(it), outkmer), (S + 1, UInt(1))
end

function Base.iterate(it::EachSkipmerIterator{SK, UT, SQ}, state::Tuple{UInt, UInt}) where
        {SK, UT, A <: NucleicAcidAlphabet{2}, SQ <: BioSequence{A}}
    pos = state[1]
    fi  = state[2]
    N = cycle_len(eltype(it))
    M = bases_per_cycle(eltype(it))
    if pos > lastindex(it.seq)
        return nothing
    end
    _consider_position!(it, pos)
    fi += 1
    fi = ifelse(fi == (N + 1), UInt(1), fi)
    fkmer = it.fkmer[fi]
    rkmer = it.rkmer[fi]
    outkmer = ifelse(fkmer < rkmer, fkmer, rkmer)
    return reinterpret(eltype(it), outkmer), (pos + 1, fi)  
end


