
struct CanonicalSkipmers{SK <: Skipmer, UT <: Unsigned, SQ <: BioSequence}
    seq::SQ
    cycle_pos::Vector{UInt8}
    last_unknown::Vector{Int64}
    fkmer::Vector{UT}
    rkmer::Vector{UT}
end

@inline Base.IteratorSize(::Type{T}) where T <: CanonicalSkipmers = Base.HasLength()
@inline Base.IteratorEltype(::Type{T}) where T <: CanonicalSkipmers = Base.HasEltype()
@inline Base.eltype(::Type{CanonicalSkipmers{SK, UT, SQ}}) where {SK <: Skipmer, UT <: Unsigned, SQ <: BioSequence} = SK
@inline Base.length(it::CanonicalSkipmers) = length(it.seq) - span(eltype(it)) + 1

@inline kmersize(x::CanonicalSkipmers) = kmersize(eltype(x))

@inline firstoffset(x::CanonicalSkipmers) = (kmersize(x) - 1) * 2

@inline function kmer_mask(x::CanonicalSkipmers{SK,UT}) where {SK <: Skipmer, UT <: Unsigned}
    return (UT(1) << (kmersize(SK) * 2)) - 1
end

function CanonicalSkipmers(::Type{SK}, seq::SQ) where {SK <: Skipmer, SQ <: BioSequence}
    checkskipmer(SK)
    if span(SK) > length(seq)
        throw(ArgumentError(string("The span of ", SK, " (", span(SK), ") is greater than the input sequence length (", length(seq), ").")))
    end
    last_unknown = Vector{Int64}(undef, cycle_len(SK))
    fkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    rkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    cycle_pos = Vector{UInt8}(undef, cycle_len(SK))
    return CanonicalSkipmers{SK, encoded_data_eltype(SK), SQ}(seq, cycle_pos, last_unknown, fkmer, rkmer)
end

@inline function init_iterator!(it::CanonicalSkipmers)
    N = cycle_len(eltype(it))
    @inbounds for i in 1:N
        it.cycle_pos[i] = N - i
        it.last_unknown[i] = -1
        it.fkmer[i] = 0
        it.rkmer[i] = 0
    end
end

@inline function _consider_position!(it::CanonicalSkipmers{SK, UT, SQ}, pos) where
        {SK, UT, A <: NucleicAcidAlphabet{2}, SQ <: BioSequence{A}}
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

function Base.iterate(it::CanonicalSkipmers{SK, UT, SQ}) where
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

function Base.iterate(it::CanonicalSkipmers{SK, UT, SQ}, state::Tuple{UInt, UInt}) where
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

function Base.iterate(it::CanonicalSkipmers{SK, UT, SQ}) where 
        {SK, UT, A <: NucleicAcidAlphabet{4}, SQ <: BioSequence{A}}
    N = cycle_len(eltype(it))
    M = bases_per_cycle(eltype(it))
    S = span(eltype(it))
    init_iterator!(it)
    pos = firstindex(it.seq)
    lastpos = lastindex(seq)
    
    while pos <= lastpos
        
        for ni in 1:N
            cycle_pos[ni] += 1
            if cycle_pos[ni] == N
                cycle_pos[ni] = 0
            end
            
            if cycle_pos[ni] < M
                println("Sequence position: ", p, ", Phase: ", ni)
                fbits = BioSequences.twobitnucs[reinterpret(UInt8, seq[pos]) + 0x01]
                if fbits == 0xFF
                    it.last_unknown[ni] = pos
                    fbits = 0x00
                end
                rbits = ~fbits & 0x03
                fkmer[ni] = ((fkmer[ni] << 2) | fbits) & kmer_mask
                rkmer[ni] = (rkmer[ni] >> 2) | (UInt64(rbits) << firstoffset)
            end
        end
        
        # If we are at pos, the skip-mer that started at pos-S is now done. 
        if pos >= S
            if pos == S
                fi = 0x01
            else
                fi += 0x01
                if fi == (N + 1)
                    fi = 0x01
                end
            end
            if last_unknown[fi] + S <= pos
                outkmer = ifelse(fkmer[fi] <= rkmer[fi], fkmer[fi], rkmer[fi])
                return reinterpret(SK, outkmer), (pos + 1, fi)
            end
        end
        
        pos += 1
        
    end
    
    return nothing
    
end


