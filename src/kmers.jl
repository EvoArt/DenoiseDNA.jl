# K-mer screens used to skip alignments between distant sequences.
# Ports of dada2's kmers.cpp with the semantics of the default SSE=2 path:
#  * kmer distance = 1 - sum(min(kv1, kv2)) / (min(len1,len2) - k + 1)
#  * ordered-kmer distance compares kmers position by position over the
#    shorter sequence (the SSE implementation does NOT reject unequal lengths).

const KMER_SIZE = 5
const N_KMER = 4^KMER_SIZE

"Kmer index (0-based) of `seq[i:i+k-1]`, or -1 if it contains a non-ACGT base."
@inline function _kmer_at(seq::AbstractVector{UInt8}, i::Int, k::Int)
    kmer = 0
    @inbounds for j in i:i+k-1
        nt = Int(seq[j]) - 1
        (0 <= nt <= 3) || return -1
        kmer = 4kmer + nt
    end
    return kmer
end

"Fill `kv` (length 4^k) with kmer counts of `seq`."
function assign_kmer!(kv::AbstractVector{UInt16}, seq::AbstractVector{UInt8}, k::Int=KMER_SIZE)
    fill!(kv, 0)
    @inbounds for i in 1:length(seq)-k+1
        km = _kmer_at(seq, i, k)
        km >= 0 && (kv[km+1] += 0x0001)
    end
    return kv
end

"Fill `kord` with the kmer index at each position (0 for invalid kmers, as in C)."
function assign_kmer_order!(kord::AbstractVector{UInt16}, seq::AbstractVector{UInt8}, k::Int=KMER_SIZE)
    @inbounds for i in 1:length(seq)-k+1
        km = _kmer_at(seq, i, k)
        kord[i] = km >= 0 ? UInt16(km) : 0x0000
    end
    return kord
end

@inline function kmer_dist(kv1::AbstractVector{UInt16}, len1::Int, kv2::AbstractVector{UInt16}, len2::Int, k::Int=KMER_SIZE)
    dotsum = 0
    @inbounds @simd for i in eachindex(kv1)
        a = kv1[i]; b = kv2[i]
        dotsum += ifelse(a < b, a, b)
    end
    return _kdist(dotsum, len1, len2, k)
end

# (uint16 accumulator in C)
@inline _kdist(dotsum, len1, len2, k) = 1.0 - Float64(dotsum % UInt16) / ((min(len1, len2) - k) + 1.0)

"""
    KmerTable(seqs)

K-mer count vectors of a set of integer-coded sequences. Counts are stored in
8 bits, as in dada2's default SSE path, so a distance streams 1 KB per
sequence and takes 32 counts per instruction; if any count exceeds 254
(possible only for long reads) 16-bit counts are used instead, so distances
are always exact.
"""
struct KmerTable
    k8::Vector{UInt8}       # N_KMER counts per sequence, concatenated
    k16::Matrix{UInt16}     # N_KMER x n, or N_KMER x 0 when every count fits in 8 bits
end

function KmerTable(seqs::AbstractVector{<:AbstractVector{UInt8}})
    n = length(seqs)
    k16 = Matrix{UInt16}(undef, N_KMER, n)
    for (j, s) in enumerate(seqs)
        assign_kmer!(view(k16, :, j), s)
    end
    maximum(k16; init=0x0000) < 0x00ff || return KmerTable(UInt8[], k16)
    return KmerTable(vec(UInt8.(k16)), zeros(UInt16, N_KMER, 0))
end

"K-mer distance between sequences `i` and `j` of the table (lengths `li`, `lj`)."
@inline function kmer_dist(t::KmerTable, i::Int, li::Int, j::Int, lj::Int, k::Int=KMER_SIZE)
    if size(t.k16, 2) > 0
        return kmer_dist(view(t.k16, :, i), li, view(t.k16, :, j), lj, k)
    end
    a = t.k8
    oi = (i - 1) * N_KMER; oj = (j - 1) * N_KMER
    acc = Vec{32,UInt16}(0)
    @inbounds for c in 1:32:N_KMER
        m = min(vload(Vec{32,UInt8}, a, oi + c), vload(Vec{32,UInt8}, a, oj + c))
        acc += convert(Vec{32,UInt16}, m)
    end
    return _kdist(sum(acc), li, lj, k)
end

@inline function kord_dist(kord1::AbstractVector{UInt16}, len1::Int, kord2::AbstractVector{UInt16}, len2::Int, k::Int=KMER_SIZE)
    klen = min(len1, len2) - k + 1
    dotsum = 0
    @inbounds @simd for i in 1:klen
        dotsum += kord1[i] == kord2[i]
    end
    return _kdist(dotsum, len1, len2, k)
end
