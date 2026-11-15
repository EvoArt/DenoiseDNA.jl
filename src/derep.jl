# Dereplication (port of dada2::derepFastq).

"""
    Derep

Dereplicated reads of one sample.

- `sequence`: unique sequences, most abundant first (ties keep lexical order).
- `abundance`: number of reads of each unique.
- `quals`: `nunique x maxlen` matrix of mean Phred quality per position
  (`NaN` past the end of a sequence).
- `map`: for each input read, the index of its unique.
"""
struct Derep
    sequence::Vector{String}
    abundance::Vector{Int}
    quals::Matrix{Float64}
    map::Vector{Int}
end
Base.length(d::Derep) = length(d.sequence)

# Biostrings orders DNA letters by their internal codes (A<C<G<T<...<N),
# which differs from byte order only for ambiguity codes.
const _BSTRING_DNA_CODE = let t = UInt16.(0:255) .+ 0x0100   # unknown bytes sort after DNA letters
    for (c, v) in zip("ACMGRSVTWYHKDBN-+.", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 32, 64))
        t[UInt8(c)+1] = v
    end
    t
end

function _bstring_isless(a::Vector{UInt8}, b::Vector{UInt8})
    n = min(length(a), length(b))
    @inbounds for i in 1:n
        x = _BSTRING_DNA_CODE[a[i]+1]; y = _BSTRING_DNA_CODE[b[i]+1]
        x != y && return x < y
    end
    return length(a) < length(b)
end

"""
    derepFastq(file; n=1_000_000) -> Derep
    derepFastq(files::AbstractVector; n=1_000_000, progress) -> Vector{Derep}

Dereplicate a FASTQ file. Reads are processed in chunks of `n` (as R does);
within a chunk uniques are ordered lexically, new uniques from later chunks
are appended, and the result is stably sorted by decreasing abundance.
"""
function derepFastq(file::AbstractString; n::Integer=1_000_000)
    fq = read_fastq(file)
    return derep_records(fq.seq, fq.qual; n=n)
end
function derepFastq(files::AbstractVector{<:AbstractString}; progress::Bool=_progress_default(), kwargs...)
    bar = _progress(length(files), "Dereplicating"; unit=" files", enabled=progress)
    out = Derep[]
    for f in files
        push!(out, derepFastq(f; kwargs...))
        _next!(bar)
    end
    _finish!(bar; info="$(_fmt(sum(d -> sum(d.abundance), out; init=0))) reads → $(_fmt(sum(length, out; init=0))) uniques")
    return out
end

function derep_records(seqs::Vector{Vector{UInt8}}, quals::Vector{Vector{UInt8}}; n::Integer=1_000_000)
    nread = length(seqs)
    keep = findall(s -> !isempty(s), seqs)
    length(keep) == 0 && error("Only zero-length sequences detected during dereplication.")
    length(keep) < nread && @warn "Zero-length sequences detected during dereplication. They will be ignored."

    index = Dict{Vector{UInt8},Int}()
    useqs = Vector{Vector{UInt8}}()
    counts = Int[]
    qsums = Vector{Vector{Float64}}()
    rmap = zeros(Int, nread)

    for chunk in Iterators.partition(keep, n)
        # uniques first seen in this chunk, in lexical order
        newseqs = Vector{UInt8}[]
        for r in chunk
            s = seqs[r]
            if !haskey(index, s)
                index[s] = 0
                push!(newseqs, s)
            end
        end
        sort!(newseqs; lt=_bstring_isless)
        for s in newseqs
            push!(useqs, s); push!(counts, 0); push!(qsums, zeros(Float64, length(s)))
            index[s] = length(useqs)
        end
        for r in chunk
            u = index[seqs[r]]
            counts[u] += 1
            q = qsums[u]; qr = quals[r]
            @inbounds for k in eachindex(q)
                q[k] += qr[k]
            end
            rmap[r] = u
        end
    end

    ord = sortperm(counts; rev=true, alg=Base.Sort.DEFAULT_STABLE)
    rank = invperm(ord)
    maxlen = maximum(length, useqs)
    qm = fill(NaN, length(useqs), maxlen)
    for (k, u) in enumerate(ord)
        q = qsums[u]; c = counts[u]
        @inbounds for p in eachindex(q)
            qm[k, p] = q[p] / c
        end
    end
    map = [u == 0 ? 0 : rank[u] for u in rmap]
    return Derep(String.(copy.(useqs[ord])), counts[ord], qm, map)
end
