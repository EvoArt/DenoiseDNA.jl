# Small utilities, ports of R dada2's misc.R, sequenceIO.R and paired.R.

# ------------------------------------------------------------------ sequences

"""
    getSequences(x; collapse=false) -> Vector{String}

The sequences of `x`, upper-cased: a FASTA/FASTQ file path, a vector of
sequences, or anything [`getUniques`](@ref) accepts. With `collapse`,
duplicates are removed (as R's `getSequences`).
"""
function getSequences(x; collapse::Bool=false)
    uppercase.(getUniques(x; collapse).sequence)
end
function getSequences(path::AbstractString; collapse::Bool=false)
    isfile(path) || return getSequences([path]; collapse)
    _isfasta(path) ? uppercase.(read_fasta(path).seq) : String.(read_fastq(path).seq)
end
getSequences(seqs::AbstractVector{<:AbstractString}; collapse::Bool=false) =
    uppercase.(collapse ? unique(seqs) : seqs)

"Total abundance of `x` (R's `dada2:::getAbund`)."
getAbund(x) = sum(getUniques(x).abundance)

function _isfasta(path)
    data = _readbytes(path)
    i = findfirst(c -> !isspace(Char(c)), data)
    return i !== nothing && data[i] == UInt8('>')
end

"""
    read_fasta(path) -> (id, seq)

Read a (possibly gzipped) FASTA file. Multi-line sequences are joined; `id`
is the header line without the '>'.
"""
function read_fasta(path::AbstractString)
    ids = String[]; seqs = String[]
    buf = IOBuffer()
    for line in eachline(IOBuffer(_readbytes(path)))
        line = rstrip(line, '\r')
        if startswith(line, '>')
            isempty(ids) || push!(seqs, String(take!(buf)))
            push!(ids, line[2:end])
        elseif !isempty(ids)
            print(buf, strip(line))
        end
    end
    isempty(ids) || push!(seqs, String(take!(buf)))
    return (id=ids, seq=seqs)
end

"""
    uniquesToFasta(unqs, fout; ids=nothing, mode="w", width=20000)

Write the unique sequences of `unqs` to a FASTA file, with ids
`sq<i>;size=<abundance>;` unless `ids` are given (the uchime format).
Sequences are wrapped at `width` letters; `mode="a"` appends.
"""
function uniquesToFasta(unqs, fout::AbstractString; ids=nothing, mode::AbstractString="w",
                        width::Integer=20000)
    u = getUniques(unqs; collapse=false)
    ids === nothing && (ids = ["sq$i;size=$a;" for (i, a) in enumerate(u.abundance)])
    length(ids) == length(u.sequence) || throw(ArgumentError("need one id per sequence"))
    open(fout, mode) do io
        for (id, s) in zip(ids, u.sequence)
            println(io, '>', id)
            for k in 1:width:max(length(s), 1)
                println(io, SubString(s, k, min(k + width - 1, length(s))))
            end
        end
    end
    return nothing
end

# ------------------------------------------------------------------ alignment helpers

"""
    nwextract(query, ref; kwargs...) -> String

The part of `ref` that aligns to `query` (between the first and last aligned
query bases), without gaps. Keywords go to [`nwalign`](@ref).
"""
function nwextract(query::AbstractString, ref::AbstractString; kwargs...)
    a0, a1 = nwalign(query, ref; kwargs...)
    isbase(c) = c in ('A', 'C', 'G', 'T')
    lo = findfirst(isbase, a0)
    lo === nothing && return ""
    hi = findlast(isbase, a0)
    return replace(a1[lo:hi], '-' => "")
end

"""
    strdiff(s1, s2) -> DataFrame(pos, nt0, nt1)

Positions where two (equal-length) strings differ.
"""
function strdiff(s1::AbstractString, s2::AbstractString)
    x = collect(s1); y = collect(s2)
    n = isempty(x) || isempty(y) ? 0 : max(length(x), length(y))   # R recycles the shorter one
    pos = [i for i in 1:n if x[mod1(i, length(x))] != y[mod1(i, length(y))]]
    return DataFrame(pos=pos, nt0=string.(x[mod1.(pos, length(x))]), nt1=string.(y[mod1.(pos, length(y))]))
end

# ------------------------------------------------------------------ misc

"""
    checkConvergence(x) -> Vector{Float64}

Sum of absolute differences between the final error matrix and each input
error matrix of a [`learnErrors`](@ref) or self-consistent [`dada`](@ref)
result. It reaches zero when the estimate converged.
"""
checkConvergence(x::Union{ErrorModel,DadaResult}) = [sum(abs.(x.err_out .- e)) for e in x.err_in]

"""
    sameOrder(fnF, fnR) -> Bool

Whether two FASTQ files list their reads in the same order (ids compared up
to the first space).
"""
function sameOrder(fnF::AbstractString, fnR::AbstractString)
    idF = _id_stem.(read_fastq(fnF).id)
    idR = _id_stem.(read_fastq(fnR).id)
    matched = true
    # R compares chunks of 1e6 reads, recycling the shorter in the last chunk
    chunk = 1_000_000
    for start in 1:chunk:min(length(idF), length(idR))
        f = idF[start:min(start + chunk - 1, end)]
        r = idR[start:min(start + chunk - 1, end)]
        n = max(length(f), length(r))
        matched &= all(f[mod1(k, length(f))] == r[mod1(k, length(r))] for k in 1:n)
    end
    return matched
end

# ShortRead's trimTails(id, 1, " "): cut at the first character <= ' '.
function _id_stem(id::AbstractString)
    i = findfirst(c -> c <= ' ', id)
    return i === nothing ? id : id[1:prevind(id, i)]
end
