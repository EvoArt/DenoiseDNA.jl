# Minimal, fast FASTQ reading/writing (4-line records, optionally gzipped).

"""
    FastqRecords

Reads from a FASTQ file. `seq[i]` and `qual[i]` (Phred scores, offset removed)
are byte vectors; `id[i]` is the header without the leading '@'.
"""
struct FastqRecords
    id::Vector{String}
    seq::Vector{Vector{UInt8}}
    qual::Vector{Vector{UInt8}}
end
Base.length(f::FastqRecords) = length(f.seq)

_isgz(path::AbstractString) = endswith(lowercase(path), ".gz")

function _readbytes(path::AbstractString)
    data = read(path)
    if length(data) >= 2 && data[1] == 0x1f && data[2] == 0x8b
        data = transcode(GzipDecompressor, data)
    end
    return data
end

@inline function _line_end(data::Vector{UInt8}, pos::Int)
    n = length(data)
    e = something(findnext(==(UInt8('\n')), data, pos), n + 1)
    stop = e - 1
    stop >= pos && data[stop] == UInt8('\r') && (stop -= 1)
    return stop, e + 1
end

"""
    read_fastq(path; offset=33) -> FastqRecords

Read a (possibly gzipped) FASTQ file.
"""
function read_fastq(path::AbstractString; offset::Integer=33)
    data = _readbytes(path)
    ids = String[]; seqs = Vector{UInt8}[]; quals = Vector{UInt8}[]
    n = length(data); pos = 1; off = UInt8(offset)
    while pos <= n
        # skip blank lines
        if data[pos] == UInt8('\n') || data[pos] == UInt8('\r')
            pos += 1; continue
        end
        data[pos] == UInt8('@') || error("Malformed FASTQ record in $path at byte $pos")
        hstart = pos + 1
        hstop, pos = _line_end(data, pos)
        sstart = pos; sstop, pos = _line_end(data, pos)
        _, pos = _line_end(data, pos)                  # '+' line
        qstart = pos; qstop, pos = _line_end(data, pos)
        (qstop - qstart) == (sstop - sstart) ||
            error("Sequence and quality lengths differ in $path (record $(length(seqs)+1))")
        push!(ids, String(data[hstart:hstop]))
        s = data[sstart:sstop]
        @inbounds for k in eachindex(s)
            c = s[k]
            (UInt8('a') <= c <= UInt8('z')) && (s[k] = c - 0x20)
        end
        push!(seqs, s)
        q = data[qstart:qstop]
        @inbounds for k in eachindex(q)
            q[k] -= off
        end
        push!(quals, q)
    end
    return FastqRecords(ids, seqs, quals)
end

"""
    FastqWriter(path; compress=endswith(path, ".gz"))

Streaming FASTQ writer. Use `write_record(w, id, seq, qual)` then `close(w)`.
"""
struct FastqWriter
    io::IO
    raw::IO
end
function FastqWriter(path::AbstractString; compress::Bool=_isgz(path))
    raw = open(path, "w")
    io = compress ? GzipCompressorStream(raw) : raw
    return FastqWriter(io, raw)
end
function write_record(w::FastqWriter, id::AbstractString, seq::AbstractVector{UInt8},
                      qual::AbstractVector{UInt8}; offset::Integer=33)
    io = w.io
    write(io, UInt8('@')); write(io, id); write(io, UInt8('\n'))
    write(io, seq); write(io, "\n+\n")
    off = UInt8(offset)
    for q in qual
        write(io, q + off)
    end
    write(io, UInt8('\n'))
    return nothing
end
Base.close(w::FastqWriter) = (close(w.io); w.io === w.raw || close(w.raw); nothing)
