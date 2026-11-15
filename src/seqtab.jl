# Sequence tables (port of dada2 makeSequenceTable / collapseNoMismatch).

"""
    SequenceTable

Samples (rows) by sequences (columns) count table, from
[`makeSequenceTable`](@ref).

- `samples`: sample names.
- `sequences`: the sequences, one per column.
- `table`: `nsamples x nsequences` read counts.

`st[:, keep]` keeps the columns where `keep` is `true`; `DataFrame(st)` gives
the long form (`sample, sequence, abundance`, zeros dropped).
"""
struct SequenceTable
    samples::Vector{String}
    sequences::Vector{String}
    table::Matrix{Int}
end
Base.size(st::SequenceTable) = size(st.table)

"""
    getUniques(x; collapse=true) -> (sequence, abundance)

Unique sequences and abundances of a dada result, derep, merged-pairs
DataFrame, NamedTuple, sequence table or FASTQ file. With `collapse`,
duplicated sequences are summed (and, as in R, then ordered alphabetically).
"""
getUniques(df::AbstractDataFrame; collapse::Bool=true) =
    _uniques(String.(df.sequence), Int.(df.abundance), collapse)
getUniques(nt::NamedTuple{(:sequence, :abundance)}; collapse::Bool=true) =
    _uniques(String.(nt.sequence), Int.(nt.abundance), collapse)
getUniques(st::SequenceTable; collapse::Bool=true) =
    _uniques(st.sequences, vec(sum(st.table; dims=1)), collapse)
getUniques(path::AbstractString; collapse::Bool=true) = getUniques(derepFastq(path); collapse)

_uniques(seqs, ab, collapse) = collapse ? _collapse(seqs, ab) : (sequence=seqs, abundance=ab)

function _collapse(seqs::Vector{String}, ab::Vector{Int})
    allunique(seqs) && return (sequence=seqs, abundance=ab)
    d = Dict{String,Int}()
    for (s, a) in zip(seqs, ab)
        d[s] = get(d, s, 0) + a
    end
    k = sort!(collect(keys(d)))      # R's tapply orders by factor levels
    return (sequence=k, abundance=[d[s] for s in k])
end

"""
    makeSequenceTable(samples; orderBy="abundance", names=nothing) -> SequenceTable

Build a sample-by-sequence table from per-sample results (dada results,
merged-pair data frames, ...). `samples` may be a vector, or a vector of
`name => result` pairs / a `Dict`.
"""
function makeSequenceTable(samples::AbstractVector; orderBy::Union{Nothing,AbstractString}="abundance",
                           names::Union{Nothing,AbstractVector}=nothing)
    if eltype(samples) <: Pair
        names = [string(first(p)) for p in samples]
        samples = [last(p) for p in samples]
    end
    unqs = [getUniques(s) for s in samples]
    index = Dict{String,Int}(); seqs = String[]
    for u in unqs, s in u.sequence
        if !haskey(index, s)
            push!(seqs, s); index[s] = length(seqs)
        end
    end
    tab = zeros(Int, length(unqs), length(seqs))
    for (i, u) in enumerate(unqs), (s, a) in zip(u.sequence, u.abundance)
        tab[i, index[s]] = a
    end
    if orderBy !== nothing
        key = orderBy == "abundance" ? vec(sum(tab; dims=1)) :
              orderBy == "nsamples" ? vec(sum(tab .> 0; dims=1)) :
              error("orderBy must be \"abundance\", \"nsamples\" or nothing")
        ord = sortperm(key; rev=true, alg=Base.Sort.DEFAULT_STABLE)
        tab = tab[:, ord]; seqs = seqs[ord]
    end
    snames = names === nothing ? [string(i) for i in 1:length(unqs)] : string.(names)
    return SequenceTable(snames, seqs, tab)
end
makeSequenceTable(d::AbstractDict; kwargs...) = makeSequenceTable(collect(d); kwargs...)

"Keep only the given sequence columns."
function Base.getindex(st::SequenceTable, ::Colon, keep::AbstractVector{Bool})
    SequenceTable(st.samples, st.sequences[keep], st.table[:, keep])
end

"Long-format DataFrame (sample, sequence, abundance) of the nonzero entries."
function DataFrames.DataFrame(st::SequenceTable)
    s = String[]; q = String[]; a = Int[]
    for i in eachindex(st.samples), j in eachindex(st.sequences)
        if st.table[i, j] > 0
            push!(s, st.samples[i]); push!(q, st.sequences[j]); push!(a, st.table[i, j])
        end
    end
    return DataFrame(sample=s, sequence=q, abundance=a)
end

"""
    collapseNoMismatch(st; minOverlap=20, orderBy="abundance", identicalOnly=false,
                       vec=true, band=-1, verbose=false)

Merge sequences that are identical up to shifts or length variation.
"""
function collapseNoMismatch(st::SequenceTable; minOverlap::Integer=20, orderBy="abundance",
                            identicalOnly::Bool=false, vec::Bool=true, band::Integer=-1,
                            verbose::Bool=false)
    # collapse duplicates
    if !allunique(st.sequences)
        idx = Dict{String,Int}(); seqs = String[]
        for s in st.sequences
            haskey(idx, s) || (push!(seqs, s); idx[s] = length(seqs))
        end
        tab = zeros(Int, length(st.samples), length(seqs))
        for (j, s) in enumerate(st.sequences)
            tab[:, idx[s]] .+= st.table[:, j]
        end
        st = SequenceTable(st.samples, seqs, tab)
    end
    identicalOnly && return st
    tot = vec_sum(st.table)
    ord = sortperm(tot; rev=true, alg=Base.Sort.DEFAULT_STABLE)
    col = Dict(s => j for (j, s) in enumerate(st.sequences))
    collapsed = zeros(Int, size(st.table))
    out = String[]
    for qi in ord
        query = st.sequences[qi]
        added = false
        prefix = first(query, minOverlap)
        for ref in out
            pref = first(ref, minOverlap)
            if occursin(prefix, ref) || occursin(pref, query)
                if nwhamming(query, ref; vec=vec, band=band) == 0
                    collapsed[:, col[ref]] .+= st.table[:, qi]
                    added = true
                    break
                end
            end
        end
        if !added
            collapsed[:, qi] .= st.table[:, qi]
            push!(out, query)
        end
    end
    keep = [s in Set(out) for s in st.sequences]
    res = SequenceTable(st.samples, st.sequences[keep], collapsed[:, keep])
    o = sortperm(vec_sum(res.table); rev=true, alg=Base.Sort.DEFAULT_STABLE)
    res = SequenceTable(res.samples, res.sequences[o], res.table[:, o])
    verbose && println("Output $(length(res.sequences)) collapsed sequences out of $(length(st.sequences)) input sequences.")
    return res
end
vec_sum(t) = vec(sum(t; dims=1))

"""
    nwhamming(s1, s2; band=-1, vec=false, kwargs...) -> Int

Number of mismatches and internal indels between `s1` and `s2` after an
ends-free alignment ([`nwalign`](@ref)).
"""
function nwhamming(s1::AbstractString, s2::AbstractString; band::Integer=-1, vec::Bool=false, kwargs...)
    a0, a1 = nwalign(s1, s2; band=band, vec=vec, kwargs...)
    ev = eval_pair(Vector{UInt8}(codeunits(a0)), Vector{UInt8}(codeunits(a1)))
    return ev.mismatch + ev.indel
end

"""
    mergeSequenceTables(tables...; repeats="error", orderBy="abundance", tryRC=false) -> SequenceTable

Combine sequence tables (for example from different sequencing runs) into
one. Sample names repeated across tables are an error unless `repeats` is
something else (e.g. `"sum"`), in which case their counts are summed. With
`tryRC`, a sequence whose reverse complement appeared earlier is turned
around and merged with it. Columns are ordered by `orderBy` (`"abundance"`,
`"nsamples"` or `nothing`).
"""
function mergeSequenceTables(tables::AbstractVector{SequenceTable}; repeats::AbstractString="error",
                             orderBy::Union{Nothing,AbstractString}="abundance", tryRC::Bool=false)
    length(tables) >= 2 || throw(ArgumentError("At least two sequence tables are expected"))
    names = reduce(vcat, [t.samples for t in tables])
    if !allunique(names)
        dups = unique(n for (i, n) in enumerate(names) if n in names[1:i-1])
        repeats == "error" &&
            error("Duplicated sample names detected in the sequence table row names: ", join(dups, ", "))
        names = unique(names)
        @info "Duplicated sample names detected in the sequence table row names."
    end
    seqs = unique(reduce(vcat, [t.sequences for t in tables]))
    if tryRC && length(seqs) > 1
        seen = Set{String}()
        rcseqs = Set{String}()
        for s in seqs
            rc(s) in seen && push!(rcseqs, s)
            push!(seen, s)
        end
        if !isempty(rcseqs)
            @info "Reverse complemented sequences detected and re-oriented."
            tables = [_reorient(t, rcseqs) for t in tables]
            filter!(s -> !(s in rcseqs), seqs)
        end
    end
    row = Dict(n => i for (i, n) in enumerate(names))
    col = Dict(s => j for (j, s) in enumerate(seqs))
    tab = zeros(Int, length(names), length(seqs))
    for t in tables, (j, s) in enumerate(t.sequences), (i, n) in enumerate(t.samples)
        tab[row[n], col[s]] += t.table[i, j]
    end
    if orderBy == "abundance" || orderBy == "nsamples"
        key = orderBy == "abundance" ? vec(sum(tab; dims=1)) : vec(sum(tab .> 0; dims=1))
        ord = sortperm(key; rev=true, alg=Base.Sort.DEFAULT_STABLE)
        tab = tab[:, ord]; seqs = seqs[ord]
    end
    return SequenceTable(names, seqs, tab)
end
mergeSequenceTables(t1::SequenceTable, t2::SequenceTable, more::SequenceTable...; kwargs...) =
    mergeSequenceTables([t1, t2, more...]; kwargs...)

function _reorient(t::SequenceTable, rcseqs)
    any(in(rcseqs), t.sequences) || return t
    seqs = [s in rcseqs ? rc(s) : s for s in t.sequences]
    return collapseNoMismatch(SequenceTable(t.samples, seqs, t.table); identicalOnly=true)
end
