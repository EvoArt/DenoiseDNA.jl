# Read filtering and trimming (port of dada2 filterAndTrim / fastqFilter /
# fastqPairedFilter). Every filter in R is applied per read, and a pair is kept
# iff both of its reads pass, so each read is processed independently here.

const _PHIX = Ref{String}("")
function _phix_genome()
    if isempty(_PHIX[])
        lines = readlines(joinpath(@__DIR__, "..", "data", "phix_genome.fa"))
        _PHIX[] = uppercase(join(filter(l -> !startswith(l, ">"), lines)))
    end
    return _PHIX[]
end

# Words of up to 32 nt are 2-bit coded into a UInt64; words containing a
# non-ACGT character get the code `_NOWORD` and never match.
const _NOWORD = typemax(UInt64)

# Code of the word ending at each position of `s` (`_NOWORD` until a full,
# all-ACGT word has been seen).
function _word_codes!(codes::Vector{UInt64}, s, w::Int)
    b = codeunits(s)
    resize!(codes, length(b))
    mask = w == 32 ? typemax(UInt64) : (UInt64(1) << 2w) - 1
    code = UInt64(0); run = 0
    @inbounds for (i, c) in enumerate(b)
        nt = c == UInt8('A') ? 0 : c == UInt8('C') ? 1 : c == UInt8('G') ? 2 : c == UInt8('T') ? 3 : -1
        if nt < 0
            run = 0
        else
            code = ((code << 2) | UInt64(nt)) & mask
            run += 1
        end
        codes[i] = run >= w ? code : _NOWORD
    end
    return codes
end

function _ref_words(ref::String, w::Int)
    codes = _word_codes!(UInt64[], ref * ref[1:w], w)    # the genome is circular
    return Set(c for c in codes[w:length(ref)+w-1] if c != _NOWORD)
end

# Port of C_matchRef: number of reference words found in each sequence.
# With `non_overlapping`, a hit at position j resumes the scan at j + w + 1.
function _match_ref(seqs::AbstractVector, words::Set{UInt64}, w::Int, non_overlapping::Bool)
    hits = zeros(Int, length(seqs))
    codes = UInt64[]
    for (k, s) in enumerate(seqs)
        L = ncodeunits(s)
        L < w && continue
        _word_codes!(codes, s, w)
        j = 1
        while j <= L - w + 1
            c = codes[j+w-1]
            if c != _NOWORD && c in words
                hits[k] += 1
                non_overlapping && (j += w)
            end
            j += 1
        end
    end
    return hits
end

const _PHIX_WORDS = Dict{Int,Tuple{Set{UInt64},Set{UInt64}}}()
const _PHIX_LOCK = ReentrantLock()
_phix_words(w) = lock(_PHIX_LOCK) do
    get!(_PHIX_WORDS, w) do
        g = _phix_genome()
        (_ref_words(g, w), _ref_words(rc(g), w))
    end
end

"""
    isPhiX(seqs; wordSize=16, minMatches=2, nonOverlapping=true) -> Vector{Bool}

Flag sequences sharing at least `minMatches` `wordSize`-mers with the phiX genome
(either strand).
"""
function isPhiX(seqs::AbstractVector{<:AbstractString}; wordSize::Integer=16, minMatches::Integer=2,
                nonOverlapping::Bool=true)
    1 <= wordSize <= 32 || throw(ArgumentError("wordSize must be between 1 and 32"))
    fwd, rev = _phix_words(wordSize)
    h = _match_ref(seqs, fwd, wordSize, nonOverlapping)
    hrc = _match_ref(seqs, rev, wordSize, nonOverlapping)
    return (h .>= minMatches) .| (hrc .>= minMatches)
end

"""
    seqComplexity(seqs; kmerSize=2, window=nothing, by=5) -> Vector{Float64}

Effective number of kmers (exp of Shannon entropy of kmer frequencies).
"""
function seqComplexity(seqs::AbstractVector{<:AbstractString}; kmerSize::Integer=2,
                       window::Union{Nothing,Integer}=nothing, by::Integer=5)
    window !== nothing && kmerSize >= window &&
        error("The window over which kmer frequency is calculated must be larger than the kmerSize.")
    function sindex(s)
        counts = zeros(Int, 4^kmerSize)
        codes = encode_nt(s)
        for i in 1:length(codes)-kmerSize+1
            km = _kmer_at(codes, i, kmerSize)
            km >= 0 && (counts[km+1] += 1)
        end
        y = counts ./ sum(counts)
        y = y[y .> 0]
        return exp(sum(-y .* log.(y)))
    end
    window === nothing && return [sindex(s) for s in seqs]
    si = fill(Float64(4^kmerSize), length(seqs))
    maxw = maximum(length, seqs)
    for i in 1:by:(maxw - window)
        for (k, s) in enumerate(seqs)
            length(s) >= i + window - 1 && (si[k] = min(si[k], sindex(s[i:i+window-1])))
        end
    end
    return si
end

struct _FilterParams
    truncQ::Int
    start::Int
    stop::Int          # required length after trimLeft (0 = none)
    maxLen::Float64
    minLen::Int
    trimRight::Int
    maxN::Int
    minQ::Int
    maxEE::Float64
end

function _filter_params(truncQ, truncLen, trimLeft, trimRight, maxLen, minLen, maxN, minQ, maxEE)
    start = max(1, trimLeft + 1)
    stop = truncLen < start ? 0 : truncLen - start + 1
    return _FilterParams(truncQ, start, stop, Float64(maxLen), minLen, trimRight, maxN, minQ, Float64(maxEE))
end

# Returns the kept (start, stop) range within the read, or nothing.
function _filter_read(seq::Vector{UInt8}, q::Vector{UInt8}, F::_FilterParams, paired::Bool)
    w = length(seq)
    w <= F.maxLen || return nothing
    w >= F.start || return nothing
    a = F.start; b = w
    if F.trimRight > 0
        (b - a + 1) > F.trimRight || return nothing
        b -= F.trimRight
    end
    # truncate at the first quality <= truncQ
    for i in a:b
        if q[i] <= F.truncQ
            b = i - 1
            break
        end
    end
    len = b - a + 1
    paired && len <= 0 && return nothing
    if F.stop > 0
        len >= F.stop || return nothing
        b = a + F.stop - 1; len = F.stop
    end
    len >= F.minLen || return nothing
    nother = 0
    for i in a:b
        c = seq[i]
        (c == UInt8('A') || c == UInt8('C') || c == UInt8('G') || c == UInt8('T')) || (nother += 1)
    end
    nother <= F.maxN || return nothing
    if F.minQ > F.truncQ && len > 0
        minimum(view(q, a:b)) > F.minQ || return nothing
    end
    if F.maxEE < Inf
        ee = 0.0
        for i in a:b
            ee += 10.0^(-Int(q[i]) / 10.0)
        end
        ee <= F.maxEE || return nothing
    end
    return a:b
end

_pair(x) = x isa Union{Tuple,AbstractVector} ? (x[1], x[end]) : (x, x)

function _orient!(seqs, quals, fwd::String)
    L = length(fwd)
    keep = trues(length(seqs))
    for k in eachindex(seqs)
        s = String(copy(seqs[k]))
        if length(s) >= L && s[1:L] == fwd
            continue
        end
        r = rc(s)
        if length(r) >= L && r[1:L] == fwd
            seqs[k] = Vector{UInt8}(codeunits(r)); quals[k] = reverse(quals[k])
        else
            keep[k] = false
        end
    end
    return keep
end

"""
    filterAndTrim(fwd, filt, rev=nothing, filt_rev=nothing; truncQ=2, truncLen=0,
                  trimLeft=0, trimRight=0, maxLen=Inf, minLen=20, maxN=0, minQ=0,
                  maxEE=Inf, rm_phix=true, rm_lowcomplex=0, orient_fwd=nothing,
                  matchIDs=false, id_sep=r"\\s", id_field=nothing, compress=true,
                  multithread=false, verbose=false, progress) -> DataFrame

Filter and trim FASTQ files (port of R's `filterAndTrim`). For paired reads,
per-direction parameters can be given as 2-tuples/vectors `(fwd, rev)`.
Returns a DataFrame with `reads_in` and `reads_out` per input file.

`progress=false` hides the progress bar (default: [`DenoiseDNA.setprogress!`](@ref)).
"""
function filterAndTrim(fwd, filt, rev=nothing, filt_rev=nothing; truncQ=2, truncLen=0, trimLeft=0,
                       trimRight=0, maxLen=Inf, minLen=20, maxN=0, minQ=0, maxEE=Inf,
                       rm_phix=true, rm_lowcomplex=0, orient_fwd=nothing, matchIDs::Bool=false,
                       id_sep=r"\s", id_field=nothing, compress::Bool=true, multithread=false,
                       verbose::Bool=false, progress::Bool=_progress_default())
    fwd = fwd isa AbstractString ? [fwd] : collect(fwd)
    filt = filt isa AbstractString ? [filt] : collect(filt)
    if length(filt) == 1 && length(fwd) > 1
        filt = joinpath.(filt[1], basename.(fwd))
    end
    length(fwd) == length(filt) || error("Every input file must have a corresponding output file.")
    paired = rev !== nothing
    if paired
        rev = rev isa AbstractString ? [rev] : collect(rev)
        filt_rev === nothing && error("Output files for the reverse reads are required.")
        filt_rev = filt_rev isa AbstractString ? [filt_rev] : collect(filt_rev)
        length(filt_rev) == 1 && length(rev) > 1 && (filt_rev = joinpath.(filt_rev[1], basename.(rev)))
        length(rev) == length(fwd) == length(filt_rev) || error("Paired forward and reverse input files must correspond.")
    end
    for f in vcat(filt, paired ? filt_rev : String[])
        mkpath(dirname(abspath(f)))
    end
    allunique(vcat(filt, paired ? filt_rev : String[])) || error("All output files must be distinct.")
    ins = zeros(Int, length(fwd)); outs = zeros(Int, length(fwd))
    getp(x) = _pair(x)
    tq = getp(truncQ); tl = getp(truncLen); lft = getp(trimLeft); rgt = getp(trimRight)
    mx = getp(maxLen); mn = getp(minLen); nn = getp(maxN); mq = getp(minQ); ee = getp(maxEE)
    ph = getp(rm_phix); lc = getp(rm_lowcomplex)
    FF = _filter_params(tq[1], tl[1], lft[1], rgt[1], mx[1], mn[1], nn[1], mq[1], ee[1])
    FR = _filter_params(tq[2], tl[2], lft[2], rgt[2], mx[2], mn[2], nn[2], mq[2], ee[2])
    bar = _progress(length(fwd), "Filtering"; unit=paired ? " pairs" : " files", enabled=progress)
    run_file = function (i)
        if paired
            ins[i], outs[i] = _filter_pair(fwd[i], rev[i], filt[i], filt_rev[i], FF, FR,
                                           ph, lc, orient_fwd, matchIDs, id_sep, id_field, compress)
        else
            ins[i], outs[i] = _filter_single(fwd[i], filt[i], FF, ph[1], lc[1], orient_fwd, compress)
        end
        verbose && _progress_println(bar, "Read in $(ins[i]), output $(outs[i]) filtered sequences ($(basename(fwd[i]))).")
        _next!(bar)
    end
    if _nthreads(multithread) > 1
        Threads.@threads :dynamic for i in eachindex(fwd)
            run_file(i)
        end
    else
        foreach(run_file, eachindex(fwd))
    end
    _finish!(bar; info="$(_fmt(sum(outs))) of $(_fmt(sum(ins))) reads kept ($(_pct(sum(outs) / max(sum(ins), 1))))", ok=any(>(0), outs))
    all(==(0), outs) && @warn "No reads passed the filter. Please revisit your filtering parameters."
    return DataFrame(file=basename.(fwd), reads_in=ins, reads_out=outs)
end

function _apply_extra(seqs::Vector{String}, phix::Bool, lowc::Real)
    keep = trues(length(seqs))
    phix && (keep .&= .!isPhiX(seqs))
    lowc > 0 && (keep .&= seqComplexity(seqs) .>= lowc)
    return keep
end

function _filter_single(fin, fout, F::_FilterParams, phix::Bool, lowc, orient_fwd, compress)
    fq = read_fastq(fin)
    n = length(fq)
    seqs = fq.seq; quals = fq.qual
    okeep = orient_fwd === nothing ? trues(n) : _orient!(seqs, quals, String(orient_fwd))
    ranges = Vector{Union{Nothing,UnitRange{Int}}}(nothing, n)
    for k in 1:n
        okeep[k] && (ranges[k] = _filter_read(seqs[k], quals[k], F, false))
    end
    idx = findall(!isnothing, ranges)
    trimmed = String[String(seqs[k][ranges[k]]) for k in idx]
    extra = _apply_extra(trimmed, phix, lowc)
    idx = idx[extra]
    _write_kept(fout, fq, idx, ranges, compress)
    return n, length(idx)
end

function _filter_pair(finF, finR, foutF, foutR, FF, FR, ph, lc, orient_fwd, matchIDs, id_sep, id_field, compress)
    fF = read_fastq(finF); fR = read_fastq(finR)
    nin = length(fF)
    if matchIDs
        field = id_field === nothing ? _detect_id_field(fF.id[1], id_sep) : id_field
        key(id) = split(id, id_sep)[field]
        kF = [split(key(x), '#')[1] for x in fF.id]; kR = [split(key(x), '#')[1] for x in fR.id]
        sR = Set(kR); sF = Set(kF)
        selF = findall(in(sR), kF); selR = findall(in(sF), kR)
        fF = FastqRecords(fF.id[selF], fF.seq[selF], fF.qual[selF])
        fR = FastqRecords(fR.id[selR], fR.seq[selR], fR.qual[selR])
    else
        length(fF) == length(fR) || error("Mismatched forward and reverse sequence files: $(length(fF)), $(length(fR)).")
    end
    n = length(fF)
    if orient_fwd !== nothing
        L = length(orient_fwd)
        for k in 1:n
            sf = String(copy(fF.seq[k]))
            if !(length(sf) >= L && sf[1:L] == orient_fwd)
                sr = String(copy(fR.seq[k]))
                if length(sr) >= L && sr[1:L] == orient_fwd
                    fF.seq[k], fR.seq[k] = fR.seq[k], fF.seq[k]
                    fF.qual[k], fR.qual[k] = fR.qual[k], fF.qual[k]
                    fF.id[k], fR.id[k] = fR.id[k], fF.id[k]
                else
                    fF.seq[k] = UInt8[]   # dropped below
                end
            end
        end
    end
    rF = Vector{Union{Nothing,UnitRange{Int}}}(nothing, n); rR = similar(rF)
    for k in 1:n
        isempty(fF.seq[k]) && orient_fwd !== nothing && continue
        a = _filter_read(fF.seq[k], fF.qual[k], FF, true)
        a === nothing && continue
        b = _filter_read(fR.seq[k], fR.qual[k], FR, true)
        b === nothing && continue
        rF[k] = a; rR[k] = b
    end
    idx = findall(!isnothing, rF)
    keep = trues(length(idx))
    sF = String[String(fF.seq[k][rF[k]]) for k in idx]; sR = String[String(fR.seq[k][rR[k]]) for k in idx]
    if ph[1] || ph[2]
        isphi = falses(length(idx))
        ph[1] && (isphi .|= isPhiX(sF))
        ph[2] && (isphi .|= isPhiX(sR))
        keep .&= .!isphi
    end
    if lc[1] > 0 || lc[2] > 0
        lowc = falses(length(idx))
        lc[1] > 0 && (lowc .|= seqComplexity(sF) .< lc[1])
        lc[2] > 0 && (lowc .|= seqComplexity(sR) .< lc[2])
        keep .&= .!lowc
    end
    idx = idx[keep]
    _write_kept(foutF, fF, idx, rF, compress)
    _write_kept(foutR, fR, idx, rR, compress)
    return nin, length(idx)
end

function _detect_id_field(id::AbstractString, id_sep)
    fields = split(id, id_sep)
    ncolon = [count(==(':'), f) for f in fields]
    count(==(6), ncolon) == 1 && maximum(ncolon) == 6 && return findfirst(==(6), ncolon)
    count(==(4), ncolon) == 1 && maximum(ncolon) == 4 && return findfirst(==(4), ncolon)
    error("Couldn't automatically detect the sequence identifier field in the fastq id string.")
end

function _write_kept(fout, fq::FastqRecords, idx, ranges, compress)
    isfile(fout) && rm(fout)
    isempty(idx) && return
    w = FastqWriter(fout; compress=compress)
    try
        for k in idx
            r = ranges[k]
            write_record(w, fq.id[k], view(fq.seq[k], r), view(fq.qual[k], r))
        end
    finally
        close(w)
    end
end

"""
    fastqFilter(fn, fout; kwargs...) -> (reads_in, reads_out)

Filter and trim one FASTQ file; the single-file form of
[`filterAndTrim`](@ref), taking the same keywords.
"""
function fastqFilter(fn::AbstractString, fout::AbstractString; kwargs...)
    t = filterAndTrim(fn, fout; progress=false, kwargs...)
    return (t.reads_in[1], t.reads_out[1])
end

"""
    fastqPairedFilter(fn, fout; kwargs...) -> (reads_in, reads_out)

Filter and trim one pair of FASTQ files, `fn = (fwd, rev)` into `fout =
(filt_fwd, filt_rev)`; the single-pair form of [`filterAndTrim`](@ref).
Per-read keywords may be given as `(fwd, rev)` pairs.
"""
function fastqPairedFilter(fn, fout; kwargs...)
    length(fn) == 2 && length(fout) == 2 || error("Two paired input and two output files are required.")
    t = filterAndTrim(fn[1], fout[1], fn[2], fout[2]; progress=false, kwargs...)
    return (t.reads_in[1], t.reads_out[1])
end
