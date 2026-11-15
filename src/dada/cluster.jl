# The DADA2 divisive partitioning algorithm.
#
# Faithful port of dada2's containers.cpp / cluster.cpp / pval.cpp. Iteration
# orders, tie-breaking and floating-point evaluation order are preserved so that
# the partition is identical to R's.

struct Comparison
    i::Int32          # cluster
    index::Int32      # raw
    lambda::Float64
    hamming::Int32
end
const NO_COMP = Comparison(0, 0, 0.0, 0)

mutable struct Raw
    seq::Vector{UInt8}
    qual::Vector{UInt8}       # rounded qualities (empty when absent)
    prior::Bool
    length::Int
    reads::Int
    index::Int
    p::Float64
    E_minmax::Float64
    comp::Comparison
    lock::Bool
    correct::Bool
end

mutable struct Bi
    center::Int               # raw index (0 = none)
    raws::Vector{Int}
    reads::Int
    i::Int
    update_e::Bool
    check_locks::Bool
    birth_type::Char
    birth_from::Int
    birth_pval::Float64
    birth_fold::Float64
    birth_e::Float64
    birth_comp::Comparison
    comp::Vector{Comparison}
end
Bi() = Bi(0, Int[], 0, 0, true, true, 'I', 0, 0.0, 1.0, 0.0, NO_COMP, Comparison[])

"Algorithm parameters for a single partitioning run (C-level options)."
Base.@kwdef struct ClusterParams
    match::Int = 5
    mismatch::Int = -4
    gap_p::Int = -8
    homo_gap_p::Int = -8
    use_kmers::Bool = true
    kdist_cutoff::Float64 = 0.42
    band::Int = 16
    omegaA::Float64 = 1e-40
    omegaP::Float64 = 1e-4
    detect_singletons::Bool = false
    max_clust::Int = 0
    min_fold::Float64 = 1.0
    min_hamming::Int = 1
    min_abund::Int = 1
    use_quals::Bool = true
    vectorized_alignment::Bool = true
    gapless::Bool = true
    greedy::Bool = true
end

mutable struct B
    nraw::Int
    reads::Int
    raw::Vector{Raw}
    bi::Vector{Bi}
    kmer::KmerTable           # k-mer counts of the raws
    kord::Matrix{UInt16}      # maxlen x nraw
    ident::Matrix{Float64}    # maxlen x nraw: P(no error) at each position (see identity_table!)
end

nclust(b::B) = length(b.bi)

# ---------------------------------------------------------------- containers

function bi_add_raw!(bi::Bi, b::B, index::Int)
    push!(bi.raws, index)
    bi.reads += b.raw[index].reads
    bi.update_e = true
    return length(bi.raws)
end

function bi_pop_raw!(bi::Bi, b::B, r::Int)
    index = bi.raws[r]
    bi.raws[r] = bi.raws[end]
    pop!(bi.raws)
    bi.reads -= b.raw[index].reads
    bi.update_e = true
    return index
end

function bi_assign_center!(bi::Bi, b::B)
    bi.center = 0
    max_reads = 0
    @inbounds for index in bi.raws
        raw = b.raw[index]
        raw.lock = false
        if raw.reads > max_reads
            bi.center = index
            max_reads = raw.reads
        end
    end
    bi.check_locks = true
    return bi
end

function b_add_bi!(b::B, bi::Bi)
    push!(b.bi, bi)
    bi.i = length(b.bi)
    return bi.i
end

function b_new(raws::Vector{Raw}, maxlen::Int)
    nraw = length(raws)
    kord = zeros(UInt16, max(maxlen, 1), nraw)
    reads = 0
    for (index, raw) in enumerate(raws)
        raw.index = index
        reads += raw.reads
        assign_kmer_order!(view(kord, :, index), raw.seq)
    end
    kmer = KmerTable([raw.seq for raw in raws])
    b = B(nraw, reads, raws, Bi[], kmer, kord, zeros(0, 0))
    bi0 = Bi()
    bi0.birth_e = reads
    b_add_bi!(b, bi0)
    for index in 1:nraw
        bi_add_raw!(bi0, b, index)
    end
    bi_assign_center!(bi0, b)
    return b
end

# ---------------------------------------------------------------- alignment

"""
Align raw `j` against raw `i` (the reference) using the kmer and gapless
screens of `raw_align`. Returns `true` and fills `sub` on success, `false`
if the pair was screened out by the kmer distance.
"""
function sub_new!(sub::Sub, buf::AlignBuffer, b::B, i::Int, j::Int, P::ClusterParams,
                  use_kmers::Bool, kdist_cutoff::Float64)
    r0 = b.raw[i]; r1 = b.raw[j]
    kdist = 0.0
    kodist = -1.0
    if use_kmers
        kdist = kmer_dist(b.kmer, i, r0.length, j, r1.length)
        if P.gapless
            kodist = kord_dist(view(b.kord, :, i), r0.length, view(b.kord, :, j), r1.length)
        end
        kdist > kdist_cutoff && return false
    end
    if P.band == 0 || (P.gapless && kodist == kdist)
        nwalign_gapless!(buf, r0.seq, r1.seq)
    elseif P.vectorized_alignment
        nwalign_vectorized!(buf, r0.seq, r1.seq, P.match, P.mismatch, P.gap_p, 0, P.band)
    elseif P.homo_gap_p != P.gap_p && P.homo_gap_p <= 0
        nwalign_endsfree_homo!(buf, r0.seq, r1.seq, P.match, P.mismatch, P.gap_p, P.homo_gap_p, P.band)
    else
        nwalign_endsfree!(buf, r0.seq, r1.seq, P.match, P.mismatch, P.gap_p, P.band)
    end
    al2subs!(sub, buf)
    return true
end

# ---------------------------------------------------------------- lambda

# lambda: the probability that a raw is produced from a reference by a given
# set of substitutions, the product over the raw's positions of
# err[transition, quality]. At unsubstituted positions the factor depends
# only on the raw and the error model, so it is tabulated once per run.

"Fill `b.ident[pos, r]` = err[nt -> nt, quality] for every raw `r`."
function identity_table!(b::B, err::Matrix{Float64}, use_quals::Bool)
    maxlen = size(b.kord, 1)
    b.ident = Matrix{Float64}(undef, maxlen, b.nraw)
    ncol = size(err, 2)
    @inbounds for (r, raw) in enumerate(b.raw)
        for pos in 1:raw.length
            nti = Int(raw.seq[pos]) - 1
            (0 <= nti <= 3) || error("Non-ACGT sequences in compute_lambda.")
            q = use_quals ? Int(raw.qual[pos]) : 0
            q < ncol || error("Rounded quality exceeded range of err lookup table.")
            b.ident[pos, r] = err[5nti+1, q+1]
        end
    end
    return b
end

@inline _sub_row(sub::Sub, s) = 4 * (Int(sub.nt0[s]) - 1) + Int(sub.nt1[s])

# Product over positions in sequence order (as R does): the tabulated
# identity factor, except at the `n` substituted positions `pos(s)` (which
# increase with s), where the factor is err[row(s), quality].
@inline function _lambda(b::B, r::Int, err, use_quals::Bool, n::Int, pos, row)
    raw = b.raw[r]; qual = raw.qual; ident = b.ident
    lambda = 1.0
    s = 1
    next = n >= 1 ? Int(pos(1)) : 0
    @inbounds for p in 1:raw.length
        if p == next
            lambda *= err[row(s), use_quals ? Int(qual[p]) + 1 : 1]
            s += 1
            next = s <= n ? Int(pos(s)) : 0
        else
            lambda *= ident[p, r]
        end
    end
    s == n + 1 || error("compute_lambda: bad substitution position")
    (0.0 <= lambda <= 1.0) || error("Bad lambda.")
    return lambda
end

"""
    compute_lambda(b, r, sub, err, use_quals)

Probability that raw `r` is produced from the reference by the substitutions
in `sub` under error model `err` (16 x nq, rows A2A,A2C,...,T2T; column Q+1).
"""
compute_lambda(b::B, r::Int, sub::Sub, err, use_quals::Bool) =
    _lambda(b, r, err, use_quals, sub.nsubs, s -> sub.map[sub.pos[s]], s -> _sub_row(sub, s))

# Same, from a compact list of (raw position, transition row) substitutions.
compute_lambda(b::B, r::Int, pos1::AbstractVector{UInt16}, rows::AbstractVector{UInt8}, err, use_quals::Bool) =
    _lambda(b, r, err, use_quals, length(pos1), s -> pos1[s], s -> Int(rows[s]))

# ---------------------------------------------------------------- compare

# Raws aligned together in one batch (Int16 lanes; an AVX-512 register holds 32).
const LANES = 32

"""
Substitutions of raws aligned to one cluster center. They depend only on the
two sequences, not on the error model, so across the rounds of
`learnErrors` each (center, raw) pair is aligned once.
"""
mutable struct CenterAlignments
    has::BitVector
    start::Vector{Int32}
    nsub::Vector{Int32}
    pos1::Vector{UInt16}      # substituted positions in the raw
    row::Vector{UInt8}        # transition (error matrix row) at each
end
CenterAlignments(nraw) = CenterAlignments(falses(nraw), zeros(Int32, nraw), zeros(Int32, nraw), UInt16[], UInt8[])

"Alignment cache of one sample: center raw index => its alignments."
const AlignCache = Dict{Int,CenterAlignments}

# Per-thread working memory for comparisons.
struct Workspace
    buf::AlignBuffer
    sub::Sub
    vec::VecBuffer{LANES}
    queries::Vector{Vector{UInt8}}
    pos1::Vector{UInt16}      # substitutions found in this thread's batches
    row::Vector{UInt8}
end
Workspace() = Workspace(AlignBuffer(), Sub(), VecBuffer{LANES}(),
                              Vector{UInt8}[], UInt16[], UInt8[])

# Scratch for b_compare!, one entry per raw.
struct CompareScratch
    lambda::Vector{Float64}
    ham::Vector{Int32}
    queued::BitVector         # needs a full alignment
    thread::Vector{Int32}     # where its substitutions were recorded (for the cache)
    start::Vector{Int32}
end
CompareScratch(n) = CompareScratch(zeros(n), zeros(Int32, n), falses(n), zeros(Int32, n), zeros(Int32, n))

"""
First pass over one raw: apply the greedy and k-mer screens, and handle the
cheap cases (gapless alignment, cached alignment, non-vectorized aligners).
Returns `(lambda, hamming, queued)`; queued raws need a vectorized alignment.
"""
function _screen(b::B, center, creads, index, err, P::ClusterParams, kdist_cutoff, ws::Workspace, cached)
    raw = b.raw[index]
    P.greedy && (raw.reads > creads || raw.lock) && return (0.0, Int32(-1), false)
    r0 = b.raw[center]
    kdist = 0.0; kodist = -1.0
    if P.use_kmers
        kdist = kmer_dist(b.kmer, center, r0.length, index, raw.length)
        P.gapless && (kodist = kord_dist(view(b.kord, :, center), r0.length, view(b.kord, :, index), raw.length))
        kdist > kdist_cutoff && return (0.0, Int32(-1), false)
    end
    if P.band == 0 || (P.gapless && kodist == kdist)
        nwalign_gapless!(ws.buf, r0.seq, raw.seq)
    elseif cached !== nothing && cached.has[index]
        r = cached.start[index] .+ (0:cached.nsub[index]-1)
        λ = compute_lambda(b, index, view(cached.pos1, r), view(cached.row, r), err, P.use_quals)
        return (λ, cached.nsub[index], false)
    elseif P.vectorized_alignment
        return (0.0, Int32(-1), true)
    elseif P.homo_gap_p != P.gap_p && P.homo_gap_p <= 0
        nwalign_endsfree_homo!(ws.buf, r0.seq, raw.seq, P.match, P.mismatch, P.gap_p, P.homo_gap_p, P.band)
    else
        nwalign_endsfree!(ws.buf, r0.seq, raw.seq, P.match, P.mismatch, P.gap_p, P.band)
    end
    al2subs!(ws.sub, ws.buf)
    return (compute_lambda(b, index, ws.sub, err, P.use_quals), Int32(ws.sub.nsubs), false)
end

# (a function barrier, so the loop is compiled for the concrete type of `cached`)
function _screen_range!(cs::CompareScratch, range, b, center, creads, err, P, kdist_cutoff, ws, cached)
    @inbounds for index in range
        l, h, q = _screen(b, center, creads, index, err, P, kdist_cutoff, ws, cached)
        cs.lambda[index] = l; cs.ham[index] = h; cs.queued[index] = q
    end
end

# Align one batch of equal-length raws to the center, `LANES` at a time.
function _align_batch!(b::B, center, batch, err, P::ClusterParams, ws::Workspace, cs::CompareScratch,
                       t, record::Bool)
    empty!(ws.queries)
    for j in batch
        push!(ws.queries, b.raw[j].seq)
    end
    load!(ws.vec, b.raw[center].seq, ws.queries)
    fill_dp!(ws.vec, P.match, P.mismatch, P.gap_p, 0, P.band)
    for (lane, j) in enumerate(batch)
        traceback!(ws.buf, ws.vec, lane)
        al2subs!(ws.sub, ws.buf)
        cs.lambda[j] = compute_lambda(b, j, ws.sub, err, P.use_quals)
        cs.ham[j] = ws.sub.nsubs
        if record
            cs.thread[j] = t; cs.start[j] = length(ws.pos1) + 1
            for s in 1:ws.sub.nsubs
                push!(ws.pos1, ws.sub.map[ws.sub.pos[s]]); push!(ws.row, _sub_row(ws.sub, s))
            end
        end
    end
end

"""
    member_subs(b, P, ws) -> Vector{Union{Nothing,Sub}}

Alignment of every raw to the center of its cluster (no k-mer screen), as
used for the output transition counts. Members are aligned in batches.
"""
function member_subs(b::B, P::ClusterParams, ws::Workspace)
    subs = Vector{Union{Nothing,Sub}}(nothing, b.nraw)
    batched = P.vectorized_alignment && P.band != 0
    for bi in b.bi
        if !batched
            for index in bi.raws
                ok = sub_new!(ws.sub, ws.buf, b, bi.center, index, P, false, 1.0)
                subs[index] = ok ? copy(ws.sub) : nothing
            end
            continue
        end
        members = falses(b.nraw)
        members[bi.raws] .= true
        q, batches = _batches(b, members)
        for r in batches
            empty!(ws.queries)
            for j in view(q, r)
                push!(ws.queries, b.raw[j].seq)
            end
            load!(ws.vec, b.raw[bi.center].seq, ws.queries)
            fill_dp!(ws.vec, P.match, P.mismatch, P.gap_p, 0, P.band)
            for (lane, j) in enumerate(view(q, r))
                traceback!(ws.buf, ws.vec, lane)
                subs[j] = copy(al2subs!(ws.sub, ws.buf))
            end
        end
    end
    return subs
end

# Split the queued raws into batches of at most LANES raws of equal length.
function _batches(b::B, queued::BitVector)
    q = findall(queued)
    sort!(q; by=j -> b.raw[j].length, alg=Base.Sort.DEFAULT_STABLE)
    batches = UnitRange{Int}[]
    lo = 1
    while lo <= length(q)
        len = b.raw[q[lo]].length
        hi = lo
        while hi < length(q) && hi - lo + 1 < LANES && b.raw[q[hi+1]].length == len
            hi += 1
        end
        push!(batches, lo:hi)
        lo = hi + 1
    end
    return q, batches
end

"""
Compare the center of cluster `i` with every raw and store the comparisons
that could recruit the raw (port of `b_compare`). Lambdas are computed in
parallel and stored serially in raw order, so the outcome does not depend on
the number of threads or on batching.
"""
function b_compare!(b::B, i::Int, err::Matrix{Float64}, P::ClusterParams, kdist_cutoff::Float64,
                    wss::Vector{Workspace}, cs::CompareScratch, cache::Union{Nothing,AlignCache}=nothing)
    bi = b.bi[i]
    center = bi.center
    center_reads = b.raw[center].reads
    nraw = b.nraw
    cached = cache === nothing ? nothing : get(cache, center, nothing)
    nt = length(wss)
    parallel = nt > 1 && nraw >= 64

    # 1. screens and cheap cases
    if parallel
        Threads.@threads :static for t in 1:nt
            _screen_range!(cs, div((t - 1) * nraw, nt)+1:div(t * nraw, nt), b, center, center_reads,
                           err, P, kdist_cutoff, wss[t], cached)
        end
    else
        _screen_range!(cs, 1:nraw, b, center, center_reads, err, P, kdist_cutoff, wss[1], cached)
    end

    # 2. batched vectorized alignments
    q, batches = _batches(b, cs.queued)
    record = cache !== nothing
    if parallel && length(batches) > 1
        Threads.@threads :static for t in 1:nt
            for k in t:nt:length(batches)
                _align_batch!(b, center, view(q, batches[k]), err, P, wss[t], cs, t, record)
            end
        end
    else
        for r in batches
            _align_batch!(b, center, view(q, r), err, P, wss[1], cs, 1, record)
        end
    end
    if record && !isempty(q)
        entry = get!(() -> CenterAlignments(nraw), cache, center)
        for j in q
            ws = wss[cs.thread[j]]
            n = Int(cs.ham[j])
            entry.has[j] = true
            entry.start[j] = length(entry.pos1) + 1
            entry.nsub[j] = n
            r = cs.start[j] .+ (0:n-1)
            append!(entry.pos1, view(ws.pos1, r)); append!(entry.row, view(ws.row, r))
        end
        for ws in wss
            empty!(ws.pos1); empty!(ws.row)
        end
    end

    # 3. store the comparisons, in raw order
    lambdas = cs.lambda; hams = cs.ham
    breads = Float64(b.reads)
    creads = Float64(center_reads)
    @inbounds for index in 1:nraw
        raw = b.raw[index]
        lambda = lambdas[index]
        if lambda * breads > raw.E_minmax
            if lambda * creads > raw.E_minmax
                raw.E_minmax = lambda * creads
            end
            comp = Comparison(i, index, lambda, hams[index])
            push!(bi.comp, comp)
            if i == 1 || index == center
                raw.comp = comp
            end
        end
    end
    return b
end

# ---------------------------------------------------------------- p-values

function get_pA(raw::Raw, bi::Bi, detect_singletons::Bool)
    lambda = raw.comp.lambda
    hamming = raw.comp.hamming
    if raw.reads == 1 && !raw.prior && !detect_singletons
        return 1.0
    elseif hamming == 0
        return 1.0
    elseif lambda == 0.0
        return 0.0
    else
        return calc_pA(raw.reads, lambda * bi.reads, raw.prior || detect_singletons)
    end
end

function b_p_update!(b::B, greedy::Bool, detect_singletons::Bool)
    for bi in b.bi
        if bi.update_e
            for index in bi.raws
                raw = b.raw[index]
                raw.p = get_pA(raw, bi, detect_singletons)
            end
            bi.update_e = false
        end
        if greedy && bi.check_locks
            creads = b.raw[bi.center].reads
            for index in bi.raws
                raw = b.raw[index]
                E_center = creads * raw.comp.lambda
                E_center > raw.reads && (raw.lock = true)
                index == bi.center && (raw.lock = true)
            end
            bi.check_locks = false
        end
    end
    return b
end

# ---------------------------------------------------------------- shuffle

function b_shuffle!(b::B, compmax::Vector{Comparison}, emax::Vector{Float64})
    shuffled = false
    bi1 = b.bi[1]
    r1 = Float64(bi1.reads)
    @inbounds for index in 1:b.nraw
        c = bi1.comp[index]
        compmax[index] = c
        emax[index] = c.lambda * r1
    end
    @inbounds for i in 2:nclust(b)
        bi = b.bi[i]
        reads = Float64(bi.reads)
        for c in bi.comp
            e = c.lambda * reads
            if e > emax[c.index]
                compmax[c.index] = c
                emax[c.index] = e
            end
        end
    end
    for i in 1:nclust(b)
        bi = b.bi[i]
        r = length(bi.raws)
        while r >= 1
            index = bi.raws[r]
            best = compmax[index]
            if best.i != i && index != bi.center
                bi_pop_raw!(bi, b, r)
                bi_add_raw!(b.bi[best.i], b, index)
                b.raw[index].comp = best
                shuffled = true
            end
            r -= 1
        end
    end
    return shuffled
end

# ---------------------------------------------------------------- bud

function b_bud!(b::B, P::ClusterParams)
    mini = 0; minr = 0; mini_prior = 0; minr_prior = 0
    minraw = b.raw[b.bi[1].center]
    minraw_prior = minraw
    for i in 1:nclust(b)
        bi = b.bi[i]
        bireads = bi.reads
        for r in 2:length(bi.raws)
            raw = b.raw[bi.raws[r]]
            raw.reads < P.min_abund && continue
            hamming = raw.comp.hamming
            lambda = raw.comp.lambda
            if hamming >= P.min_hamming
                if P.min_fold <= 1 || Float64(raw.reads) >= P.min_fold * lambda * bireads
                    if raw.p < minraw.p || (raw.p == minraw.p && raw.reads > minraw.reads)
                        mini = i; minr = r; minraw = raw
                    end
                    if raw.prior && (raw.p < minraw_prior.p || (raw.p == minraw_prior.p && raw.reads > minraw_prior.reads))
                        mini_prior = i; minr_prior = r; minraw_prior = raw
                    end
                end
            end
        end
    end
    pA = minraw.p * b.nraw
    pP = minraw_prior.p
    if pA < P.omegaA && mini >= 1
        return _bud_from!(b, mini, minr, minraw, 'A', pA)
    elseif pP < P.omegaP && mini_prior >= 1
        return _bud_from!(b, mini_prior, minr_prior, minraw_prior, 'P', pP)
    end
    return 0
end

function _bud_from!(b::B, from::Int, r::Int, minraw::Raw, type::Char, pval::Float64)
    expected = minraw.comp.lambda * b.bi[from].reads
    index = bi_pop_raw!(b.bi[from], b, r)
    bi = Bi()
    i = b_add_bi!(b, bi)
    bi.birth_type = type
    bi.birth_from = from
    bi.birth_pval = pval
    bi.birth_fold = b.raw[index].reads / expected
    bi.birth_e = expected
    bi.birth_comp = minraw.comp
    bi_add_raw!(bi, b, index)
    bi_assign_center!(bi, b)
    return i
end

# ---------------------------------------------------------------- driver

const MAX_SHUFFLE = 10

"""
    run_dada(raws, err, P; maxlen, nthreads, onbud=nothing)

Partition `raws` under error model `err`. Returns the clustering `B`.
`onbud(nclust)` is called after each new cluster (for progress display).
"""
function run_dada(raws::Vector{Raw}, err::Matrix{Float64}, P::ClusterParams;
                  maxlen::Int=maximum(r -> r.length, raws), nthreads::Int=1, onbud=nothing,
                  cache::Union{Nothing,AlignCache}=nothing)
    b = b_new(raws, maxlen)
    identity_table!(b, err, P.use_quals)
    wss = [Workspace() for _ in 1:max(nthreads, 1)]
    cs = CompareScratch(b.nraw)
    # everyone is aligned to the initial cluster (no kmer screen)
    b_compare!(b, 1, err, P, 1.0, wss, cs, cache)
    b_p_update!(b, P.greedy, P.detect_singletons)
    max_clust = P.max_clust < 1 ? b.nraw : P.max_clust
    compmax = Vector{Comparison}(undef, b.nraw)
    emax = Vector{Float64}(undef, b.nraw)
    while nclust(b) < max_clust
        newi = b_bud!(b, P)
        newi == 0 && break
        b_compare!(b, newi, err, P, P.kdist_cutoff, wss, cs, cache)
        nshuffle = 0
        while true
            shuffled = b_shuffle!(b, compmax, emax)
            nshuffle += 1
            (shuffled && nshuffle < MAX_SHUFFLE) || break
        end
        b_p_update!(b, P.greedy, P.detect_singletons)
        onbud === nothing || onbud(nclust(b))
    end
    return b
end
