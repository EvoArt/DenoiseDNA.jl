# Primer removal (port of dada2::removePrimers) and the parts of Biostrings'
# pattern matching it relies on: `vmatchPattern` (mismatches only, matches
# may run off either end of the read) and `matchPattern(with.indels=TRUE)`
# (Biostrings' "best local matches" under edit distance).

# Biostrings' internal DNA codes: one bit per base, so IUPAC codes are unions.
const _DNA_BITS = let t = zeros(UInt8, 256)
    for (c, v) in zip("ACMGRSVTWYHKDBN-+.", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 32, 64))
        t[UInt8(c)+1] = v
        t[UInt8(lowercase(c))+1] = v
    end
    t
end
_dnabits(s::AbstractString) = _dnabits(codeunits(s))
_dnabits(v::AbstractVector{UInt8}) = [_DNA_BITS[c+1] for c in v]

# fixed: letters match iff equal; otherwise iff their base sets intersect.
@inline _lmatch(x::UInt8, y::UInt8, fixed::Bool) = fixed ? x == y : (x & y) != 0x00

# Mismatches of pattern P with its first letter at 0-based offset `shift` in S;
# letters off the ends of S count as mismatches. Stops once past `maxmis`.
function _nmismatch_at(P, S, shift, maxmis, fixed)
    nmis = 0
    @inbounds for i in eachindex(P)
        j = shift + i
        1 <= j <= length(S) && _lmatch(P[i], S[j], fixed) && continue
        nmis += 1
        nmis > maxmis && break
    end
    return nmis
end

"""
    _match_pattern(P, S, maxmis, indels, fixed) -> Vector{UnitRange{Int}}

Matches of pattern P in subject S (both Biostrings codes), as Biostrings'
`matchPattern(P, S, max.mismatch=maxmis, with.indels=indels, fixed=fixed)`.
Without indels, ranges may extend past either end of S.
"""
function _match_pattern(P, S, maxmis, indels, fixed)
    hits = UnitRange{Int}[]
    lp = length(P); ls = length(S)
    lp > 0 || throw(ArgumentError("empty pattern"))
    maxmis < lp - ls && return hits
    if indels && lp > maxmis
        _match_indels!(hits, P, S, maxmis, fixed)
    else
        lo = lp <= maxmis ? 1 - lp : -maxmis
        for shift in lo:(ls - lo - lp)
            _nmismatch_at(P, S, shift, maxmis, fixed) <= maxmis && push!(hits, shift+1:shift+lp)
        end
    end
    return hits
end

# Biostrings' _match_pattern_indels: report only "best local matches", holding
# a provisory match until a later one ends further right.
function _match_indels!(hits, P, S, maxmis, fixed)
    lp = length(P)
    first_offset = fill(-1, 256)          # first (0-based) letter of P matching each code
    for y in 0x00:0xff
        n = findfirst(p -> _lmatch(p, y, fixed), P)
        n === nothing || (first_offset[y+1] = n - 1)
    end
    prov_start = prov_stop = 0; prov_nedit = -1
    rows = (zeros(Int, 2maxmis + 1), zeros(Int, 2maxmis + 1))
    for j0 in 0:length(S)-1
        i0 = first_offset[S[j0+1]+1]
        i0 < 0 && continue
        P1 = view(P, i0+2:lp)
        m1 = maxmis - i0
        m1 < 0 && continue
        if m1 == 0
            nedit = _nmismatch_at(P1, S, j0 + 1, 0, fixed); width = length(P1)
        else
            nedit, width = _nedit_for_Ploffset(P1, S, j0 + 1, m1, fixed, rows)
        end
        nedit <= m1 || continue
        start = j0 + 1; stop = start + width; nedit += i0
        if prov_nedit != -1
            if stop > prov_stop
                push!(hits, prov_start:prov_stop)
            elseif nedit > prov_nedit
                continue
            end
        end
        prov_start, prov_stop, prov_nedit = start, stop, nedit
    end
    prov_nedit != -1 && push!(hits, prov_start:prov_stop)
    return hits
end

# Biostrings' _nedit_for_Ploffset: smallest edit distance between P and a
# substring of S starting at 0-based offset `off`, with a band of half-width
# `maxnedit` and early bailout. Returns (nedit, width of the shortest such substring).
function _nedit_for_Ploffset(P, S, off, maxnedit, fixed, rows)
    lp = length(P)
    lp == 0 && return 0, 0
    bailout = maxnedit + 1
    maxnedit = min(maxnedit, lp)
    rowlen = 2maxnedit + 1
    prev, curr = rows
    propagate!(curr, prev, B, Si, Pc) = _propagate_nedit!(curr, prev, B, S, Si, Pc, fixed, rowlen)
    minSi = off
    for B in maxnedit:rowlen-1                       # stage 0
        curr[B+1] = B - maxnedit
    end
    a = 1; Pi = 0
    while a < maxnedit                               # stage 1
        prev, curr = curr, prev
        B = maxnedit - a
        curr[B+1] = a
        for (k, BB) in enumerate(B+1:rowlen-1)
            propagate!(curr, prev, BB, minSi + k - 1, P[Pi+1])
        end
        a += 1; Pi += 1
    end
    prev, curr = curr, prev                          # stage 2
    curr[1] = a; minnedit = a; width = 0
    for (k, B) in enumerate(1:rowlen-1)
        Si = minSi + k - 1
        propagate!(curr, prev, B, Si, P[Pi+1])
        curr[B+1] < minnedit && (minnedit = curr[B+1]; width = Si - off + 1)
    end
    a += 1; Pi += 1
    while Pi < lp                                    # stage 3, with bailout
        prev, curr = curr, prev
        minnedit = a; width = 0
        for B in 0:rowlen-1
            Si = minSi + B
            propagate!(curr, prev, B, Si, P[Pi+1])
            curr[B+1] < minnedit && (minnedit = curr[B+1]; width = Si - off + 1)
        end
        minnedit >= bailout && break
        Pi += 1; a += 1; minSi += 1
    end
    return minnedit, width
end

# One cell of the banded edit-distance row (Biostrings' PROPAGATE_NEDIT); B and Si are 0-based.
@inline function _propagate_nedit!(curr, prev, B, S, Si, Pc, fixed, rowlen)
    @inbounds begin
        nedit = prev[B+1] + !(0 <= Si < length(S) && _lmatch(Pc, S[Si+1], fixed))
        B >= 1 && curr[B] + 1 < nedit && (nedit = curr[B] + 1)
        B + 1 < rowlen && prev[B+2] + 1 < nedit && (nedit = prev[B+2] + 1)
        curr[B+1] = nedit
    end
    return nothing
end

"""
    removePrimers(fn, fout, primer_fwd, primer_rev=nothing; max_mismatch=2,
                  allow_indels=false, trim_fwd=true, trim_rev=true, orient=true,
                  compress=true, verbose=false) -> DataFrame

Keep only reads that contain the forward primer (and, if given, the reverse
primer, in read orientation: e.g. `rc(R1492)`), trim the primers off and
write the reads to `fout`. With `orient`, reads in which only the reverse
complement holds the forward primer are turned around. IUPAC codes in a
primer match any base they stand for. Port of R's `removePrimers` (used for
PacBio full-length amplicons). Returns reads in and out per file.
"""
function removePrimers(fn, fout, primer_fwd::AbstractString, primer_rev=nothing; max_mismatch::Integer=2,
                       allow_indels::Bool=false, trim_fwd::Bool=true, trim_rev::Bool=true,
                       orient::Bool=true, compress::Bool=true, verbose::Bool=false)
    fn = fn isa AbstractString ? [fn] : collect(fn)
    fout = fout isa AbstractString ? [fout] : collect(fout)
    length(fn) == length(fout) || error("Every input file must have a corresponding output file.")
    allow_indels && @info "Primer matching with indels allowed is slower."
    all(isfile, fn) || error("Some input files do not exist.")
    fn = abspath.(fn); fout = abspath.(fout)
    allunique(fout) || error("All output files must be distinct.")
    any(in(fn), fout) && error("Output files must be distinct from the input files.")
    foreach(f -> mkpath(dirname(f)), fout)
    has_rev = primer_rev !== nothing
    Pf = _dnabits(primer_fwd); fixed_f = isACGT(primer_fwd)
    Pr = has_rev ? _dnabits(primer_rev) : UInt8[]; fixed_r = has_rev && isACGT(primer_rev)
    ins = zeros(Int, length(fn)); outs = zeros(Int, length(fn))
    multi_msg = true
    for i in eachindex(fn)
        fq = read_fastq(fn[i])
        ins[i] = n = length(fq)
        hits(P, fixed, codes) = [_match_pattern(P, c, max_mismatch, allow_indels, fixed) for c in codes]
        codes = _dnabits.(fq.seq)
        mf = hits(Pf, fixed_f, codes)
        mr = has_rev ? hits(Pr, fixed_r, codes) : nothing
        if orient
            seqs_rc = revcomp.(fq.seq)
            codes_rc = _dnabits.(seqs_rc)
            mf_rc = hits(Pf, fixed_f, codes_rc)
            mr_rc = has_rev ? hits(Pr, fixed_r, codes_rc) : nothing
        end
        hits_f = length.(mf); hits_r = has_rev ? length.(mr) : nothing
        # R gives up on a file before orienting if no read holds the primer as is
        if sum(hits_f) == 0 || (has_rev && sum(hits_r) == 0)
            verbose && println("Read in $n, output 0 filtered sequences.")
            continue
        end
        multi = any(>(1), hits_f) || (has_rev && any(>(1), hits_r)) ||
                (orient && (any(m -> length(m) > 1, mf_rc) || (has_rev && any(m -> length(m) > 1, mr_rc))))
        if multi && verbose && multi_msg
            @info "Multiple matches to the primer(s) in some sequences. Using the longest possible match."
            multi_msg = false
        end
        # forward primer: first match; reverse primer: last match
        pick_f(m) = isempty(m) ? nothing : first(m)
        pick_r(m) = isempty(m) ? nothing : last(m)
        seqs = fq.seq; quals = fq.qual
        f = pick_f.(mf); r = has_rev ? pick_r.(mr) : nothing
        if orient
            flip = (hits_f .== 0) .& (length.(mf_rc) .> 0)
            verbose && any(flip) && println(sum(flip), " sequences out of ", n, " are being reverse-complemented.")
            seqs = [flip[k] ? seqs_rc[k] : seqs[k] for k in 1:n]
            quals = [flip[k] ? reverse(quals[k]) : quals[k] for k in 1:n]
            f = [flip[k] ? pick_f(mf_rc[k]) : f[k] for k in 1:n]
            has_rev && (r = [flip[k] ? pick_r(mr_rc[k]) : r[k] for k in 1:n])
        end
        isfile(fout[i]) && rm(fout[i])
        w = FastqWriter(fout[i]; compress)
        try
            for k in 1:n
                (f[k] === nothing || (has_rev && r[k] === nothing)) && continue
                lo = trim_fwd ? last(f[k]) + 1 : 1
                hi = has_rev && trim_rev ? first(r[k]) - 1 : length(seqs[k])
                hi > lo || continue
                write_record(w, fq.id[k], view(seqs[k], lo:hi), view(quals[k], lo:hi))
                outs[i] += 1
            end
        finally
            close(w)
        end
        verbose && println("Read in $n, output $(outs[i]) filtered sequences.")
    end
    if all(==(0), outs)
        @warn "No reads passed the primer detection."
    elseif any(==(0), outs)
        @info "Some input samples had no reads pass the primer detection."
    end
    return DataFrame(file=basename.(fn), reads_in=ins, reads_out=outs)
end
