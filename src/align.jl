# Pairwise alignment used throughout DADA2.
#
# These are faithful ports of dada2's C++ aligners (src/nwalign_vectorized.cpp,
# src/nwalign_endsfree.cpp). Tie-breaking, banding and the anti-diagonal
# layout are reproduced exactly so that alignments (and therefore every
# downstream decision) match R's dada2 bit-for-bit.
#
# Sequences are integer coded: A=1, C=2, G=3, T=4, N=5. The gap glyph in an
# alignment is `GAP` ('-').

const GAP = UInt8('-')

"""
    AlignBuffer()

Reusable scratch space for the aligners. Results of the last alignment are in
`al0[1:len]` / `al1[1:len]` (forward order).
"""
mutable struct AlignBuffer
    p::Vector{Int16}          # traceback pointers
    d::Vector{Int16}          # scores (anti-diagonal aligner)
    diag::Vector{Int16}       # diagonal scores of one anti-diagonal
    d32::Vector{Int32}        # scores (row-major aligners)
    al0::Vector{UInt8}
    al1::Vector{UInt8}
    len::Int
end
AlignBuffer() = AlignBuffer(Int16[], Int16[], Int16[], Int32[], UInt8[], UInt8[], 0)

@inline function _grow!(v::Vector, n::Integer)
    length(v) < n && resize!(v, n)
    return v
end

# --------------------------------------------------------------------------
# Gapless alignment (used when the ordered-kmer screen says the sequences
# can be aligned without gaps).
# --------------------------------------------------------------------------
function nwalign_gapless!(buf::AlignBuffer, s1::AbstractVector{UInt8}, s2::AbstractVector{UInt8})
    len1 = length(s1); len2 = length(s2)
    len_al = max(len1, len2)
    _grow!(buf.al0, len_al); _grow!(buf.al1, len_al)
    al0 = buf.al0; al1 = buf.al1
    @inbounds for i in 1:len_al
        al0[i] = i <= len1 ? s1[i] : GAP
        al1[i] = i <= len2 ? s2[i] : GAP
    end
    buf.len = len_al
    return buf
end

# --------------------------------------------------------------------------
# Single anti-diagonal alignment (port of nwalign_vectorized2). Batches of
# equal-length queries go through the SIMD kernel in align_vec.jl, which
# follows exactly this recurrence in each lane; this scalar version is used
# for one-off alignments, where a one-lane SIMD batch would be slower.
# --------------------------------------------------------------------------
"""
    nwalign_vectorized!(buf, s1, s2, match, mismatch, gap_p, end_gap_p, band)

Align `s2` against `s1` with dada2's vectorized banded ends-free aligner;
the alignment is left in `buf`.
"""
function nwalign_vectorized!(buf::AlignBuffer, s1in::AbstractVector{UInt8}, s2in::AbstractVector{UInt8},
                             match::Integer, mismatch::Integer, gap_p::Integer, end_gap_p::Integer, band::Integer)
    match = Int16(match); mismatch = Int16(mismatch)
    gap_p = Int16(gap_p); end_gap_p = Int16(end_gap_p)
    s1 = s1in; s2 = s2in
    len1 = length(s1); len2 = length(s2)
    swap = false
    if len1 > len2
        s1, s2 = s2, s1
        len1, len2 = len2, len1
        swap = true
    end
    band = Int(band)
    band < 0 && (band = len2)

    start_col = 1 + (1 + min(band, len1)) ÷ 2
    ncol = 2 + start_col + min(len2 - len1 + band, len2) ÷ 2
    nrow = len1 + len2 + 1
    _grow!(buf.d, ncol * nrow); _grow!(buf.p, ncol * nrow); _grow!(buf.diag, ncol)
    d = buf.d; p = buf.p; diag_buf = buf.diag

    # 0-based (row, col) -> 1-based linear index
    @inline idx(r, c) = r * ncol + c + 1

    fill_val = Int16(Int(typemin(Int16)) - min(min(Int(mismatch), Int(gap_p)), min(Int(match), 0)))
    @inbounds for row in 0:nrow-1
        d[idx(row, 0)] = fill_val
        d[idx(row, 1)] = fill_val
        d[idx(row, ncol - 2)] = fill_val
        d[idx(row, ncol - 1)] = fill_val
    end

    d[idx(0, start_col)] = Int16(0)
    p[idx(0, start_col)] = Int16(0)

    # "left" column
    row = 1; col = start_col - 1; d_free = end_gap_p
    lim_left = min(band, len1)
    @inbounds while row < 1 + lim_left
        d[idx(row, col)] = d_free
        p[idx(row, col)] = Int16(3)
        row % 2 == 0 && (col -= 1)
        row += 1
        d_free += end_gap_p
    end
    # "top" row
    row = 1; col = start_col; d_free = end_gap_p
    lim_top = min(band + len2 - len1, len2)
    @inbounds while row < 1 + lim_top
        d[idx(row, col)] = d_free
        p[idx(row, col)] = Int16(2)
        row % 2 == 1 && (col += 1)
        row += 1
        d_free += end_gap_p
    end

    row = 2
    col_min = start_col; col_max = start_col
    i_max = 0; j_min = 0
    even = 1
    recalc_left = false; recalc_right = false
    @inbounds while row <= len1 + len2
        # diagonal scores for this anti-diagonal
        i = i_max; j = j_min
        rowm2 = (row - 2) * ncol + 1
        for c in col_min:col_max
            diag_buf[c+1] = d[rowm2 + c] + (s1[i+1] == s2[j+1] ? match : mismatch)
            i -= 1; j += 1
        end
        left0 = (row - 1) * ncol + col_min - even + 1   # 1-based linear index into d
        up0 = (row - 1) * ncol + col_min + 1 - even + 1
        out0 = row * ncol + col_min + 1
        n = col_max - col_min + 1
        if swap
            for t in 0:n-1
                left = d[left0 + t] + gap_p
                dg = diag_buf[col_min + t + 1]
                up = d[up0 + t] + gap_p
                entry = left >= up ? left : up
                pentry = left >= up ? Int16(2) : Int16(3)
                pentry = entry >= dg ? pentry : Int16(1)
                entry = entry >= dg ? entry : dg
                d[out0 + t] = entry
                p[out0 + t] = pentry
            end
        else
            for t in 0:n-1
                left = d[left0 + t] + gap_p
                dg = diag_buf[col_min + t + 1]
                up = d[up0 + t] + gap_p
                entry = up >= left ? up : left
                pentry = up >= left ? Int16(3) : Int16(2)
                pentry = entry >= dg ? pentry : Int16(1)
                entry = entry >= dg ? entry : dg
                d[out0 + t] = entry
                p[out0 + t] = pentry
            end
        end

        if row == min(band, len1)
            col_min -= 1; i_max += 1; j_min -= 1
        end
        if row == min(band + len2 - len1, len2)
            col_max += 1
        end

        # ends-free recalculation of the lower boundary
        if end_gap_p > gap_p
            if recalc_left
                dfree = d[left0] + end_gap_p
                k = idx(row, col_min)
                if dfree > d[k]
                    d[k] = dfree; p[k] = Int16(2)
                elseif !swap && dfree == d[k] && p[k] == 1
                    p[k] = Int16(2)
                elseif swap && dfree == d[k] && p[k] != 2
                    p[k] = Int16(2)
                end
            end
            i_max == len1 - 1 && (recalc_left = true)
            if recalc_right
                dfree = d[up0 + (col_max - col_min)] + end_gap_p
                k = idx(row, col_max)
                if dfree > d[k]
                    d[k] = dfree; p[k] = Int16(3)
                elseif !swap && dfree == d[k] && p[k] != 3
                    p[k] = Int16(3)
                elseif swap && dfree == d[k] && p[k] == 1
                    p[k] = Int16(3)
                end
            end
            _usz((row + 1) ÷ 2 + col_max - start_col) == len2 && (recalc_right = true)
        end

        # update the indices
        if row < band && row < len1
            even == 1 && (col_min -= 1)
            i_max += 1
        elseif i_max < len1 - 1
            if band % 2 == 0
                if even == 1
                    j_min += 1
                else
                    i_max += 1
                end
            else
                if even == 1
                    col_min -= 1; i_max += 1
                else
                    col_min += 1; j_min += 1
                end
            end
        else
            if even == 0
                col_min += 1
            end
            j_min += 1
        end

        if row < min(band + len2 - len1, len2)
            even == 0 && (col_max += 1)
        elseif _usz((row + 1) ÷ 2 + col_max - start_col) < len2
            if (band + len2 - len1) % 2 == 0
                if even == 1
                    col_max -= 1
                else
                    col_max += 1
                end
            end
        else
            even == 1 && (col_max -= 1)
        end

        row += 1
        even = 1 - even
    end

    # traceback
    _grow!(buf.al0, nrow); _grow!(buf.al1, nrow)
    a0 = swap ? buf.al1 : buf.al0   # a0 holds the (possibly swapped) s1 column
    a1 = swap ? buf.al0 : buf.al1
    len_al = 0
    i = len1; j = len2
    @inbounds while i > 0 || j > 0
        pij = p[idx(i + j, (2 * start_col + j - i) ÷ 2)]
        len_al += 1
        if pij == 1
            i -= 1; j -= 1
            a0[len_al] = s1[i+1]; a1[len_al] = s2[j+1]
        elseif pij == 2
            j -= 1
            a0[len_al] = GAP; a1[len_al] = s2[j+1]
        elseif pij == 3
            i -= 1
            a0[len_al] = s1[i+1]; a1[len_al] = GAP
        else
            error("N-W alignment out of range (len1=$len1, len2=$len2, i=$i, j=$j, p=$pij)")
        end
    end
    reverse!(view(buf.al0, 1:len_al)); reverse!(view(buf.al1, 1:len_al))
    buf.len = len_al
    return buf
end

# --------------------------------------------------------------------------
# Standard banded ends-free Needleman-Wunsch (port of nwalign_endsfree).
# Used by mergePairs (via nwalign) and when VECTORIZED_ALIGNMENT=FALSE.
# band < 0 disables banding.
# --------------------------------------------------------------------------
function nwalign_endsfree!(buf::AlignBuffer, s1::AbstractVector{UInt8}, s2::AbstractVector{UInt8},
                           match::Integer, mismatch::Integer, gap_p::Integer, band::Integer;
                           endsfree::Bool=true)
    len1 = length(s1); len2 = length(s2)
    ncol = len2 + 1; nrow = len1 + 1
    _grow!(buf.d32, ncol * nrow); _grow!(buf.p, ncol * nrow)
    d = buf.d32; p = buf.p
    @inline idx(i, j) = i * ncol + j + 1
    m32 = Int32(match); mm32 = Int32(mismatch); g32 = Int32(gap_p)

    if endsfree
        @inbounds for i in 0:len1
            d[idx(i, 0)] = 0; p[idx(i, 0)] = 3
        end
        @inbounds for j in 0:len2
            d[idx(0, j)] = 0; p[idx(0, j)] = 2
        end
    else
        d[1] = 0; p[1] = 0
        @inbounds for i in 1:len1
            d[idx(i, 0)] = d[idx(i - 1, 0)] + g32; p[idx(i, 0)] = 3
        end
        @inbounds for j in 1:len2
            d[idx(0, j)] = d[idx(0, j - 1)] + g32; p[idx(0, j)] = 2
        end
    end

    if len2 > len1
        lband = band; rband = band + len2 - len1
    elseif len1 > len2
        lband = band + len1 - len2; rband = band
    else
        lband = band; rband = band
    end
    if band >= 0 && (band < len1 || band < len2)
        @inbounds for i in 0:len1
            i - lband - 1 >= 0 && (d[idx(i, i - lband - 1)] = -9999)
            i + rband + 1 <= len2 && (d[idx(i, i + rband + 1)] = -9999)
        end
    end

    @inbounds for i in 1:len1
        if band >= 0
            l = max(i - lband, 1); r = min(i + rband, len2)
        else
            l = 1; r = len2
        end
        a = s1[i]
        for j in l:r
            left = d[idx(i, j - 1)]
            (!endsfree || i != len1) && (left += g32)
            up = d[idx(i - 1, j)]
            (!endsfree || j != len2) && (up += g32)
            dg = d[idx(i - 1, j - 1)] + (a == s2[j] ? m32 : mm32)
            k = idx(i, j)
            if up >= dg && up >= left
                d[k] = up; p[k] = 3
            elseif left >= dg
                d[k] = left; p[k] = 2
            else
                d[k] = dg; p[k] = 1
            end
        end
    end

    _grow!(buf.al0, len1 + len2); _grow!(buf.al1, len1 + len2)
    al0 = buf.al0; al1 = buf.al1
    len_al = 0; i = len1; j = len2
    @inbounds while i > 0 || j > 0
        pij = p[idx(i, j)]
        len_al += 1
        if pij == 1
            i -= 1; j -= 1; al0[len_al] = s1[i+1]; al1[len_al] = s2[j+1]
        elseif pij == 2
            j -= 1; al0[len_al] = GAP; al1[len_al] = s2[j+1]
        elseif pij == 3
            i -= 1; al0[len_al] = s1[i+1]; al1[len_al] = GAP
        else
            error("N-W alignment out of range.")
        end
    end
    reverse!(view(al0, 1:len_al)); reverse!(view(al1, 1:len_al))
    buf.len = len_al
    return buf
end

"Mark positions lying in homopolymer runs of length >= 3."
function _homopolymer_mask(s::AbstractVector{UInt8})
    n = length(s)
    h = falses(n)
    i = 1
    for j in 1:n
        if j == n || s[j] != s[j+1]
            j - i >= 2 && (h[i:j] .= true)
            i = j + 1
        end
    end
    return h
end

# Port of nwalign_endsfree_homo: gaps inside homopolymers cost homo_gap_p.
function nwalign_endsfree_homo!(buf::AlignBuffer, s1::AbstractVector{UInt8}, s2::AbstractVector{UInt8},
                                match::Integer, mismatch::Integer, gap_p::Integer, homo_gap_p::Integer,
                                band::Integer)
    len1 = length(s1); len2 = length(s2)
    homo1 = _homopolymer_mask(s1); homo2 = _homopolymer_mask(s2)
    ncol = len2 + 1; nrow = len1 + 1
    _grow!(buf.d32, ncol * nrow); _grow!(buf.p, ncol * nrow)
    d = buf.d32; p = buf.p
    @inline idx(i, j) = i * ncol + j + 1
    @inbounds for i in 0:len1
        d[idx(i, 0)] = 0; p[idx(i, 0)] = 3
    end
    @inbounds for j in 0:len2
        d[idx(0, j)] = 0; p[idx(0, j)] = 2
    end
    if len2 > len1
        lband = band; rband = band + len2 - len1
    elseif len1 > len2
        lband = band + len1 - len2; rband = band
    else
        lband = band; rband = band
    end
    if band >= 0 && (band < len1 || band < len2)
        @inbounds for i in 0:len1
            i - lband - 1 >= 0 && (d[idx(i, i - lband - 1)] = -9999)
            i + rband + 1 <= len2 && (d[idx(i, i + rband + 1)] = -9999)
        end
    end
    @inbounds for i in 1:len1
        l, r = band >= 0 ? (max(i - lband, 1), min(i + rband, len2)) : (1, len2)
        for j in l:r
            left = i == len1 ? d[idx(i, j - 1)] :
                   homo2[j] ? d[idx(i, j - 1)] + homo_gap_p : d[idx(i, j - 1)] + gap_p
            up = j == len2 ? d[idx(i - 1, j)] :
                 homo1[i] ? d[idx(i - 1, j)] + homo_gap_p : d[idx(i - 1, j)] + gap_p
            dg = d[idx(i - 1, j - 1)] + (s1[i] == s2[j] ? match : mismatch)
            k = idx(i, j)
            if up >= dg && up >= left
                d[k] = up; p[k] = 3
            elseif left >= dg
                d[k] = left; p[k] = 2
            else
                d[k] = dg; p[k] = 1
            end
        end
    end
    _grow!(buf.al0, len1 + len2); _grow!(buf.al1, len1 + len2)
    al0 = buf.al0; al1 = buf.al1
    len_al = 0; i = len1; j = len2
    @inbounds while i > 0 || j > 0
        pij = p[idx(i, j)]
        len_al += 1
        if pij == 1
            i -= 1; j -= 1; al0[len_al] = s1[i+1]; al1[len_al] = s2[j+1]
        elseif pij == 2
            j -= 1; al0[len_al] = GAP; al1[len_al] = s2[j+1]
        else
            i -= 1; al0[len_al] = s1[i+1]; al1[len_al] = GAP
        end
    end
    reverse!(view(al0, 1:len_al)); reverse!(view(al1, 1:len_al))
    buf.len = len_al
    return buf
end

# --------------------------------------------------------------------------
# Substitutions of al1 relative to al0 (port of al2subs).
# --------------------------------------------------------------------------

"""
    Sub

Compressed alignment: for each position of the reference (`len0` long),
`map[pos0]` is the aligned position in the query (or `GAP_GLYPH`), and the
`nsubs` substitutions are at reference positions `pos[1:nsubs]`.
Positions are 1-based.
"""
mutable struct Sub
    nsubs::Int
    len0::Int
    map::Vector{UInt16}
    pos::Vector{UInt16}
    nt0::Vector{UInt8}
    nt1::Vector{UInt8}
end
Sub() = Sub(0, 0, UInt16[], UInt16[], UInt8[], UInt8[])

const GAP_GLYPH = 0xffff

@inline _isnt(x::UInt8) = (x >= 0x01) & (x <= 0x05)

"Fill `sub` from the alignment currently held in `buf`."
function al2subs!(sub::Sub, buf::AlignBuffer)
    al0 = buf.al0; al1 = buf.al1; L = buf.len
    len0 = 0; nsubs = 0
    @inbounds for i in 1:L
        a0 = al0[i]; a1 = al1[i]
        n0 = _isnt(a0); n1 = _isnt(a1)
        len0 += n0
        nsubs += (n0 & n1 & (a0 != a1) & (a0 != 0x05) & (a1 != 0x05))
    end
    _grow!(sub.map, len0); _grow!(sub.pos, nsubs); _grow!(sub.nt0, nsubs); _grow!(sub.nt1, nsubs)
    sub.len0 = len0
    s = 0; i0 = 0; i1 = 0
    @inbounds for i in 1:L
        a0 = al0[i]; a1 = al1[i]
        n0 = _isnt(a0); n1 = _isnt(a1)
        i0 += n0; i1 += n1
        if n0
            sub.map[i0] = n1 ? UInt16(i1) : GAP_GLYPH
            if n1 && a0 != a1 && a0 != 0x05 && a1 != 0x05
                s += 1
                sub.pos[s] = UInt16(i0); sub.nt0[s] = a0; sub.nt1[s] = a1
            end
        end
    end
    sub.nsubs = s
    return sub
end

Base.copy(s::Sub) = Sub(s.nsubs, s.len0, s.map[1:s.len0], s.pos[1:s.nsubs], s.nt0[1:s.nsubs], s.nt1[1:s.nsubs])

# --------------------------------------------------------------------------
# Alignment evaluation helpers used by mergePairs / chimera code
# (ports of C_eval_pair and C_pair_consensus). Operate on aligned byte
# strings where gaps are '-'.
# --------------------------------------------------------------------------

"Count (match, mismatch, indel) in the internal (non end-gapped) part of an alignment."
function eval_pair(a0::AbstractVector{UInt8}, a1::AbstractVector{UInt8})
    n = length(a0)
    n == length(a1) || throw(ArgumentError("aligned strings differ in length"))
    g0 = true; g1 = true; start = 0
    while true
        start += 1
        start > n && break
        g0 = g0 && a0[start] == GAP
        g1 = g1 && a1[start] == GAP
        (g0 || g1) || break
    end
    g0 = true; g1 = true; stop = n + 1
    while true
        stop -= 1
        stop < start && break
        g0 = g0 && a0[stop] == GAP
        g1 = g1 && a1[stop] == GAP
        (g0 || g1) || break
    end
    match = 0; mismatch = 0; indel = 0
    @inbounds for i in start:stop
        if a0[i] == GAP || a1[i] == GAP
            indel += 1
        elseif a0[i] == a1[i]
            match += 1
        else
            mismatch += 1
        end
    end
    return (match=match, mismatch=mismatch, indel=indel)
end

"Consensus of two aligned sequences; `prefer` (1 or 2) wins mismatches."
function pair_consensus(a0::AbstractVector{UInt8}, a1::AbstractVector{UInt8}, prefer::Integer, trim_overhang::Bool)
    n = length(a0)
    o = Vector{UInt8}(undef, n)
    @inbounds for i in 1:n
        x = a0[i]; y = a1[i]
        o[i] = x == y ? x : y == GAP ? x : x == GAP ? y : prefer == 1 ? x : prefer == 2 ? y : UInt8('N')
    end
    if trim_overhang
        for i in 1:n
            a0[i] != GAP && break
            o[i] = GAP
        end
        for i in n:-1:1
            a1[i] != GAP && break
            o[i] = GAP
        end
    end
    return filter!(!=(GAP), o)
end
