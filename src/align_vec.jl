# Anti-diagonal ("vectorized") banded Needleman-Wunsch: port of dada2's
# nwalign_vectorized2, the aligner R uses by default.
#
# Each DP cell holds a SIMD vector of W scores, so W query sequences of equal
# length are aligned against the same reference at once (inter-sequence SIMD,
# as in SWIPE/HMMER). The geometry of the band depends only on the two
# lengths, so it is shared by all lanes, and every lane follows exactly the
# recurrence and tie-breaking of a single alignment (W = 1).



"""
    Lanes{W,T}

Cells of `W` lanes stored lane-fastest in a plain `Vector{T}`; indexing loads
or stores a whole cell as a `Vec{W,T}` (unaligned, so no alignment
requirement on the memory). `a[c, s]` indexes column `s` of `stride` cells.
"""
mutable struct Lanes{W,T}
    data::Vector{T}
    stride::Int
    guard::Int        # padding cells before cell 1 (and after the last)
end
Lanes{W,T}() where {W,T} = Lanes{W,T}(T[], 0, 0)
@inline _pos(a::Lanes{W}, k) where {W} = (k - 1 + a.guard) * W + 1
Base.@propagate_inbounds Base.getindex(a::Lanes{W,T}, k::Integer) where {W,T} = vload(Vec{W,T}, a.data, _pos(a, k))
Base.@propagate_inbounds Base.setindex!(a::Lanes{W,T}, v::Vec{W,T}, k::Integer) where {W,T} =
    vstore(v, a.data, _pos(a, k))
Base.@propagate_inbounds Base.getindex(a::Lanes, c::Integer, s::Integer) = a[c+(s-1)*a.stride]
Base.@propagate_inbounds Base.setindex!(a::Lanes, v, c::Integer, s::Integer) = (a[c+(s-1)*a.stride] = v)
Base.eltype(::Lanes{W,T}) where {W,T} = Vec{W,T}
"Element `l` (lane) of cell `k`."
Base.@propagate_inbounds _lane(a::Lanes, k::Integer, l::Integer) = a.data[_pos(a, k)+l-1]
Base.@propagate_inbounds _setlane!(a::Lanes, x, k::Integer, l::Integer) = (a.data[_pos(a, k)+l-1] = x)

# Make room for `n` cells (in columns of `stride` cells) with `guard` spare
# cells on each side. The anti-diagonal recurrence, like the C code it ports,
# reads a few positions just outside the sequences for cells outside the band
# (their values are never used); the guards keep those reads in memory we own.
function _fit!(a::Lanes{W,T}, n, stride=0, guard=2) where {W,T}
    m = (n + 2guard) * W
    length(a.data) < m && resize!(a.data, max(m, 2length(a.data)))
    a.stride = stride; a.guard = guard
    return a
end

"""
    VecBuffer{W}()

Working memory for `W`-lane anti-diagonal alignment.
"""
mutable struct VecBuffer{W}
    d::Lanes{W,Int16}     # ncol x 3: scores of the last three anti-diagonals (row mod 3)
    p::Lanes{W,UInt8}     # nrow*ncol traceback pointers (1 diag, 2 left, 3 up)
    dg::Lanes{W,Int16}    # diagonal scores of the current anti-diagonal
    s1::Lanes{W,UInt8}    # len1 (the shorter sequence), one byte per lane
    s2::Lanes{W,UInt8}    # len2
    len1::Int
    len2::Int
    swap::Bool            # true when the query (lanes) is the longer sequence
    start_col::Int
    ncol::Int
end
VecBuffer{W}() where {W} = VecBuffer{W}(Lanes{W,Int16}(), Lanes{W,UInt8}(), Lanes{W,Int16}(),
                                       Lanes{W,UInt8}(), Lanes{W,UInt8}(), 0, 0, false, 0, 0)
lanes(::VecBuffer{W}) where {W} = W

# C size_t arithmetic: a negative intermediate would wrap to a huge value.
@inline _usz(x::Int) = x < 0 ? typemax(Int) : x

@inline _cell(r, c, ncol) = r * ncol + c + 1      # 0-based (row, col) -> pointer cell
@inline _slot(r) = r % 3 + 1                      # row -> score buffer

"""
    load!(vb, ref, queries)

Set up a batch: align each of `queries` (all the same length; at most W, the
remaining lanes repeat the last query) against `ref`.
"""
function load!(vb::VecBuffer{W}, ref::AbstractVector{UInt8}, queries) where {W}
    nq = length(queries)
    1 <= nq <= W || throw(ArgumentError("need 1 to $W queries"))
    lq = length(first(queries))
    all(q -> length(q) == lq, queries) || throw(ArgumentError("queries must have equal length"))
    lr = length(ref)
    swap = lr > lq                       # s1 must be the shorter sequence
    len1, len2 = swap ? (lq, lr) : (lr, lq)
    g = len2 + 2                         # generous: the band never reaches further
    _fit!(vb.s1, len1, 0, g); _fit!(vb.s2, len2, 0, g)
    a, b = swap ? (vb.s1, vb.s2) : (vb.s2, vb.s1)   # a: query lanes, b: reference
    @inbounds for l in 1:W
        q = queries[min(l, nq)]
        for i in 1:lq
            _setlane!(a, q[i], i, l)
        end
    end
    @inbounds for i in 1:lr, l in 1:W
        _setlane!(b, ref[i], i, l)
    end
    vb.len1 = len1; vb.len2 = len2; vb.swap = swap
    return vb
end

# Set the cells of row `r` that are fixed before the recurrence runs: the band
# boundaries and, near the start, the free left column / top row.
function _init_row!(d, p, r, ncol, start_col, lim_left, lim_top, fill_val::Int16, end_gap_p::Int16)
    s = _slot(r)
    V = eltype(d)
    @inbounds for c in (0, 1, ncol - 2, ncol - 1)
        d[c+1, s] = V(fill_val)
    end
    P = eltype(p)
    if r == 0
        @inbounds d[start_col+1, s] = V(0); p[_cell(0, start_col, ncol)] = P(0)
    end
    if 1 <= r < 1 + lim_left
        c = start_col - 1 - (r - 1) ÷ 2
        @inbounds d[c+1, s] = V(Int16(r * end_gap_p)); p[_cell(r, c, ncol)] = P(3)
    end
    if 1 <= r < 1 + lim_top
        c = start_col + r ÷ 2
        @inbounds d[c+1, s] = V(Int16(r * end_gap_p)); p[_cell(r, c, ncol)] = P(2)
    end
    return nothing
end

"""
    fill_dp!(vb, match, mismatch, gap_p, end_gap_p, band)

Fill the DP for the loaded batch. Scores are kept only for the three
anti-diagonals the recurrence reads; traceback pointers for all cells.
"""
function fill_dp!(vb::VecBuffer{W}, match, mismatch, gap_p, end_gap_p, band) where {W}
    match = Int16(match); mismatch = Int16(mismatch); gap_p = Int16(gap_p); end_gap_p = Int16(end_gap_p)
    len1 = vb.len1; len2 = vb.len2; swap = vb.swap
    band = Int(band)
    band < 0 && (band = len2)
    start_col = 1 + (1 + min(band, len1)) ÷ 2
    ncol = 2 + start_col + min(len2 - len1 + band, len2) ÷ 2
    nrow = len1 + len2 + 1
    vb.start_col = start_col; vb.ncol = ncol
    _fit!(vb.d, 3ncol, ncol); _fit!(vb.p, ncol * nrow); _fit!(vb.dg, ncol)
    d = vb.d; p = vb.p; dg = vb.dg; s1 = vb.s1; s2 = vb.s2
    fill_val = Int16(Int(typemin(Int16)) - min(min(Int(mismatch), Int(gap_p)), min(Int(match), 0)))
    lim_left = min(band, len1); lim_top = min(band + len2 - len1, len2)
    vmatch = Vec{W,Int16}(match); vmismatch = Vec{W,Int16}(mismatch)
    vgap = Vec{W,Int16}(gap_p); vend = Vec{W,Int16}(end_gap_p)
    P1 = Vec{W,UInt8}(0x01); P2 = Vec{W,UInt8}(0x02); P3 = Vec{W,UInt8}(0x03)
    _init_row!(d, p, 0, ncol, start_col, lim_left, lim_top, fill_val, end_gap_p)
    _init_row!(d, p, 1, ncol, start_col, lim_left, lim_top, fill_val, end_gap_p)

    row = 2
    col_min = start_col; col_max = start_col
    i_max = 0; j_min = 0
    even = 1
    recalc_left = false; recalc_right = false
    @inbounds while row <= len1 + len2
        _init_row!(d, p, row, ncol, start_col, lim_left, lim_top, fill_val, end_gap_p)
        s0 = _slot(row); sm1 = _slot(row - 1); sm2 = _slot(row - 2)
        # diagonal scores of this anti-diagonal
        i = i_max; j = j_min
        for c in col_min:col_max
            dg[c+1] = d[c+1, sm2] + vifelse(s1[i+1] == s2[j+1], vmatch, vmismatch)
            i -= 1; j += 1
        end
        left0 = col_min - even + 1                 # buffer column of the first "left" cell
        out0 = _cell(row, col_min, ncol)
        for t in 0:col_max-col_min
            left = d[left0+t, sm1] + vgap
            diag = dg[col_min+t+1]
            up = d[left0+t+1, sm1] + vgap
            # precedence up > left > diag (left > up > diag when swapped)
            takeup = swap ? ~(left >= up) : up >= left
            entry = vifelse(takeup, up, left)
            isdiag = ~(entry >= diag)
            d[col_min+t+1, s0] = vifelse(isdiag, diag, entry)
            p[out0+t] = vifelse(isdiag, P1, vifelse(takeup, P3, P2))
        end

        if row == min(band, len1)
            col_min -= 1; i_max += 1; j_min -= 1
        end
        if row == min(band + len2 - len1, len2)
            col_max += 1
        end

        # ends-free recalculation of the lower boundary: a strictly better free
        # gap replaces the cell, an equal one only takes over the pointer
        if end_gap_p > gap_p
            if recalc_left
                k = _cell(row, col_min, ncol)
                dfree = d[left0, sm1] + vend
                cur = d[col_min+1, s0]; pc = p[k]
                better = dfree > cur
                tie = (dfree == cur) & (swap ? pc != P2 : pc == P1)
                d[col_min+1, s0] = vifelse(better, dfree, cur)
                p[k] = vifelse(better | tie, P2, pc)
            end
            i_max == len1 - 1 && (recalc_left = true)
            if recalc_right
                k = _cell(row, col_max, ncol)
                dfree = d[left0+1+(col_max-col_min), sm1] + vend
                cur = d[col_max+1, s0]; pc = p[k]
                better = dfree > cur
                tie = (dfree == cur) & (swap ? pc == P1 : pc != P3)
                d[col_max+1, s0] = vifelse(better, dfree, cur)
                p[k] = vifelse(better | tie, P3, pc)
            end
            _usz((row + 1) ÷ 2 + col_max - start_col) == len2 && (recalc_right = true)
        end

        # advance the band
        if row < band && row < len1
            even == 1 && (col_min -= 1)
            i_max += 1
        elseif i_max < len1 - 1
            if band % 2 == 0
                even == 1 ? (j_min += 1) : (i_max += 1)
            else
                even == 1 ? (col_min -= 1; i_max += 1) : (col_min += 1; j_min += 1)
            end
        else
            even == 0 && (col_min += 1)
            j_min += 1
        end
        if row < min(band + len2 - len1, len2)
            even == 0 && (col_max += 1)
        elseif _usz((row + 1) ÷ 2 + col_max - start_col) < len2
            if (band + len2 - len1) % 2 == 0
                even == 1 ? (col_max -= 1) : (col_max += 1)
            end
        else
            even == 1 && (col_max -= 1)
        end
        row += 1
        even = 1 - even
    end
    return vb
end

"""
    traceback!(buf::AlignBuffer, vb, lane)

Write the alignment of lane `lane` (reference first, query second) into `buf`.
"""
function traceback!(buf, vb::VecBuffer, lane)
    len1 = vb.len1; len2 = vb.len2; ncol = vb.ncol; start_col = vb.start_col
    p = vb.p; s1 = vb.s1; s2 = vb.s2
    _grow!(buf.al0, len1 + len2); _grow!(buf.al1, len1 + len2)
    a0, a1 = vb.swap ? (buf.al1, buf.al0) : (buf.al0, buf.al1)   # a0 receives s1's column
    n = 0; i = len1; j = len2
    @inbounds while i > 0 || j > 0
        pij = _lane(p, (i + j) * ncol + (2start_col + j - i) ÷ 2 + 1, lane)
        n += 1
        if pij == 1
            i -= 1; j -= 1
            a0[n] = _lane(s1, i + 1, lane); a1[n] = _lane(s2, j + 1, lane)
        elseif pij == 2
            j -= 1
            a0[n] = GAP; a1[n] = _lane(s2, j + 1, lane)
        elseif pij == 3
            i -= 1
            a0[n] = _lane(s1, i + 1, lane); a1[n] = GAP
        else
            error("N-W alignment out of range (len1=$len1, len2=$len2, i=$i, j=$j, p=$pij)")
        end
    end
    reverse!(view(buf.al0, 1:n)); reverse!(view(buf.al1, 1:n))
    buf.len = n
    return buf
end

