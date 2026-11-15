# The previous single-alignment implementation of the vectorized aligner
# (commit 0622bb6), kept only as a speed baseline for align_micro.jl.
mutable struct OldBuf
    d::Vector{Int16}
    p::Vector{Int16}
    diag::Vector{Int16}
    al0::Vector{UInt8}
    al1::Vector{UInt8}
    len::Int
end
OldBuf() = OldBuf(Int16[], Int16[], Int16[], UInt8[], UInt8[], 0)

function old_nwalign_vectorized!(buf::OldBuf, s1in::AbstractVector{UInt8}, s2in::AbstractVector{UInt8},
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
