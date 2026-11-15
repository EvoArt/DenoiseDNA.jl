# Displays of the objects the package returns.
#
# `show(io, x)` is one line (for vectors, tuples, string interpolation);
# `show(io, MIME"text/plain"(), x)` — what the REPL prints — is a summary:
# a headline, a few key-value lines, and the top rows of the data, sized to
# the terminal. Helpers are in style.jl.

# ------------------------------------------------------------------ helpers

# abundance-weighted mean of each column of `quals` (NaN = past the end)
function _mean_qual(quals::AbstractMatrix{<:Real}, w::AbstractVector{<:Real})
    out = fill(NaN, size(quals, 2))
    for p in axes(quals, 2)
        s = 0.0; c = 0.0
        for i in axes(quals, 1)
            q = quals[i, p]
            isnan(q) && continue
            s += w[i] * q; c += w[i]
        end
        c > 0 && (out[p] = s / c)
    end
    return out
end

# "Q37 → Q29": mean quality of the first and last 10 positions
function _qtrend(io, q)
    v = filter(!isnan, q)
    isempty(v) && return ""
    k = min(10, length(v))
    a, b = sum(v[1:k]) / k, sum(v[end-k+1:end]) / k
    return _sty(io, "Q$(round(Int, a))"; color=_qcolor(a)) * _note(io, " → ") * _sty(io, "Q$(round(Int, b))"; color=_qcolor(b))
end

# "240 bp" or "251–254 bp · 253 bp in 92%" (weighted by `w`)
function _lengths(io, seqs::AbstractVector{<:AbstractString}, w::AbstractVector{<:Real})
    isempty(seqs) && return "–"
    lens = length.(seqs)
    lo, hi = extrema(lens)
    lo == hi && return "$lo bp"
    d = Dict{Int,Float64}()
    for (l, x) in zip(lens, w)
        d[l] = get(d, l, 0.0) + x
    end
    mode = argmax(d)
    return "$(lo)–$(hi) bp" * _note(io, " · $mode bp in $(_pct(d[mode] / sum(w)))")
end

# width left for a sequence column after `used` characters of other columns
_seqw(io, used) = max(12, _termwidth(io) - used - 2)

function _changed(o::DadaOptions)
    d = DadaOptions()
    return [(f, getfield(o, f), getfield(d, f)) for f in fieldnames(DadaOptions) if !isequal(getfield(o, f), getfield(d, f))]
end

_optval(v::Float64) = isinteger(v) && abs(v) < 1e6 ? string(Int(v)) : _g(v)
_optval(v) = string(v)

function _options_line(io, o::DadaOptions)
    ch = _changed(o)
    isempty(ch) && return _note(io, "defaults")
    return join((string(f, "=", _sty(io, _optval(v); bold=true)) for (f, v, _) in ch), _note(io, ", "))
end

function _errmodel_line(io, err_in, err_out)
    err_out === nothing && return _note(io, "none estimated")
    nq = size(err_out, 2)
    if err_in isa AbstractVector
        conv = any(==(err_out), err_in)
        return "self-consistent, 16 × $nq " *
               (conv ? _note(io, "(converged after $(length(err_in)) rounds)") :
                _sty(io, "(not converged after $(length(err_in)) rounds)"; color=_PAL.warn))
    end
    return "16 × $nq " * _note(io, "(given)")
end

# ------------------------------------------------------------------ FastqRecords

Base.show(io::IO, f::FastqRecords) = print(io, "FastqRecords(", _fmt(length(f)), " reads)")

function Base.show(io::IO, ::MIME"text/plain", f::FastqRecords)
    n = length(f)
    lens = length.(f.seq)
    _title(io, "FastqRecords", "$(_fmt(n)) reads" * (n == 0 ? "" : _note(io, " · $(_fmt(sum(lens))) bases")))
    n == 0 && return
    lo, hi = extrema(lens)
    _kv(io, "length", lo == hi ? "$lo bp" : "$(lo)–$(hi) bp")
    q = zeros(hi); c = zeros(Int, hi)
    for r in f.qual, (p, x) in enumerate(r)
        q[p] += x; c[p] += 1
    end
    mq = [c[p] > 0 ? q[p] / c[p] : NaN for p in 1:hi]
    W = _termwidth(io)
    _kv(io, "quality", _qspark(io, mq, W - _KEY_W - 16) * "  " * _qtrend(io, mq))
    k = _nshow(io, n, 5)
    _section(io, "First reads", min(W, 80))
    # Illumina ids share a long instrument/run prefix; the tile:x:y that tells
    # reads apart is at the end of the first word, so that is what is kept
    ids = [first(split(f.id[i]; limit=2), 1)[1] for i in 1:k]
    idw = min(24, maximum(length, ids))
    for i in 1:k
        id = length(ids[i]) > idw ? "…" * last(ids[i], idw - 1) : ids[i]
        println(io, "  ", _note(io, rpad(id, idw)), "  ", _seq(io, String(copy(f.seq[i])), _seqw(io, idw + 4)))
    end
    _more(io, n, k, "reads")
end

# ------------------------------------------------------------------ Derep

Base.show(io::IO, d::Derep) = print(io, "Derep(", _fmt(sum(d.abundance)), " reads, ", _fmt(length(d)), " uniques)")

function Base.show(io::IO, ::MIME"text/plain", d::Derep)
    nr = sum(d.abundance); nu = length(d)
    _title(io, "Derep", "$(_fmt(nr)) reads in $(_fmt(nu)) unique sequences")
    nu == 0 && return
    _kv(io, "length", _lengths(io, d.sequence, d.abundance))
    ns = count(==(1), d.abundance)
    _kv(io, "singletons", "$(_fmt(ns)) uniques" * _note(io, " · $(_pct(ns / nr)) of reads"))
    W = _termwidth(io)
    if size(d.quals, 2) > 0
        mq = _mean_qual(d.quals, d.abundance)
        _kv(io, "quality", _qspark(io, mq, W - _KEY_W - 16) * "  " * _qtrend(io, mq))
    end
    k = _nshow(io, nu, 8)
    _section(io, "Most abundant", min(W, 100))
    nw = ndigits(nu)
    rows = [[string(i), _fmt(d.abundance[i]), _pct(d.abundance[i] / nr)] for i in 1:k]
    used = 2 + nw + 2 + length(_fmt(d.abundance[1])) + 2 + 6 + 2
    for r in rows
        push!(r, _seq(io, d.sequence[parse(Int, r[1])], _seqw(io, used)))
    end
    _table(io, ["#", "reads", "%", "sequence"], rows; align=[:r, :r, :r, :l])
    _more(io, nu, k, "uniques")
end

function Base.show(io::IO, ::MIME"text/plain", ds::Vector{Derep})
    n = length(ds)
    nr = sum(d -> sum(d.abundance), ds; init=0)
    _title(io, "$n-element Vector{Derep}", "$(_fmt(nr)) reads in $(_fmt(sum(length, ds; init=0))) uniques")
    n == 0 && return
    W = _termwidth(io)
    k = _nshow(io, n, 20)
    rows = Vector{String}[]
    qw = clamp(W - 55, 10, 40)
    for i in 1:k
        d = ds[i]
        r = sum(d.abundance)
        lens = isempty(d.sequence) ? (0, 0) : extrema(length, d.sequence)
        q = size(d.quals, 2) > 0 ? _qspark(io, _mean_qual(d.quals, d.abundance), qw) : ""
        push!(rows, [string(i), _fmt(r), _fmt(length(d)), _pct(count(==(1), d.abundance) / max(r, 1)),
                     lens[1] == lens[2] ? "$(lens[1])" : "$(lens[1])–$(lens[2])", q])
    end
    _table(io, ["#", "reads", "uniques", "singletons", "length", "quality"], rows;
           align=[:r, :r, :r, :r, :r, :l])
    _more(io, n, k, "samples")
end

# ------------------------------------------------------------------ DadaResult

Base.show(io::IO, r::DadaResult) =
    print(io, "DadaResult(", _fmt(length(r.sequence)), " ASVs from ", _fmt(length(r.map)), " uniques)")

function Base.show(io::IO, ::MIME"text/plain", r::DadaResult)
    na = length(r.sequence); nu = length(r.map)
    _title(io, "DadaResult", "$(_fmt(na)) ASVs from $(_fmt(nu)) unique sequences")
    nr = sum(r.abundance)
    na == 0 && return
    _kv(io, "reads", "$(_fmt(nr)) in ASVs" * _note(io, " · top ASV $(_pct(r.abundance[1] / nr))"))
    nunc = count(==(0), r.map)
    nunc > 0 && _kv(io, "uncorrected", "$(_fmt(nunc)) uniques" * _note(io, " not assigned (p < OMEGA_C)"))
    cl = r.clustering
    if na > 1 && hasproperty(cl, :birth_ham)
        hams = cl.birth_ham[2:end]
        fold = cl.birth_fold[2:end]
        med = sort(hams)[cld(length(hams), 2)]
        n1 = count(==(1), hams)
        _kv(io, "births", "distance to parent: median $med" *
                         (n1 > 0 ? ", $(_fmt(n1)) at 1" : "") *
                         _note(io, " · min fold $(_g(minimum(fold)))"))
    end
    _kv(io, "error model", _errmodel_line(io, r.err_in, r.err_out))
    _kv(io, "options", _options_line(io, r.opts))
    W = _termwidth(io)
    k = _nshow(io, na, 8)
    _section(io, "Sequence variants", min(W, 100))
    has = hasproperty(cl, :birth_ham)
    rows = Vector{String}[]
    for i in 1:k
        row = [string(i), _fmt(r.abundance[i]), _pct(r.abundance[i] / nr)]
        if has
            push!(row, _fmt(cl.nunq[i]),
                  i == 1 || cl.birth_from[i] == 0 ? "–" : string(cl.birth_from[i]),
                  i == 1 || cl.birth_from[i] == 0 ? "–" : string(cl.birth_ham[i]),
                  i == 1 || cl.birth_from[i] == 0 ? "–" : _g(cl.birth_fold[i]))
        end
        push!(rows, row)
    end
    hdr = has ? ["#", "reads", "%", "uniques", "from", "ham", "fold"] : ["#", "reads", "%"]
    used = 2 + sum(j -> maximum(_vis(rw[j]) for rw in vcat([hdr], rows)) + 2, eachindex(hdr))
    for (i, row) in enumerate(rows)
        push!(row, _seq(io, r.sequence[i], _seqw(io, used)))
    end
    _table(io, vcat(hdr, "sequence"), rows; align=vcat(fill(:r, length(hdr)), :l))
    _more(io, na, k, "ASVs")
end

function Base.show(io::IO, ::MIME"text/plain", rs::Vector{DadaResult})
    n = length(rs)
    distinct = length(unique(s for r in rs for s in r.sequence))
    _title(io, "$n-element Vector{DadaResult}",
           "$(_fmt(distinct)) distinct ASVs" * _note(io, " · $(_fmt(sum(r -> sum(r.abundance), rs; init=0))) reads"))
    n == 0 && return
    _kv(io, "error model", _errmodel_line(io, rs[1].err_in, rs[1].err_out))
    _kv(io, "options", _options_line(io, rs[1].opts))
    println(io)
    k = _nshow(io, n, 20)
    rows = Vector{String}[]
    for i in 1:k
        r = rs[i]
        nr = sum(r.abundance)
        push!(rows, [string(i), _fmt(length(r.sequence)), _fmt(length(r.map)), _fmt(nr),
                     isempty(r.abundance) ? "–" : _pct(r.abundance[1] / nr), _fmt(count(==(0), r.map))])
    end
    _table(io, ["#", "ASVs", "uniques", "reads", "top ASV", "uncorrected"], rows)
    _more(io, n, k, "samples")
end

# ------------------------------------------------------------------ ErrorModel

const _NT = ('A', 'C', 'G', 'T')
const _OFFDIAG = [r for r in 1:16 if (r - 1) ÷ 4 != (r - 1) % 4]
_trans_name(r) = string(_NT[(r-1)÷4+1], "→", _NT[(r-1)%4+1])

_converged(e::ErrorModel) = any(==(e.err_out), e.err_in)

Base.show(io::IO, e::ErrorModel) =
    print(io, "ErrorModel(16 × ", size(e.err_out, 2), ", ", length(e.err_in), " rounds", _converged(e) ? "" : ", not converged", ")")

function Base.show(io::IO, ::MIME"text/plain", e::ErrorModel)
    err = e.err_out
    nq = size(err, 2)
    conv = _converged(e)
    nobs = sum(e.trans)
    _title(io, "ErrorModel", "16 × $nq error rates" * _note(io, " · $(_fmt(nobs)) bases observed"))
    # quality scores where something was observed
    obs = [q <= size(e.trans, 2) && sum(view(e.trans, :, q)) > 0 for q in 1:nq]
    any(obs) || (obs .= true)
    lerr = [err[r, q] > 0 ? log10(err[r, q]) : NaN for r in 1:16, q in 1:nq]
    vals = [lerr[r, q] for r in _OFFDIAG for q in 1:nq if obs[q] && !isnan(lerr[r, q])]
    lo, hi = isempty(vals) ? (-5.0, -1.0) : extrema(vals)
    W = _termwidth(io)
    sw = min(nq, W - 40)
    qcols = [q - 1 for q in (31, 41) if q <= nq]                 # Phred 30, 40
    _section(io, "Fitted rate by quality score", min(W, 4 + 4 + sw + 11 * length(qcols) + 14))
    cap = rpad("Q0", sw - length("Q$(nq-1)")) * "Q$(nq-1)"
    hdr = "        " * _note(io, cap) * join(_note(io, lpad("Q$q", 10)) for q in qcols) * _note(io, lpad("subs seen", 11))
    println(io, hdr)
    for r in _OFFDIAG
        row = lerr[r, :]
        sp = _spark(io, _bin(row, sw), lo, hi;
                    color=(k, x) -> obs[clamp(round(Int, (k - 0.5) * nq / sw + 0.5), 1, nq)] ? _PAL.spark : _PAL.tail)
        n = r <= size(e.trans, 1) ? sum(view(e.trans, r, :)) : 0
        println(io, "  ", _sty(io, _trans_name(r); bold=true), "  ", sp,
                join(lpad(@sprintf("%.1e", err[r, q+1]), 10) for q in qcols), lpad(_fmt(n), 11))
    end
    println(io, _note(io, "        log scale $(@sprintf("%.0e", 10^lo))–$(@sprintf("%.0e", 10^hi)), all rows · dim: no data"))
    println(io)
    items = Tuple{Bool,String}[]
    push!(items, conv ? (true, "converged after $(length(e.err_in)) rounds") :
                 (false, "not converged after $(length(e.err_in)) rounds"))
    qs = findall(obs)
    rising = [_trans_name(r) for r in _OFFDIAG if err[r, qs[end]] > err[r, qs[1]]]
    push!(items, isempty(rising) ? (true, "rates fall with quality") :
                 (false, join(rising, ", ") * " rise with quality"))
    _checks(io, items)
    conv || println(io, _sty(io, "   increase MAX_CONSIST, or use more reads (nbases)"; color=_PAL.warn))
end

# ------------------------------------------------------------------ SequenceTable

Base.show(io::IO, st::SequenceTable) =
    print(io, "SequenceTable(", length(st.samples), " samples × ", _fmt(length(st.sequences)), " ASVs)")

# a sample's reads over the top ASVs: shade by share of the sample (log steps)
function _strip(io, row::AbstractVector{<:Integer}, tot)
    b = IOBuffer()
    for x in row
        f = tot > 0 ? x / tot : 0.0
        print(b, x == 0 ? '·' : f >= 0.1 ? '█' : f >= 0.01 ? '▓' : f >= 0.001 ? '▒' : '░')
    end
    s = String(take!(b))
    _hascolor(io) || return s
    return replace(s, "·" => _sty(io, "·"; color=_PAL.tail), r"[░▒▓█]+" => m -> _sty(io, m; color=_PAL.spark))
end

function Base.show(io::IO, ::MIME"text/plain", st::SequenceTable)
    ns, na = size(st.table)
    tot = sum(st.table)
    _title(io, "SequenceTable", "$ns samples × $(_fmt(na)) ASVs" * _note(io, " · $(_fmt(tot)) reads"))
    (ns == 0 || na == 0) && return
    colsum = vec(sum(st.table; dims=1))
    rowsum = vec(sum(st.table; dims=2))
    prev = vec(sum(st.table .> 0; dims=1))
    _kv(io, "ASV length", _lengths(io, st.sequences, ones(na)))
    sr = sort(rowsum)
    _kv(io, "reads", ns == 1 ? _fmt(tot) :
                     "$(_fmt(sr[1]))–$(_fmt(sr[end])) per sample" * _note(io, " · median $(_fmt(sr[cld(ns, 2)]))"))
    if ns > 1
        _kv(io, "prevalence", "$(_fmt(count(==(1), prev))) ASVs in 1 sample" *
                             _note(io, " · $(_fmt(count(==(ns), prev))) in all $ns"))
    end
    W = _termwidth(io)
    # samples, with a strip over the most abundant ASVs (table columns are
    # ordered by abundance unless built with orderBy=nothing)
    ord = sortperm(colsum; rev=true)
    k = _nshow(io, ns, 12)
    nmw = min(20, maximum(length, st.samples[1:k]))
    rows = Vector{String}[]
    for i in 1:k
        nm = length(st.samples[i]) > nmw ? first(st.samples[i], nmw - 1) * "…" : st.samples[i]
        push!(rows, [nm, _fmt(rowsum[i]), _fmt(count(>(0), view(st.table, i, :)))])
    end
    hdr = ["sample", "reads", "ASVs"]
    used = 2 + sum(j -> maximum(_vis(rw[j]) for rw in vcat([hdr], rows)) + 2, eachindex(hdr))
    nstrip = min(na, W - used - 1)
    for (i, row) in enumerate(rows)
        push!(row, _strip(io, view(st.table, i, ord[1:nstrip]), rowsum[i]))
    end
    _section(io, "Samples", min(W, max(used + nstrip, 60)))
    _table(io, vcat(hdr, "ASVs by total abundance →"), rows; align=[:l, :r, :r, :l])
    _more(io, ns, k, "samples")
    println(io, _note(io, "  share of reads: ") *
                join((_sty(io, g; color=_PAL.spark) * _note(io, " " * t) for (g, t) in
                      (("█", "≥10%"), ("▓", "≥1%"), ("▒", "≥0.1%"), ("░", "less"))), "  ") *
                "  " * _sty(io, "·"; color=_PAL.tail) * _note(io, " absent"))
    # the most abundant ASVs
    ka = _nshow(io, na, 8)
    _section(io, "Most abundant ASVs", min(W, 100))
    rows = [[string(ord[j]), _fmt(colsum[ord[j]]), _pct(colsum[ord[j]] / tot), string(prev[ord[j]]),
             string(length(st.sequences[ord[j]]))] for j in 1:ka]
    hdr = ["#", "reads", "%", "samples", "bp"]
    used = 2 + sum(j -> maximum(_vis(rw[j]) for rw in vcat([hdr], rows)) + 2, eachindex(hdr))
    for (j, row) in enumerate(rows)
        push!(row, _seq(io, st.sequences[ord[j]], _seqw(io, used)))
    end
    _table(io, vcat(hdr, "sequence"), rows; align=[:r, :r, :r, :r, :r, :l])
    _more(io, na, ka, "ASVs")
end

# ------------------------------------------------------------------ DadaOptions

function Base.show(io::IO, o::DadaOptions)
    ch = _changed(o)
    print(io, "DadaOptions(", join((string(f, "=", repr(v)) for (f, v, _) in ch), ", "), ")")
end

function Base.show(io::IO, ::MIME"text/plain", o::DadaOptions)
    ch = _changed(o)
    _title(io, "DadaOptions", isempty(ch) ? "defaults" : "$(length(ch)) changed from the defaults")
    fs = fieldnames(DadaOptions)
    chs = Set(first.(ch))
    cells = map(fs) do f
        v = _optval(getfield(o, f))
        f in chs ? (string(f), _sty(io, v; color=_PAL.amber, bold=true), _note(io, " (default " * _optval(getfield(DadaOptions(), f)) * ")")) :
                   (string(f), v, "")
    end
    kw = maximum(c -> length(c[1]), cells)
    cw = kw + 2 + maximum(c -> _vis(c[2]) + _vis(c[3]), cells)
    ncol = max(1, min(3, (_termwidth(io) - 2) ÷ (cw + 3)))
    nrow = cld(length(cells), ncol)
    for i in 1:nrow
        line = "  "
        for j in 1:ncol
            k = (j - 1) * nrow + i
            k > length(cells) && break
            f, v, d = cells[k]
            line *= _rpad(_note(io, rpad(f, kw)) * "  " * v * d, cw + 3)
        end
        println(io, rstrip(line))
    end
end
