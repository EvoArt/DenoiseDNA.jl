# Local regression matching R's `stats::loess` for one predictor.
#
# R's default fit (`surface = "interpolate"`) does not evaluate the local
# regression at each prediction point. It builds a k-d tree over x, fits the
# local regression (value and slope) at the tree's vertices and evaluates by
# cubic Hermite interpolation within each cell. This file ports that algorithm
# (netlib dloess: ehg126/ehg124/ehg125/ehg106/ehg127/ehg139/ehg128) for d = 1,
# gaussian family (no robustness iterations).

using LinearAlgebra: qr, svd, dot

"""
    LoessFit

A fitted 1-D loess. Evaluate with `predict(fit, x)`; points outside the range
of the data give `NaN` (as R gives `NA`).
"""
struct LoessFit
    x::Vector{Float64}
    y::Vector{Float64}
    w::Vector{Float64}
    span::Float64
    degree::Int
    surface::Symbol
    # k-d tree (interpolate surface)
    a::Vector{Int}
    xi::Vector{Float64}
    lo::Vector{Int}
    hi::Vector{Int}
    c::Matrix{Int}
    v::Vector{Float64}
    vval::Matrix{Float64}   # 2 x nv : value, slope
end

# Floyd & Rivest selection (ehg106): permute pi[il:ir] so that p[pi[k]] is the
# k-th smallest. Ported literally so neighbour ordering matches R exactly.
function _ehg106!(il::Int, ir::Int, k::Int, p::AbstractVector{Float64}, pi::Vector{Int})
    l = il; r = ir
    @inbounds while l < r
        t = p[pi[k]]
        i = l; j = r
        pi[l], pi[k] = pi[k], pi[l]
        if t < p[pi[r]]
            pi[l], pi[r] = pi[r], pi[l]
        end
        while i < j
            pi[i], pi[j] = pi[j], pi[i]
            i += 1; j -= 1
            while p[pi[i]] < t
                i += 1
            end
            while t < p[pi[j]]
                j -= 1
            end
        end
        if p[pi[l]] == t
            pi[l], pi[j] = pi[j], pi[l]
        else
            j += 1
            pi[r], pi[j] = pi[j], pi[r]
        end
        j <= k && (l = j + 1)
        k <= j && (r = j - 1)
    end
    return pi
end

# Local fit at q (ehg127). Returns (value, slope). `psi` persists across calls.
function _local_fit(q::Float64, x::Vector{Float64}, y::Vector{Float64}, rw::Vector{Float64},
                    nf::Int, span::Float64, degree::Int, psi::Vector{Int}, dist::Vector{Float64})
    n = length(x)
    @inbounds for i in 1:n
        dist[i] = (x[i] - q)^2
    end
    _ehg106!(1, n, nf, dist, psi)
    rho = dist[psi[nf]] * max(1.0, span)
    rho <= 0 && @warn "loess: span too small"
    k = degree + 1
    Bm = zeros(Float64, nf, k)
    eta = Vector{Float64}(undef, nf)
    @inbounds for i in 1:nf
        wi = sqrt(dist[psi[i]] / rho)
        wi = sqrt(rw[psi[i]] * (1 - wi^3)^3)
        xd = x[psi[i]] - q
        Bm[i, 1] = wi
        degree >= 1 && (Bm[i, 2] = wi * xd)
        degree >= 2 && (Bm[i, 3] = wi * xd^2)
        eta[i] = wi * y[psi[i]]
    end
    colnor = ones(Float64, k)
    @inbounds for j in 1:k
        scal = sqrt(sum(abs2, view(Bm, :, j)))
        if scal > 0
            Bm[:, j] ./= scal
            colnor[j] = scal
        end
    end
    F = qr(Bm)
    qty = F.Q' * eta
    R = F.R[1:k, 1:k]
    S = svd(R)
    sigma = S.S
    tol = sigma[1] * (100 * eps(Float64))
    gamma = [sigma[j] > tol ? dot(view(S.U, :, j), view(qty, 1:k)) / sigma[j] : 0.0 for j in 1:k]
    V = S.V ./ colnor                      # compensate for equilibration
    val = dot(view(V, 1, :), gamma)
    slope = k >= 2 ? dot(view(V, 2, :), gamma) : 0.0
    return val, slope
end

"""
    loess(x, y; weights=ones, span=0.75, degree=2, surface=:interpolate, cell=0.2)

Fit a 1-D loess exactly as R's `loess(y ~ x, weights=weights, span=span,
degree=degree)` with the default gaussian family.
"""
function loess(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};
               weights::AbstractVector{<:Real}=ones(length(x)), span::Real=0.75,
               degree::Integer=2, surface::Symbol=:interpolate, cell::Real=0.2)
    xx = Float64.(collect(x)); yy = Float64.(collect(y)); ww = Float64.(collect(weights))
    n = length(xx)
    (length(yy) == n && length(ww) == n) || throw(DimensionMismatch("x, y and weights must have equal length"))
    n > 0 || throw(ArgumentError("invalid 'x'"))
    0 <= degree <= 2 || throw(ArgumentError("'degree' must be 0, 1 or 2"))
    span = Float64(span)
    nf = min(n, floor(Int, n * span + 1e-5))
    nf <= 0 && error("span is too small")
    if surface === :direct
        return LoessFit(xx, yy, ww, span, degree, :direct, Int[], Float64[], Int[], Int[],
                        zeros(Int, 2, 0), Float64[], zeros(2, 0))
    end
    surface === :interpolate || throw(ArgumentError("surface must be :interpolate or :direct"))

    fc = floor(Int, n * (span * Float64(cell)))
    nvmax = max(200, n); ncmax = nvmax
    # bounding box, expanded slightly (ehg126)
    alpha, beta = extrema(xx)
    mu = 0.005 * max(beta - alpha, 1e-10 * max(abs(alpha), abs(beta)) + 1e-30)
    v = [alpha - mu, beta + mu]
    # Cells. For a cell awaiting processing lo/hi hold its point range
    # (into the permutation `pi`); once split they hold its children.
    a = [0]; xi = [0.0]; lo = [1]; hi = [n]
    cmat = Vector{Vector{Int}}([[1, 2]])   # vertex indices of each cell
    pi = collect(1:n)
    p = 1
    while p <= length(a)
        l = lo[p]; u = hi[p]
        diam = abs(v[cmat[p][2]] - v[cmat[p][1]])
        leaf = (u - l + 1 <= fc) || (diam <= 0.0)
        if !leaf
            leaf = ncmax < length(a) + 2 || nvmax < length(v) + 1
        end
        m = 0
        if !leaf
            m = (l + u) ÷ 2
            _ehg106!(l, u, m, xx, pi)
            offset = 0
            while (m + offset) < u && (m + offset) >= l
                if offset < 0
                    lower = l; check = m + offset; upper = check
                else
                    lower = m + offset + 1; check = lower; upper = u
                end
                _ehg106!(lower, upper, check, xx, pi)
                if xx[pi[m+offset]] == xx[pi[m+offset+1]]
                    offset = -offset
                    offset >= 0 && (offset += 1)
                else
                    m = m + offset
                    break
                end
            end
            leaf = v[cmat[p][1]] == xx[pi[m]] || v[cmat[p][2]] == xx[pi[m]]
        end
        if leaf
            a[p] = 0
        else
            a[p] = 1
            t = xx[pi[m]]
            xi[p] = t
            # new vertex (ehg125), reusing an identical existing one
            vm = findfirst(==(t), v)
            if vm === nothing
                push!(v, t); vm = length(v)
            end
            push!(a, 0, 0); push!(xi, 0.0, 0.0)
            push!(lo, l, m + 1); push!(hi, m, u)
            nc = length(a)
            lo_child = nc - 1; hi_child = nc
            push!(cmat, [cmat[p][1], vm], [vm, cmat[p][2]])
            # the parent's lo/hi now point at its children
            lo[p] = lo_child; hi[p] = hi_child
        end
        p += 1
    end
    # NB: for internal cells lo/hi are child indices; for leaves they are point ranges.
    nv = length(v)
    vval = zeros(Float64, 2, nv)
    psi = collect(1:n)
    dist = Vector{Float64}(undef, n)
    for l in 1:nv
        vval[1, l], vval[2, l] = _local_fit(v[l], xx, yy, ww, nf, span, degree, psi, dist)
    end
    cm = reduce(hcat, cmat)
    return LoessFit(xx, yy, ww, span, degree, :interpolate, a, xi, lo, hi, cm, v, vval)
end

function _interp(fit::LoessFit, z::Float64)
    j = 1
    while fit.a[j] != 0
        j = z <= fit.xi[j] ? fit.lo[j] : fit.hi[j]
    end
    ll = fit.c[1, j]; ur = fit.c[2, j]
    v0 = fit.v[ll]; v1 = fit.v[ur]
    h = (z - v0) / (v1 - v0)
    phi0 = (1 - h)^2 * (1 + 2h)
    phi1 = h^2 * (3 - 2h)
    psi0 = h * (1 - h)^2
    psi1 = h^2 * (h - 1)
    return phi0 * fit.vval[1, ll] + phi1 * fit.vval[1, ur] +
           (psi0 * fit.vval[2, ll] + psi1 * fit.vval[2, ur]) * (v1 - v0)
end

"""
    predict(fit::LoessFit, x)

Evaluate the loess fit. Like R, returns `NaN` outside the range of the data
for the interpolated surface.
"""
function predict(fit::LoessFit, x::AbstractVector{<:Real})
    lo, hi = extrema(fit.x)
    out = fill(NaN, length(x))
    if fit.surface === :direct
        n = length(fit.x)
        nf = min(n, floor(Int, n * fit.span + 1e-5))
        psi = collect(1:n); dist = Vector{Float64}(undef, n)
        for (k, z) in enumerate(x)
            out[k] = _local_fit(Float64(z), fit.x, fit.y, fit.w, nf, fit.span, fit.degree, psi, dist)[1]
        end
        return out
    end
    for (k, z) in enumerate(x)
        (lo <= z <= hi) && (out[k] = _interp(fit, Float64(z)))
    end
    return out
end
predict(fit::LoessFit, x::Real) = predict(fit, [x])[1]
