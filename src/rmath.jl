# Poisson tail probabilities from R's own nmath library (Rmath_jll), so that
# abundance p-values — including underflow to exactly 0.0, which drives
# tie-breaking in cluster budding — are bit-identical to R/dada2.

"Upper tail P(X > x) for X ~ Poisson(lambda), computed by R's `ppois`."
@inline ppois_upper(x::Real, lambda::Real) =
    ccall((:ppois, libRmath), Cdouble, (Cdouble, Cdouble, Cint, Cint), x, lambda, 0, 0)

const TAIL_APPROX_CUTOFF = 1e-7

"""
    calc_pA(reads, E_reads, prior)

Abundance p-value: probability of seeing at least `reads` copies when
`E_reads` are expected. Unless `prior`, the p-value is conditioned on the
sequence being observed at all (port of dada2 `calc_pA`).
"""
function calc_pA(reads::Integer, E_reads::Float64, prior::Bool)
    pval = ppois_upper(reads - 1, E_reads)
    if !prior
        norm = 1.0 - exp(-E_reads)
        if norm < TAIL_APPROX_CUTOFF
            norm = E_reads - 0.5 * E_reads * E_reads
        end
        pval = pval / norm
    end
    return pval
end
