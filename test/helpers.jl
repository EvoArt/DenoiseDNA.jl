# Shared helpers for synthetic test data.

const BASES = ['A', 'C', 'G', 'T']

randseq(rng, n) = String(rand(rng, BASES, n))

"Introduce `k` random substitutions."
function mutate(rng, s::String, k::Int)
    c = collect(s)
    for p in randperm(rng, length(c))[1:k]
        c[p] = rand(rng, filter(!=(c[p]), BASES))
    end
    String(c)
end

"Write a FASTQ file; `reads` is a vector of (id, seq, quals::Vector{Int})."
function write_fastq(path, reads)
    open(path, "w") do io
        for (id, s, q) in reads
            println(io, "@", id); println(io, s); println(io, "+")
            println(io, String(Char.(q .+ 33)))
        end
    end
    path
end

"""
Simulate amplicon reads: each true sequence `seqs[i]` contributes `n[i]` reads
with substitution errors drawn from the Phred quality of each position.
"""
function simulate_reads(rng, seqs, n; len=length(seqs[1]), qmean=35)
    reads = Tuple{String,String,Vector{Int}}[]
    k = 0
    for (s, m) in zip(seqs, n), _ in 1:m
        q = clamp.(round.(Int, qmean .- 8 .* (0:len-1) ./ len .+ 2 .* randn(rng, len)), 2, 40)
        c = collect(s)
        for p in 1:len
            if rand(rng) < 10.0^(-q[p] / 10)
                c[p] = rand(rng, filter(!=(c[p]), BASES))
            end
        end
        k += 1
        push!(reads, ("read$k", String(c), q))
    end
    reads[randperm(rng, length(reads))]
end

"A plausible error matrix (rates decreasing with quality)."
function toy_err(nq=41)
    err = zeros(16, nq)
    for q in 0:nq-1
        e = clamp(10.0^(-q / 10), 1e-5, 0.25) / 3
        for i in 0:3, j in 0:3
            err[4i+j+1, q+1] = i == j ? 1 - 3e : e
        end
    end
    err
end
