# Compare R and Julia pipeline dumps step by step.
#
#   julia --project=<pkg> compare.jl <r_out> <julia_out> [report.md]
#
# For every step reports whether outputs are identical, and if not, how they
# differ. "iso" rows use R's output from the previous step as input.

using DataFrames, Printf, CodecZlib

rdir, jdir = ARGS[1], ARGS[2]
report = length(ARGS) >= 3 ? ARGS[3] : joinpath(jdir, "parity_report.md")
samples = readlines(joinpath(rdir, "samples.txt"))
lines = String[]
out(s="") = (println(s); push!(lines, s))

rtsv(path) = begin
    L = readlines(path); h = split(L[1], '\t')
    rows = [split(l, '\t'; keepempty=true) for l in L[2:end]]
    DataFrame([Symbol(h[j]) => String[r[j] for r in rows] for j in eachindex(h)])
end
num(s) = s == "NA" || s == "" ? NaN : parse(Float64, s)
rmat(path) = reduce(vcat, [permutedims(num.(split(l, '\t'))) for l in readlines(path)])
maxrel(a, b) = begin
    m = 0.0
    for (x, y) in zip(a, b)
        (isnan(x) && isnan(y)) && continue
        (isnan(x) || isnan(y)) && return Inf
        m = max(m, abs(x - y) / max(abs(x), abs(y), 1e-300))
    end
    m
end
yn(b) = b ? "yes" : "**NO**"

function readfq(path)
    data = read(path); data = transcode(GzipDecompressor, data)
    L = split(String(data), '\n'; keepempty=false)
    [(L[i], L[i+1], L[i+3]) for i in 1:4:length(L)]
end

out("# DenoiseDNA vs R dada2 parity report")
out()
isfile(joinpath(rdir, "versions.txt")) && out("R: " * join(readlines(joinpath(rdir, "versions.txt")), ", "))
out("Samples: $(length(samples))")
out()

for mode in ("e2e", "iso")
    J = joinpath(jdir, mode)
    out("## $(mode == "e2e" ? "End-to-end (Julia from raw reads)" : "Isolated (each step on R's previous output)")")
    out()
    if mode == "e2e"
        rf = rtsv(joinpath(rdir, "filter.tsv")); jf = rtsv(joinpath(J, "filter.tsv"))
        nsame = 0; nreadsame = 0
        for s in samples
            fr = readfq(joinpath(rdir, "filtered", s * "_F_filt.fastq.gz"))
            fj = readfq(joinpath(J, "filtered", s * "_F_filt.fastq.gz"))
            rr = readfq(joinpath(rdir, "filtered", s * "_R_filt.fastq.gz"))
            rj = readfq(joinpath(J, "filtered", s * "_R_filt.fastq.gz"))
            nreadsame += (fr == fj && rr == rj)
        end
        out("- **filterAndTrim**: counts identical: $(yn(rf.reads_out == jf.reads_out)); " *
            "filtered FASTQ records (id, sequence, quality) identical in $nreadsame/$(length(samples)) samples")
    end
    # derep
    nd = 0; nq = 0; nm = 0
    for s in samples, dir in ("F", "R")
        tag = "$(dir)_$s"
        a = rtsv(joinpath(rdir, "derep_$tag.tsv")); b = rtsv(joinpath(J, "derep_$tag.tsv"))
        nd += (a.sequence == b.sequence && a.abundance == b.abundance)
        nq += maxrel(rmat(joinpath(rdir, "derepq_$tag.tsv")), rmat(joinpath(J, "derepq_$tag.tsv"))) <= 1e-15
        nm += readlines(joinpath(rdir, "derepmap_$tag.txt")) == readlines(joinpath(J, "derepmap_$tag.txt"))
    end
    N = 2length(samples)
    out("- **derepFastq**: uniques+order identical $nd/$N; mean qualities identical $nq/$N; read map identical $nm/$N")
    # errors
    for dir in ("F", "R")
        er = rmat(joinpath(rdir, "err_$dir.tsv")); ej = rmat(joinpath(J, "err_$dir.tsv"))
        tr = rmat(joinpath(rdir, "errtrans_$dir.tsv")); tj = rmat(joinpath(J, "errtrans_$dir.tsv"))
        ncr = readline(joinpath(rdir, "errnconsist_$dir.txt")); ncj = readline(joinpath(J, "errnconsist_$dir.txt"))
        out(@sprintf("- **learnErrors (%s)**: rounds R=%s J=%s; transition counts identical: %s; error rates max rel diff %.2e%s",
                     dir, ncr, ncj, yn(tr == tj), maxrel(er, ej), er == ej ? " (bit-identical)" : ""))
    end
    # dada
    for (label, prefix) in (("dada", ""), ("dada pool=TRUE", "pool"), ("dada pool=\"pseudo\"", "pseudo"))
        dirs = prefix == "" ? ("F", "R") : ("F$prefix",)
        nsame = 0; nab = 0; nmap = 0; nn = 0; asvR = 0; asvJ = 0; shared = 0; pmax = 0.0
        for s in samples, dir in dirs
            tag = "$(dir)_$s"
            a = rtsv(joinpath(rdir, "dada_$tag.tsv")); b = rtsv(joinpath(J, "dada_$tag.tsv"))
            nn += 1
            asvR += nrow(a); asvJ += nrow(b); shared += length(intersect(a.sequence, b.sequence))
            nsame += a.sequence == b.sequence
            nab += (a.sequence == b.sequence && a.abundance == b.abundance)
            nmap += readlines(joinpath(rdir, "dadamap_$tag.txt")) == readlines(joinpath(J, "dadamap_$tag.txt"))
            if a.sequence == b.sequence
                pmax = max(pmax, maxrel(num.(a.pval), num.(b.pval)))
            end
        end
        out(@sprintf("- **%s**: ASVs R=%d J=%d shared=%d; identical ASV lists %d/%d; identical abundances %d/%d; identical unique->ASV map %d/%d; cluster p-value max rel diff %.2e",
                     label, asvR, asvJ, shared, nsame, nn, nab, nn, nmap, nn, pmax))
    end
    # merge
    nsame = 0; tot = 0
    for s in samples
        a = rtsv(joinpath(rdir, "merge_$s.tsv")); b = rtsv(joinpath(J, "merge_$s.tsv"))
        cols = [:sequence, :abundance, :forward, :reverse, :nmatch, :nmismatch, :nindel, :accept]
        a.accept = uppercase.(a.accept); b.accept = uppercase.(b.accept)
        nsame += (nrow(a) == nrow(b) && a[:, cols] == b[:, cols]); tot += 1
    end
    out("- **mergePairs** (incl. rejects): identical tables $nsame/$tot")
    # seqtab / chimeras
    sr = rtsv(joinpath(rdir, "seqtab.tsv")); sj = rtsv(joinpath(J, "seqtab.tsv"))
    keyset(df) = Set(zip(df.sample, df.sequence, df.abundance))
    cr = rtsv(joinpath(rdir, "seqtab_cols.tsv")).sequence; cj = rtsv(joinpath(J, "seqtab_cols.tsv")).sequence
    out("- **makeSequenceTable**: entries identical: $(yn(keyset(sr) == keyset(sj))); column order identical: $(yn(cr == cj)) ($(length(cr)) vs $(length(cj)) sequences)")
    for m in ("consensus", "pooled")
        a = rtsv(joinpath(rdir, "nochim_$m.tsv")).sequence; b = rtsv(joinpath(J, "nochim_$m.tsv")).sequence
        out("- **removeBimeraDenovo ($m)**: kept R=$(length(a)) J=$(length(b)); identical: $(yn(a == b))")
    end
    a = rtsv(joinpath(rdir, "nochim_persample.tsv")); b = rtsv(joinpath(J, "nochim_persample.tsv"))
    out("- **removeBimeraDenovo (per-sample)**: entries R=$(nrow(a)) J=$(nrow(b)); identical: $(yn(keyset(a) == keyset(b)))")
    out()
end

# timing
tr = rtsv(joinpath(rdir, "timing.tsv")); tj = rtsv(joinpath(jdir, "e2e", "timing.tsv"))
out("## Timing (seconds, end-to-end)")
out()
out("| step | R | Julia | R/Julia |")
out("|---|---:|---:|---:|")
tjd = Dict(zip(tj.step, num.(tj.seconds)))
for (s, v) in zip(tr.step, num.(tr.seconds))
    j = get(tjd, s, NaN)
    out(@sprintf("| %s | %.2f | %.2f | %.2fx |", s, v, j, v / j))
end
out(@sprintf("| **total** | %.2f | %.2f | %.2fx |", sum(num.(tr.seconds)), sum(num.(tj.seconds)), sum(num.(tr.seconds)) / sum(num.(tj.seconds))))
write(report, join(lines, '\n') * '\n')
println("\nReport written to $report")
