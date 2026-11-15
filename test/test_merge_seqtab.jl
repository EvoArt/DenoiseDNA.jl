@testset "mergePairs and sequence tables" begin
    rng = MersenneTwister(7)
    amp = randseq(rng, 250)
    amp2 = mutate(rng, amp, 5)
    # 150 bp forward reads, 150 bp reverse reads (reverse complemented): 50 bp overlap
    fw(s) = s[1:150]; rv(s) = rc(s)[1:150]
    dir = mktempdir()
    q = fill(36, 150)
    fr = vcat([("p$k", fw(amp), q) for k in 1:200], [("q$k", fw(amp2), q) for k in 1:80])
    rr = vcat([("p$k", rv(amp), q) for k in 1:200], [("q$k", rv(amp2), q) for k in 1:80])
    dF = derepFastq(write_fastq(joinpath(dir, "F.fastq"), fr))
    dR = derepFastq(write_fastq(joinpath(dir, "R.fastq"), rr))
    err = toy_err()
    ddF = dada(dF, err); ddR = dada(dR, err)
    m = mergePairs(ddF, dF, ddR, dR)
    @test nrow(m) == 2
    @test m.sequence == [amp, amp2]
    @test m.abundance == [200, 80]
    @test all(m.nmatch .== 50) && all(m.accept)

    # minOverlap too large -> rejected
    @test nrow(mergePairs(ddF, dF, ddR, dR; minOverlap=60)) == 0
    rej = mergePairs(ddF, dF, ddR, dR; minOverlap=60, returnRejects=true)
    @test nrow(rej) == 2 && !any(rej.accept) && all(rej.sequence .== "")
    jc = mergePairs(ddF, dF, ddR, dR; justConcatenate=true)
    @test jc.sequence[1] == fw(amp) * "N"^10 * rc(rv(amp))

    st = makeSequenceTable([m, m[1:1, :]]; names=["a", "b"])
    @test st.sequences == [amp, amp2]
    @test st.table == [200 80; 200 0]
    @test st.samples == ["a", "b"]
    st2 = makeSequenceTable(["x" => m])
    @test st2.samples == ["x"]
    long = DataFrame(st)
    @test nrow(long) == 3
    # duplicated sequences are summed (and ordered alphabetically, as in R)
    u = getUniques(DataFrame(sequence=["T", "A", "T"], abundance=[1, 2, 3]))
    @test u.sequence == ["A", "T"] && u.abundance == [2, 4]

    # collapseNoMismatch merges a prefix-shorter variant into the full one
    stc = SequenceTable(["a"], [amp, amp[1:240]], [10 3])
    c = collapseNoMismatch(stc)
    @test c.sequences == [amp] && c.table == reshape([13], 1, 1)
end
