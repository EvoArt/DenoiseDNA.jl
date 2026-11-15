@testset "alignment" begin
    rng = MersenneTwister(42)
    buf = D.AlignBuffer()
    al(a, b; band=16, f=D.nwalign_vectorized!) = begin
        f(buf, Vector{UInt8}(codeunits(a)), Vector{UInt8}(codeunits(b)), 5, -4, -8, 0, band)
        (String(buf.al0[1:buf.len]), String(buf.al1[1:buf.len]))
    end
    s = randseq(rng, 120)
    @test al(s, s) == (s, s)
    # one substitution
    t = s[1:59] * (s[60] == 'A' ? "C" : "A") * s[61:end]
    a0, a1 = al(s, t)
    @test a0 == s && a1 == t
    # deletion in the query
    t = s[1:50] * s[52:end]
    a0, a1 = al(s, t)
    @test replace(a0, "-" => "") == s && replace(a1, "-" => "") == t
    @test count(==('-'), a1) == 1
    # ends-free: an overhang costs nothing
    t = s[11:end]
    a0, a1 = al(s, t)
    @test a1 == "-"^10 * t
    # both orientations of lengths give consistent results
    for _ in 1:50
        x = randseq(rng, rand(rng, 30:80)); y = mutate(rng, x, 3)[1:end-rand(rng, 0:5)]
        a0, a1 = al(x, y); b1, b0 = al(y, x)
        @test replace(a0, "-" => "") == x && replace(a1, "-" => "") == y
        @test replace(b0, "-" => "") == x && replace(b1, "-" => "") == y
    end

    # batched (multi-lane) alignment == single alignment in every lane,
    # with the reference shorter, equal and longer than the queries
    vb = D.VecBuffer{16}()
    one = D.AlignBuffer(); lanebuf = D.AlignBuffer()
    for reflen in (140, 150, 160), band in (16, 5, -1), nq in (1, 7, 16)
        ref = D.encode_nt(randseq(rng, reflen))
        qs = [D.encode_nt(mutate(rng, String(D.NT_CHAR[ref[1:min(end, 150)]]) * randseq(rng, max(0, 150 - reflen)), 4))
              for _ in 1:nq]
        for q in qs   # indels too
            p = rand(rng, 10:140); q[p] = q[p+1]
        end
        D.load!(vb, ref, qs)
        D.fill_dp!(vb, 5, -4, -8, 0, band)
        for (l, q) in enumerate(qs)
            D.traceback!(lanebuf, vb, l)
            D.nwalign_vectorized!(one, ref, q, 5, -4, -8, 0, band)
            @test lanebuf.al0[1:lanebuf.len] == one.al0[1:one.len]
            @test lanebuf.al1[1:lanebuf.len] == one.al1[1:one.len]
        end
    end
    @test_throws ArgumentError D.load!(vb, D.encode_nt("ACGTACGT"), [D.encode_nt("ACGT"), D.encode_nt("ACG")])

    # standard banded endsfree and gapless
    D.nwalign_endsfree!(buf, Vector{UInt8}(codeunits(s)), Vector{UInt8}(codeunits(s[5:end])), 5, -4, -8, -1)
    @test String(buf.al1[1:buf.len]) == "----" * s[5:end]
    D.nwalign_gapless!(buf, Vector{UInt8}(codeunits("ACGT")), Vector{UInt8}(codeunits("AC")))
    @test String(buf.al1[1:buf.len]) == "AC--"

    # substitutions
    sub = D.Sub()
    D.nwalign_vectorized!(buf, D.encode_nt("ACGTACGT"), D.encode_nt("ACCTACGA"), 5, -4, -8, 0, 16)
    D.al2subs!(sub, buf)
    @test sub.nsubs == 2
    @test sub.pos[1:2] == [3, 8]
    @test sub.map[1:8] == 1:8

    # evaluation helpers
    ev = D.eval_pair(Vector{UInt8}(codeunits("--ACGTA-")), Vector{UInt8}(codeunits("TTACCT-A")))
    @test ev == (match=3, mismatch=1, indel=1)
    @test String(D.pair_consensus(Vector{UInt8}(codeunits("ACGT--")), Vector{UInt8}(codeunits("--GTCA")), 1, false)) == "ACGTCA"
    @test String(D.pair_consensus(Vector{UInt8}(codeunits("ACGT--")), Vector{UInt8}(codeunits("--GACA")), 2, false)) == "ACGACA"

    @test nwalign("ACGTACGT", "ACGACGT")[2] == "ACG-ACGT"
    @test nwhamming("ACGTACGTAA", "ACGTACGTAA") == 0
    @test rc("AACGTN") == "NACGTT"
end
