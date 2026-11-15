# Line-level profile and allocation count of the batched DP kernel.
using DenoiseDNA, Profile, ProfileToLLM, Random, InteractiveUtils
using DenoiseDNA: VecBuffer, load!, fill_dp!
rng = MersenneTwister(1)
ref = rand(rng, 0x01:0x04, 240)
for W in (8, 32)
    vb = VecBuffer{W}()
    qs = [(c = copy(ref); c[rand(rng, 1:240)] = 0x02; c) for _ in 1:W]
    load!(vb, ref, qs); fill_dp!(vb, 5, -4, -8, 0, 16)
    a = @allocated fill_dp!(vb, 5, -4, -8, 0, 16)
    t = @elapsed for _ in 1:200; fill_dp!(vb, 5, -4, -8, 0, 16); end
    println("W=$W  allocated per call: $a bytes   $(round(1e6t/200, digits=1)) us per batch")
    if W == 32
        Profile.clear(); Profile.init(n=10^7, delay=0.0001)
        @profile for _ in 1:3000; fill_dp!(vb, 5, -4, -8, 0, 16); end
        print_profile(stdout, profile_table(Profile.fetch(); sortby=:self); max_rows=30, backend=:markdown)
        io = IOBuffer(); code_native(io, fill_dp!, typeof.((vb, 5, -4, -8, 0, 16)); debuginfo=:none)
        asm = String(take!(io))
        println("\ncalls in native code: ", [m.match for m in eachmatch(r"call\s+\S+.*", asm)][1:min(end, 40)])
    end
end
