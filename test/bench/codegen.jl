# Inspect the compiled aligner kernel: type stability and SIMD use.
#   julia --project=<benchenv> codegen.jl
using DenoiseDNA, InteractiveUtils
using DenoiseDNA: VecBuffer, load!, fill_dp!
for W in (1, 16, 32)
    vb = VecBuffer{W}()
    ref = rand(0x01:0x04, 240)
    load!(vb, ref, [copy(ref) for _ in 1:W])
    fill_dp!(vb, 5, -4, -8, 0, 16)
    io = IOBuffer()
    code_warntype(io, fill_dp!, typeof.((vb, 5, -4, -8, 0, 16)))
    wt = String(take!(io))
    unstable = [l for l in split(wt, '\n') if occursin(r"::(Any|Union|Core\.Box)", l)]
    code_llvm(io, fill_dp!, typeof.((vb, 5, -4, -8, 0, 16)); debuginfo=:none)
    ll = String(take!(io))
    vec = Dict(t => count(t, ll) for t in ("<8 x i16>", "<16 x i16>", "<32 x i16>", "<64 x i8>", "<32 x i8>", "<16 x i8>"))
    println("W=$W  unstable lines: $(length(unstable))  vector types in IR: $vec")
    foreach(l -> println("    ", strip(l)), first(unstable, 8))
end
