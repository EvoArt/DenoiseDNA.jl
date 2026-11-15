using DenoiseDNA
using DataFrames
using Random
using Test

const D = DenoiseDNA
const FIXTURES = joinpath(@__DIR__, "fixtures")
isdefined(DenoiseDNA, :setprogress!) && DenoiseDNA.setprogress!(false)

include("helpers.jl")

@testset "DenoiseDNA" begin
    include("test_fastq.jl")
    include("test_filter.jl")
    include("test_derep.jl")
    include("test_kmers.jl")
    include("test_align.jl")
    include("test_loess.jl")
    include("test_dada.jl")
    include("test_merge_seqtab.jl")
    include("test_chimera.jl")
    include("test_utils.jl")
    include("test_plotting.jl")
    isdir(FIXTURES) && include("test_r_fixtures.jl")
end
