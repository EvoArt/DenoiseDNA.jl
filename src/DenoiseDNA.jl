module DenoiseDNA
using FASTX, BioSequences, LoopVectorization, CodecZlib, BioAlignments, DataFrames
using LazyStack, Poppler_jll,StatsBase
using RCall
using StatsFuns, Loess,CairoMakie

include.(["plotting.jl",
        "fileIO.jl",
        "kmers.jl",
        "DADA.jl"])

# Write your package code here.
export quality_plot,
        dada,
        custering,
        trim,
        dereplicate,
        get_files
end
