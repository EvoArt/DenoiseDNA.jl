module DenoiseDNA
using FASTX, BioSequences, LoopVectorization, CodecZlib, BioAlignments, DataFrames#, CSV,ProfileView
using LazyStack, Poppler_jll,StatsBase#,FastHistograms,StatsBase,ProfileView#,CairoMakie
using RCall#, ReadDataStores
#using StatsFuns, SortingLab, Tullio,Polyester,UnicodePlots
using StatsFuns, Loess,CairoMakie#, MultipleTesting, SparseArrays

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
