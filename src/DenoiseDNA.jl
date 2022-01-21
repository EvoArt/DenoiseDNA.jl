module DenoiseDNA
using FASTX, BioSequences,  BioAlignments # BioJulia
import BioAlignments.pairalign
using StatsFuns, Loess,LazyStack, StatsBase,LoopVectorization, CodecZlib,DataFrames
using CairoMakie,UnicodePlots, Poppler_jll # visualisation
using RCall

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
        get_files,
        learnErrors
end
