"""
    DenoiseDNA

The DADA2 amplicon denoising pipeline: [`filterAndTrim`](@ref),
[`learnErrors`](@ref), [`dada`](@ref), [`mergePairs`](@ref),
[`makeSequenceTable`](@ref) and [`removeBimeraDenovo`](@ref), with the names,
options and results of the R package `dada2`.
"""
module DenoiseDNA

using DataFrames
using CodecZlib: GzipDecompressor, GzipCompressorStream
using Rmath_jll: libRmath
using LinearAlgebra
using SIMD: Vec, vifelse, vload, vstore
import Random

include("style.jl")
include("progress.jl")
include("fastq.jl")
include("filter.jl")
include("derep.jl")
include("kmers.jl")
include("align_vec.jl")
include("align.jl")
include("rmath.jl")
include("loess.jl")
include(joinpath("dada", "cluster.jl"))
include(joinpath("dada", "dada.jl"))
include("errors.jl")
include("merge.jl")
include("seqtab.jl")
include("chimera.jl")
include("utils.jl")
include("primers.jl")
include("plotting.jl")
include("show.jl")

export filterAndTrim, fastqFilter, fastqPairedFilter, removePrimers, isPhiX, seqComplexity,
       read_fastq, Derep, derepFastq, combineDereps,
       DadaOptions, DadaResult, dada, setDadaOpt, getDadaOpt,
       getUniques, getSequences, uniquesToFasta,
       ErrorModel, learnErrors, getErrors, loessErrfun, noqualErrfun, makeBinnedQualErrfun,
       PacBioErrfun, inflateErr, accumulateTrans, checkConvergence, getBadBases, isBadBaseFP,
       mergePairs, nwalign, nwhamming, rc,
       SequenceTable, makeSequenceTable, collapseNoMismatch, mergeSequenceTables,
       isBimera, isBimeraDenovo, isBimeraDenovoTable, removeBimeraDenovo,
       isShiftDenovo, isShift, isShiftedPair,
       qualityProfile, errorProfile, plotQualityProfile, plotErrors, plotComplexity

end
