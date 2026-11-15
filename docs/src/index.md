# DenoiseDNA.jl

Exact amplicon sequence variants from Illumina and PacBio amplicon reads: the
[DADA2](https://benjjneb.github.io/dada2/) pipeline in Julia.

DenoiseDNA reimplements the R package `dada2` (version 1.38). It uses the
same function and option names and gives the same results: on the MiSeq SOP
data every step, from `filterAndTrim` to `removeBimeraDenovo`, returns the
same sequences, abundances and tables as R. It is about 4 times faster on
one thread.

## Installation

```julia
using Pkg
Pkg.add("DenoiseDNA")
```

Plots need a Makie backend, e.g. `Pkg.add("CairoMakie")`. Taxonomic
assignment is in [AssignTaxonomy.jl](https://github.com/EvoArt/AssignTaxonomy.jl).

## Workflow

```julia
using DenoiseDNA

out    = filterAndTrim(fnFs, filtFs, fnRs, filtRs; truncLen=(240, 160), maxEE=(2, 2))
errF   = learnErrors(filtFs; multithread=true)
errR   = learnErrors(filtRs; multithread=true)
dadaFs = dada(filtFs, errF; multithread=true)
dadaRs = dada(filtRs, errR; multithread=true)
merged = mergePairs(dadaFs, filtFs, dadaRs, filtRs)
seqtab = makeSequenceTable(merged; names=sample_names)
nochim = removeBimeraDenovo(seqtab; method="consensus", multithread=true)
```

The [Tutorial](@ref "Tutorial: the DADA2 workflow") goes through it on real
data. The [API](@ref) lists every function.

## Differences from R

- Arguments with dots in R use underscores: `rm.phix` is `rm_phix`,
  `id.sep` is `id_sep`.
- Per-direction arguments are tuples: `truncLen=(240, 160)`.
- `pool="pseudo"` can also be written `pool=:pseudo`.
- Results are Julia types ([`Derep`](@ref), [`DadaResult`](@ref),
  [`ErrorModel`](@ref), [`SequenceTable`](@ref)) and `DataFrame`s rather than
  R lists and named vectors.
- Non-Phred+33 FASTQ (`qualityType`) is not supported.

## Citation

If you use this package, please cite the DADA2 paper:

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
DADA2: High-resolution sample inference from Illumina amplicon data.
*Nature Methods* 13, 581–583.
[doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
