# DenoiseDNA.jl

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://EvoArt.github.io/DenoiseDNA.jl/stable)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://EvoArt.github.io/DenoiseDNA.jl/dev)
[![CI](https://github.com/EvoArt/DenoiseDNA.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EvoArt/DenoiseDNA.jl/actions/workflows/CI.yml)

Exact amplicon sequence variants from amplicon sequencing reads: the
[DADA2](https://benjjneb.github.io/dada2/) pipeline in Julia. Same functions,
options and results as the R package `dada2` 1.38, about 4 times faster.

```julia
using Pkg
Pkg.add("DenoiseDNA")
```

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

See the [tutorial](https://EvoArt.github.io/DenoiseDNA.jl/stable/tutorial/)
for the full workflow. Taxonomy is assigned with
[AssignTaxonomy.jl](https://github.com/EvoArt/AssignTaxonomy.jl).

Please cite the DADA2 paper: Callahan et al. (2016), *Nature Methods* 13,
581–583, [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869).
