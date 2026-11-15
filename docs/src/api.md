# API

```@meta
CurrentModule = DenoiseDNA
```

```@docs
DenoiseDNA
```

## Filtering and trimming

```@docs
filterAndTrim
fastqFilter
fastqPairedFilter
removePrimers
isPhiX
seqComplexity
```

## Reading and dereplication

```@docs
read_fastq
Derep
derepFastq
combineDereps
getUniques
getSequences
uniquesToFasta
```

## Error models

```@docs
learnErrors
ErrorModel
getErrors
loessErrfun
noqualErrfun
makeBinnedQualErrfun
PacBioErrfun
inflateErr
accumulateTrans
checkConvergence
```

## Sample inference

```@docs
dada
DadaResult
DadaOptions
setDadaOpt
getDadaOpt
getBadBases
isBadBaseFP
```

## Merging and alignment

```@docs
mergePairs
nwalign
nwhamming
rc
```

## Sequence tables

```@docs
SequenceTable
makeSequenceTable
mergeSequenceTables
collapseNoMismatch
```

## Chimeras

```@docs
removeBimeraDenovo
isBimeraDenovo
isBimeraDenovoTable
isBimera
isShiftDenovo
isShift
isShiftedPair
```

## Plots

The plotting functions need a Makie backend, e.g. `using CairoMakie`.
`qualityProfile` and `errorProfile` return the plotted data and need no
backend.

```@docs
plotQualityProfile
plotErrors
plotComplexity
qualityProfile
errorProfile
```

## Utilities

Not exported.

```@docs
DenoiseDNA.setprogress!
DenoiseDNA.read_fasta
DenoiseDNA.nwextract
DenoiseDNA.strdiff
DenoiseDNA.sameOrder
```

## Index

```@index
Pages = ["api.md"]
```
