# Tutorial: the DADA2 workflow

This tutorial follows the
[DADA2 pipeline tutorial](https://benjjneb.github.io/dada2/tutorial.html)
step by step, in Julia. It turns paired-end Illumina MiSeq reads into a table
of amplicon sequence variants (ASVs) and assigns them a taxonomy. Function
names and options are those of the R package, so the R documentation applies
too.

## Setup

Start Julia with several threads (`julia -t auto`); steps called with
`multithread=true` use all of them.

```julia
using Pkg
Pkg.add(["DenoiseDNA", "AssignTaxonomy", "DataFrames", "CairoMakie"])

using DenoiseDNA, DataFrames
```

## Data

The data are the 2x250 V4 16S rRNA reads of the
[mothur MiSeq SOP](https://mothur.org/wiki/miseq_sop/): gut samples from one
mouse, taken at time points after weaning, plus a mock community of 20
known strains. Download
[miseqsopdata.zip](https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip)
and unzip it.

```julia
path = "MiSeq_SOP"
fnFs = sort(filter(endswith("_R1_001.fastq"), readdir(path; join=true)))
fnRs = sort(filter(endswith("_R2_001.fastq"), readdir(path; join=true)))
sample_names = first.(split.(basename.(fnFs), "_"))
```

Files are named `SAMPLENAME_XXX_R1_001.fastq` (forward) and `..._R2_001.fastq`
(reverse), so the sample name is everything before the first underscore.

## Inspect read quality

```julia
using CairoMakie
plotQualityProfile(fnFs[1:2])
```

Each panel is one file: a heatmap of quality scores at each cycle, with the
mean (green line) and the median and quartiles (orange lines). The forward
reads stay good to the end; trimming the last 10 cycles is a safe choice, so
they are truncated at 240.

```julia
plotQualityProfile(fnRs[1:2])
```

The reverse reads are worse, as is usual on the MiSeq, and quality falls
around cycle 160, so they are truncated there. The reads must still overlap
after truncation: this V4 amplicon is about 250 bp, and 240 + 160 leaves
plenty of overlap.

## Filter and trim

```julia
filtFs = joinpath.(path, "filtered", sample_names .* "_F_filt.fastq.gz")
filtRs = joinpath.(path, "filtered", sample_names .* "_R_filt.fastq.gz")

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs; truncLen=(240, 160), maxN=0,
                    maxEE=(2, 2), truncQ=2, rm_phix=true, compress=true,
                    multithread=true)
```

Per-direction settings are given as `(forward, reverse)`. `maxEE` sets the
number of expected errors allowed in a read, a better filter than average
quality. Reads with any `N` are dropped (DADA2 needs none), as are reads that
match the phiX genome. `out` has the reads in and out for each sample.

## Learn the error rates

```julia
errF = learnErrors(filtFs; multithread=true)
errR = learnErrors(filtRs; multithread=true)
```

`learnErrors` alternates between inferring the sample composition and
estimating the error rates until the two agree. It uses the first `1e8`
bases of the data by default (`nbases`).

```julia
plotErrors(errF; nominalQ=true)
```

Each panel shows one transition (e.g. `A2C`): the observed error frequencies
(points), the fitted rates (black line) and the rates implied by the quality
scores (red line). The fit should follow the points, and error rates should
fall as quality rises.

## Sample inference

```julia
dadaFs = dada(filtFs, errF; multithread=true)
dadaRs = dada(filtRs, errR; multithread=true)
dadaFs[1]
```

Each [`DadaResult`](@ref) holds the denoised sequences and their abundances.
By default samples are processed independently; `pool=true` pools all of
them, and `pool=:pseudo` shares information between samples at a fraction of
the cost.

## Merge paired reads

```julia
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs; verbose=true)
first(mergers[1], 6)
```

Each forward read is merged with its reverse read when they overlap by at
least 12 bases with no mismatch. The result for each sample is a `DataFrame`
of merged sequences, with abundances and the indices of the forward and
reverse ASVs. Pairs that fail to merge are dropped (`returnRejects=true`
keeps them, marked `accept = false`).

## Sequence table

```julia
seqtab = makeSequenceTable(mergers; names=sample_names)
size(seqtab)        # (20, 293)
```

A [`SequenceTable`](@ref) has a row per sample and a column per ASV;
`seqtab.table` holds the counts. The lengths of the merged sequences:

```julia
lens = length.(seqtab.sequences)
[l => count(==(l), lens) for l in sort(unique(lens))]
```

Nearly all fall within the expected range of the V4 amplicon. Sequences much
longer or shorter can come from non-specific priming and can be removed:

```julia
seqtab2 = seqtab[:, 250 .<= lens .<= 256]
```

## Remove chimeras

```julia
seqtab_nochim = removeBimeraDenovo(seqtab; method="consensus", multithread=true, verbose=true)
size(seqtab_nochim)                       # (20, 232)
sum(seqtab_nochim.table) / sum(seqtab.table)
```

Chimeras are sequences made of a left part of one more abundant sequence and
a right part of another. Here they are a large share of the ASVs but only a
few percent of the reads.

## Track reads through the pipeline

```julia
track = DataFrame(sample = sample_names,
                  input = out.reads_in,
                  filtered = out.reads_out,
                  denoisedF = [sum(d.abundance) for d in dadaFs],
                  denoisedR = [sum(d.abundance) for d in dadaRs],
                  merged = [sum(m.abundance) for m in mergers],
                  nonchim = vec(sum(seqtab_nochim.table; dims=2)))
```

No step should lose most of the reads. A large loss at merging usually means
the reads were truncated too short to overlap; at chimera removal, that
primers were left on the reads.

## Assign taxonomy

Taxonomy is assigned with
[AssignTaxonomy.jl](https://github.com/EvoArt/AssignTaxonomy.jl), a port of
R's `assignTaxonomy` (the naive Bayesian classifier of Wang et al. 2007).
Download the Silva reference files formatted for DADA2 (see
[the DADA2 reference databases](https://benjjneb.github.io/dada2/training.html)):
`silva_nr99_v138.2_toGenus_trainset.fa.gz` and
`silva_v138.2_assignSpecies.fa.gz`.

```julia
using AssignTaxonomy

taxa = assignTaxonomy(seqtab_nochim, "silva_nr99_v138.2_toGenus_trainset.fa.gz")
taxa = addSpecies(taxa, "silva_v138.2_assignSpecies.fa.gz")
DataFrame(taxa)
```

Each row is an ASV, with a column per rank; ranks without enough bootstrap
support (`minBoot=50`) are `missing`. Species are assigned only by exact
matching, which is the only safe way to assign species from short amplicons.

## Evaluate accuracy on the mock community

The sample `Mock` contains 20 known strains, whose sequences are in
`HMP_MOCK.v35.fasta`.

```julia
i = findfirst(==("Mock"), seqtab_nochim.samples)
mock = seqtab_nochim.sequences[seqtab_nochim.table[i, :] .> 0]
mock_ref = getSequences(joinpath(path, "HMP_MOCK.v35.fasta"))

length(mock)                                    # ASVs found in the mock
count(s -> any(r -> occursin(s, r), mock_ref), mock)   # exact matches
```

All 20 strains are recovered exactly, with no spurious variants.

## Next steps

The sequence table, taxonomy table and sample data are what most community
analyses start from. `DataFrame(seqtab_nochim)` gives the counts in long
form (`sample`, `sequence`, `abundance`), and `DataFrame(taxa)` the
taxonomy, ready to join on the sequence.
