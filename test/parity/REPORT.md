# DenoiseDNA vs R dada2 parity report

R: 1.38.0, R version 4.5.3 (2026-03-11)
Samples: 20

## End-to-end (Julia from raw reads)

- **filterAndTrim**: counts identical: yes; filtered FASTQ records (id, sequence, quality) identical in 20/20 samples
- **derepFastq**: uniques+order identical 40/40; mean qualities identical 40/40; read map identical 40/40
- **learnErrors (F)**: rounds R=6 J=6; transition counts identical: yes; error rates max rel diff 6.17e-15
- **learnErrors (R)**: rounds R=6 J=6; transition counts identical: yes; error rates max rel diff 3.28e-14
- **dada**: ASVs R=4180 J=4180 shared=4180; identical ASV lists 40/40; identical abundances 40/40; identical unique->ASV map 40/40; cluster p-value max rel diff 1.99e-12
- **dada pool=TRUE**: ASVs R=6484 J=6484 shared=6484; identical ASV lists 20/20; identical abundances 20/20; identical unique->ASV map 20/20; cluster p-value max rel diff 5.24e-13
- **dada pool="pseudo"**: ASVs R=3582 J=3582 shared=3582; identical ASV lists 20/20; identical abundances 20/20; identical unique->ASV map 20/20; cluster p-value max rel diff 7.35e-13
- **mergePairs** (incl. rejects): identical tables 20/20
- **makeSequenceTable**: entries identical: yes; column order identical: yes (293 vs 293 sequences)
- **removeBimeraDenovo (consensus)**: kept R=232 J=232; identical: yes
- **removeBimeraDenovo (pooled)**: kept R=229 J=229; identical: yes
- **removeBimeraDenovo (per-sample)**: entries R=1641 J=1641; identical: yes

## Isolated (each step on R's previous output)

- **derepFastq**: uniques+order identical 40/40; mean qualities identical 40/40; read map identical 40/40
- **learnErrors (F)**: rounds R=6 J=6; transition counts identical: yes; error rates max rel diff 6.17e-15
- **learnErrors (R)**: rounds R=6 J=6; transition counts identical: yes; error rates max rel diff 3.28e-14
- **dada**: ASVs R=4180 J=4180 shared=4180; identical ASV lists 40/40; identical abundances 40/40; identical unique->ASV map 40/40; cluster p-value max rel diff 4.65e-16
- **dada pool=TRUE**: ASVs R=6484 J=6484 shared=6484; identical ASV lists 20/20; identical abundances 20/20; identical unique->ASV map 20/20; cluster p-value max rel diff 4.11e-16
- **dada pool="pseudo"**: ASVs R=3582 J=3582 shared=3582; identical ASV lists 20/20; identical abundances 20/20; identical unique->ASV map 20/20; cluster p-value max rel diff 7.35e-13
- **mergePairs** (incl. rejects): identical tables 20/20
- **makeSequenceTable**: entries identical: yes; column order identical: yes (293 vs 293 sequences)
- **removeBimeraDenovo (consensus)**: kept R=232 J=232; identical: yes
- **removeBimeraDenovo (pooled)**: kept R=229 J=229; identical: yes
- **removeBimeraDenovo (per-sample)**: entries R=1641 J=1641; identical: yes

## Timing (seconds, end-to-end)

| step | R | Julia | R/Julia |
|---|---:|---:|---:|
| filter | 31.06 | 14.92 | 2.08x |
| derep | 2.60 | 1.19 | 2.18x |
| learn_F | 125.94 | 18.26 | 6.90x |
| learn_R | 74.97 | 10.81 | 6.93x |
| dada_F | 22.42 | 6.28 | 3.57x |
| dada_R | 13.22 | 3.32 | 3.98x |
| dada_pool | 48.10 | 14.26 | 3.37x |
| dada_pseudo | 56.54 | 12.31 | 4.59x |
| merge | 1.88 | 2.20 | 0.85x |
| seqtab | 0.00 | 0.36 | 0.01x |
| nochim | 0.53 | 0.73 | 0.73x |
| **total** | 377.26 | 84.64 | 4.46x |
