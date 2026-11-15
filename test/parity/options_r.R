# R half of the parity check for options not covered by the MiSeq SOP run
# (HANDOVER section B): filterAndTrim options, dada options, error functions,
# PacBio / variable-length reads, mergePairs options and chimera options.
#   Rscript options_r.R <outdir> <MiSeq_SOP dir> <results_nt1/r dir>
# Then: julia --project=. test/parity/options_julia.jl <outdir> <MiSeq_SOP dir>
suppressMessages({library(dada2); library(ShortRead)})
a <- commandArgs(TRUE)
out <- a[1]; sop <- a[2]; rres <- a[3]
dir.create(file.path(out, "filt"), showWarnings = FALSE, recursive = TRUE)
f <- function(x) file.path(out, x)
w <- function(x, name, ...) {
  # mergePairs(justConcatenate=TRUE) returns the sequences as a list column
  if (is.data.frame(x)) x[] <- lapply(x, function(col) if (is.list(col)) vapply(col, paste, "", collapse = "") else col)
  write.table(x, f(name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", ...)
}
wmat <- function(m, name) write.table(matrix(sprintf("%.17g", m), nrow(m)), f(name), sep = "\t", quote = FALSE,
                                      row.names = FALSE, col.names = FALSE)
ext <- system.file("extdata", package = "dada2")
writeLines(ext, f("extdata.txt"))
tm <- list()
tick <- function(nm, expr) { t0 <- Sys.time(); force(expr); tm[[nm]] <<- as.numeric(Sys.time() - t0, units = "secs") }

# ---------------------------------------------------------------- filterAndTrim
rawF <- file.path(sop, "F3D0_S188_L001_R1_001.fastq"); rawR <- file.path(sop, "F3D0_S188_L001_R2_001.fastq")
# a reverse file in a different order and missing some reads, for matchIDs
set.seed(3)
fqR <- readFastq(rawR)
keep <- sample(length(fqR), length(fqR) - 300)
unlink(f("R2_shuffled.fastq.gz"))
writeFastq(fqR[keep], f("R2_shuffled.fastq.gz"), compress = TRUE)
filt <- list(
  trim     = list(trimLeft = c(10, 5), trimRight = c(5, 10), maxEE = c(2, 3), minLen = 50),
  lens     = list(maxLen = c(245, 240), minLen = c(200, 150), minQ = c(10, 5), truncQ = c(11, 2), truncLen = c(0, 0)),
  minq     = list(minLen = c(200, 150), minQ = c(10, 5), truncQ = c(11, 2)),
  orient   = list(orient.fwd = "TACGG", truncLen = c(240, 160), maxEE = 2),
  matchids = list(matchIDs = TRUE, truncLen = c(240, 160), maxEE = 2, rev = f("R2_shuffled.fastq.gz")),
  idfield  = list(matchIDs = TRUE, id.field = 1, truncLen = c(240, 160), rev = f("R2_shuffled.fastq.gz")),
  maxn     = list(maxN = 1, truncQ = 0, rm.phix = FALSE),
  single   = list(single = TRUE, truncLen = 200, maxEE = 1),
  lowcomp  = list(single = TRUE, rm.lowcomplex = 8, truncLen = 150)
)
cnt <- character(0)
for (nm in names(filt)) {
  args <- filt[[nm]]
  rev <- if (is.null(args$rev)) rawR else args$rev
  single <- isTRUE(args$single)
  args$rev <- NULL; args$single <- NULL
  oF <- f(file.path("filt", paste0(nm, "_F.fastq.gz"))); oR <- f(file.path("filt", paste0(nm, "_R.fastq.gz")))
  r <- if (single) do.call(filterAndTrim, c(list(rawF, oF), args)) else
         do.call(filterAndTrim, c(list(rawF, oF, rev, oR), args))
  cnt <- c(cnt, sprintf("%s\t%d\t%d", nm, r[1, 1], r[1, 2]))
}
writeLines(cnt, f("filt_counts.txt"))
sr <- sread(readFastq(rawF))
writeLines(sprintf("%.17g", seqComplexity(sr)), f("complexity_k2.txt"))
writeLines(sprintf("%.17g", seqComplexity(sr, kmerSize = 3, window = 50, by = 7)), f("complexity_k3w50.txt"))

# ---------------------------------------------------------------- dada options
fnF <- file.path(ext, c("sam1F.fastq.gz", "sam2F.fastq.gz")); fnR <- file.path(ext, c("sam1R.fastq.gz", "sam2R.fastq.gz"))
errF <- learnErrors(fnF, multithread = FALSE, verbose = FALSE)
errR <- learnErrors(fnR, multithread = FALSE, verbose = FALSE)
wmat(getErrors(errF), "errF.tsv"); wmat(getErrors(errR), "errR.tsv"); wmat(errF$trans, "transF.tsv")
drpF <- derepFastq(fnF); drpR <- derepFastq(fnR)
base <- dada(drpF, err = errF, multithread = FALSE, verbose = FALSE)
priors <- getSequences(base[[2]])[1:20]
writeLines(priors, f("priors.txt"))
dopts <- list(
  default    = list(),
  omega      = list(OMEGA_A = 1e-20, OMEGA_C = 1e-10),
  singletons = list(DETECT_SINGLETONS = TRUE),
  priors     = list(priors = priors, OMEGA_P = 1e-2, OMEGA_A = 1e-60),
  homo       = list(HOMOPOLYMER_GAP_PENALTY = -1),
  novec      = list(VECTORIZED_ALIGNMENT = FALSE),
  band32     = list(BAND_SIZE = 32),
  band0      = list(BAND_SIZE = 0),
  unbanded   = list(BAND_SIZE = -1),
  minfold    = list(MIN_FOLD = 2, MIN_HAMMING = 2, MIN_ABUNDANCE = 3),
  maxclust   = list(MAX_CLUST = 5),
  nokmers    = list(USE_KMERS = FALSE),
  kdist      = list(KDIST_CUTOFF = 0.2),
  nogapless  = list(GAPLESS = FALSE),
  nogreedy   = list(GREEDY = FALSE),
  gaps       = list(GAP_PENALTY = -4L, MATCH = 4L, MISMATCH = -5L)
)
wdada <- function(dd, nm) for (i in seq_along(dd)) {
  w(data.frame(sequence = dd[[i]]$sequence, abundance = dd[[i]]$denoised), sprintf("dada_%s_%d.tsv", nm, i))
  writeLines(as.character(dd[[i]]$map), f(sprintf("dadamap_%s_%d.txt", nm, i)))
}
for (nm in names(dopts)) tick(paste0("dada_", nm), {
  dd <- do.call(dada, c(list(drpF, err = errF, multithread = FALSE, verbose = FALSE), dopts[[nm]]))
  wdada(dd, nm)
})
# USE_QUALS = FALSE: learn the no-quality error model, then denoise with it
tick("learn_noqual", errNQ <- learnErrors(fnF, USE_QUALS = FALSE, multithread = FALSE, verbose = FALSE))
wmat(getErrors(errNQ), "errNQ.tsv")
wdada(dada(drpF, err = errNQ, USE_QUALS = FALSE, multithread = FALSE, verbose = FALSE), "noqual")

# ---------------------------------------------------------------- error functions (isolated, on R's trans)
trans <- errF$trans
wmat(loessErrfun(trans), "errfun_loess.tsv")
wmat(noqualErrfun(trans), "errfun_noqual.tsv")
wmat(makeBinnedQualErrfun(c(2, 11, 25, 37, 40))(trans), "errfun_binned.tsv")
wmat(inflateErr(getErrors(errF), 3), "errfun_inflate.tsv")
wmat(inflateErr(getErrors(errF), 2, inflateSelfTransitions = TRUE), "errfun_inflate_self.tsv")

# ---------------------------------------------------------------- PacBio, variable-length reads
fnPB <- file.path(ext, "samPB.fastq.gz")
drpPB <- derepFastq(fnPB)
tick("learn_pb", errPB <- learnErrors(drpPB, errorEstimationFunction = PacBioErrfun, BAND_SIZE = 32,
                                      multithread = FALSE, verbose = FALSE))
wmat(getErrors(errPB), "errPB.tsv"); wmat(errPB$trans, "transPB.tsv")
wmat(PacBioErrfun(errPB$trans), "errfun_pacbio.tsv")
tick("dada_pb", ddPB <- dada(drpPB, err = errPB, BAND_SIZE = 32, multithread = FALSE, verbose = FALSE))
wdada(list(ddPB), "pacbio")

# ---------------------------------------------------------------- mergePairs options
ddR <- dada(drpR, err = errR, multithread = FALSE, verbose = FALSE)
mopts <- list(
  mismatch = list(maxMismatch = 2, minOverlap = 20, returnRejects = TRUE),
  concat   = list(justConcatenate = TRUE),
  overhang = list(trimOverhang = TRUE, returnRejects = TRUE),
  propcol  = list(propagateCol = c("n0", "birth_ham", "nunq"))
)
for (nm in names(mopts)) {
  mm <- do.call(mergePairs, c(list(base, drpF, ddR, drpR), mopts[[nm]]))
  for (i in seq_along(mm)) w(mm[[i]], sprintf("merge_%s_%d.tsv", nm, i))
}

# ---------------------------------------------------------------- chimera options (MiSeq SOP table)
cols <- read.delim(file.path(rres, "seqtab_cols.tsv"))$sequence
lg <- read.delim(file.path(rres, "seqtab.tsv"))
st <- matrix(0L, length(unique(lg$sample)), length(cols), dimnames = list(unique(lg$sample), cols))
st[cbind(match(lg$sample, rownames(st)), match(lg$sequence, cols))] <- lg$abundance
wl <- function(x, name) writeLines(as.character(as.integer(x)), f(name))
wl(isBimeraDenovo(st, allowOneOff = TRUE), "chim_oneoff.txt")
wl(isBimeraDenovo(st, allowOneOff = TRUE, minOneOffParentDistance = 2, maxShift = 8), "chim_oneoff2.txt")
wl(isBimeraDenovoTable(st, minSampleFraction = 0.5, ignoreNNegatives = 2), "chim_table.txt")
wl(isBimeraDenovoTable(st, allowOneOff = TRUE, minFoldParentOverAbundance = 1), "chim_table_oneoff.txt")
wl(colnames(st) %in% colnames(removeBimeraDenovo(st, method = "consensus", minSampleFraction = 0.6,
                                                 ignoreNNegatives = 0)), "chim_remove_kept.txt")

w(data.frame(step = names(tm), seconds = unlist(tm)), "timing.tsv")
cat("R done\n")
