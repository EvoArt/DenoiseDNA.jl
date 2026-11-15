# Reference run of the dada2 MiSeq SOP pipeline in R.
# Dumps every intermediate at full precision for test/parity/compare.jl.
#
#   Rscript run_r.R <data_dir> <out_dir> [nthreads]
#
# <data_dir> holds the raw *_R1_001.fastq / *_R2_001.fastq files.

suppressMessages(library(dada2))
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]; out <- args[2]
nthreads <- if (length(args) >= 3) as.integer(args[3]) else 1L
mt <- if (nthreads > 1) nthreads else FALSE
dir.create(out, recursive = TRUE, showWarnings = FALSE)
filt_dir <- file.path(out, "filtered"); dir.create(filt_dir, showWarnings = FALSE)

timings <- list()
timeit <- function(label, expr) {
  t0 <- Sys.time(); val <- force(expr)
  timings[[label]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[time] %-10s %8.2f s\n", label, timings[[label]])); val
}
g17 <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
wtsv <- function(df, name) write.table(df, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)

fnFs <- sort(list.files(data_dir, pattern = "_R1_001.fastq$", full.names = TRUE), method = "radix")
fnRs <- sort(list.files(data_dir, pattern = "_R2_001.fastq$", full.names = TRUE), method = "radix")
n_env <- Sys.getenv("N_SAMPLES")
if (nzchar(n_env)) { k <- as.integer(n_env); fnFs <- fnFs[1:k]; fnRs <- fnRs[1:k] }
samples <- gsub("_R1_001.fastq", "", basename(fnFs))
writeLines(samples, file.path(out, "samples.txt"))
filtFs <- file.path(filt_dir, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(samples, "_R_filt.fastq.gz"))

# 1. filter
ft <- timeit("filter", filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160), maxN = 0,
                                     maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, compress = TRUE,
                                     multithread = mt))
wtsv(data.frame(sample = samples, reads_in = ft[, 1], reads_out = ft[, 2]), "filter.tsv")

# 2. derep
drF <- timeit("derep", derepFastq(filtFs)); drR <- derepFastq(filtRs)
if (length(samples) == 1) { drF <- list(drF); drR <- list(drR) }
names(drF) <- names(drR) <- samples
dump_derep <- function(d, tag) {
  wtsv(data.frame(sequence = names(d$uniques), abundance = as.integer(d$uniques)), paste0("derep_", tag, ".tsv"))
  q <- d$quals; q <- matrix(g17(q), nrow = nrow(q))
  write.table(q, file.path(out, paste0("derepq_", tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(as.character(d$map), file.path(out, paste0("derepmap_", tag, ".txt")))
}
for (s in samples) { dump_derep(drF[[s]], paste0("F_", s)); dump_derep(drR[[s]], paste0("R_", s)) }

# 3. learn errors
dump_err <- function(e, tag) {
  m <- getErrors(e)
  write.table(matrix(g17(m), nrow = 16), file.path(out, paste0("err_", tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(e$trans, file.path(out, paste0("errtrans_", tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(as.character(length(e$err_in)), file.path(out, paste0("errnconsist_", tag, ".txt")))
}
errF <- timeit("learn_F", learnErrors(filtFs, multithread = mt))
errR <- timeit("learn_R", learnErrors(filtRs, multithread = mt))
dump_err(errF, "F"); dump_err(errR, "R")

# 4. dada
dump_dada <- function(dd, tag) {
  cl <- dd$clustering
  for (c in c("pval", "birth_pval", "birth_fold", "birth_qave")) cl[[c]] <- g17(cl[[c]])
  wtsv(cl, paste0("dada_", tag, ".tsv"))
  writeLines(ifelse(is.na(dd$map), "0", as.character(dd$map)), file.path(out, paste0("dadamap_", tag, ".txt")))
  writeLines(g17(dd$pval), file.path(out, paste0("dadapval_", tag, ".txt")))
}
ddF <- timeit("dada_F", dada(drF, err = errF, multithread = mt))
ddR <- timeit("dada_R", dada(drR, err = errR, multithread = mt))
if (length(samples) == 1) { ddF <- list(ddF); ddR <- list(ddR) }
names(ddF) <- names(ddR) <- samples
for (s in samples) { dump_dada(ddF[[s]], paste0("F_", s)); dump_dada(ddR[[s]], paste0("R_", s)) }

# 4b. dada with the alternative pooling modes (forward reads only)
ddFp <- timeit("dada_pool", dada(drF, err = errF, pool = TRUE, multithread = mt))
ddFs <- timeit("dada_pseudo", dada(drF, err = errF, pool = "pseudo", multithread = mt))
if (length(samples) == 1) { ddFp <- list(ddFp); ddFs <- list(ddFs) }
names(ddFp) <- names(ddFs) <- samples
for (s in samples) { dump_dada(ddFp[[s]], paste0("Fpool_", s)); dump_dada(ddFs[[s]], paste0("Fpseudo_", s)) }

# 5. merge
mg <- timeit("merge", mergePairs(ddF, drF, ddR, drR, returnRejects = TRUE))
if (length(samples) == 1) mg <- list(mg)
names(mg) <- samples
for (s in samples) wtsv(mg[[s]], paste0("merge_", s, ".tsv"))
mg_acc <- lapply(mg, function(m) m[m$accept, ])

# 6-7. seqtab + chimeras (all three methods)
st <- timeit("seqtab", makeSequenceTable(mg_acc))
long <- function(tab) {
  rows <- which(tab > 0, arr.ind = TRUE)
  data.frame(sample = rownames(tab)[rows[, 1]], sequence = colnames(tab)[rows[, 2]], abundance = tab[rows])
}
wtsv(data.frame(sequence = colnames(st)), "seqtab_cols.tsv")
wtsv(long(st), "seqtab.tsv")
nc <- timeit("nochim", removeBimeraDenovo(st, method = "consensus", multithread = mt))
wtsv(data.frame(sequence = colnames(nc)), "nochim_consensus.tsv")
ncp <- removeBimeraDenovo(st, method = "pooled", multithread = mt)
wtsv(data.frame(sequence = colnames(ncp)), "nochim_pooled.tsv")
ncs <- removeBimeraDenovo(st, method = "per-sample", multithread = mt)
wtsv(long(ncs), "nochim_persample.tsv")

wtsv(data.frame(step = names(timings), seconds = unlist(timings)), "timing.tsv")
writeLines(c(as.character(packageVersion("dada2")), R.version.string), file.path(out, "versions.txt"))
cat("R pipeline complete\n")
