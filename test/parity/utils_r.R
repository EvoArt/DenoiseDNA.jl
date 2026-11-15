# R half of the parity check for the small dada2 utilities (getSequences,
# uniquesToFasta, nwextract, strdiff, checkConvergence, sameOrder,
# getBadBases, isBadBaseFP, mergeSequenceTables, isShiftDenovo,
# collapseNoMismatch). Uses dada2's own inst/extdata.
#   Rscript utils_r.R <outdir>
# Then: julia --project=. test/parity/utils_julia.jl <outdir>
suppressMessages(library(dada2))
out <- commandArgs(TRUE)[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ext <- system.file("extdata", package = "dada2")
writeLines(ext, file.path(out, "extdata.txt"))
f <- function(name) file.path(out, name)
wl <- function(x, name) writeLines(as.character(x), f(name))
wnum <- function(x, name) writeLines(sprintf("%.17g", x), f(name))
wdf <- function(x, name) write.table(x, f(name), sep = "\t", quote = FALSE, row.names = FALSE)
wst <- function(st, name) {
  df <- data.frame(sample = rownames(st), st, check.names = FALSE)
  write.table(df, f(name), sep = "\t", quote = FALSE, row.names = FALSE)
}
set.seed(1)

fnF <- file.path(ext, c("sam1F.fastq.gz", "sam2F.fastq.gz"))
fnR <- file.path(ext, c("sam1R.fastq.gz", "sam2R.fastq.gz"))

# getSequences
wl(getSequences(file.path(ext, "example_seqs.fa")), "getseq_fa.txt")
wl(getSequences(fnF[1]), "getseq_fq.txt")

# sameOrder
wl(c(dada2:::sameOrder(fnF[1], fnR[1]), dada2:::sameOrder(fnF[1], fnR[2])), "sameorder.txt")

# learnErrors -> checkConvergence
errF <- learnErrors(fnF, multithread = FALSE, verbose = FALSE)
errR <- learnErrors(fnR, multithread = FALSE, verbose = FALSE)
wnum(dada2:::checkConvergence(errF), "checkconv_F.txt")
wnum(dada2:::checkConvergence(errR), "checkconv_R.txt")

# dada -> uniquesToFasta, getBadBases, isBadBaseFP
drpF <- derepFastq(fnF); drpR <- derepFastq(fnR)
ddF <- dada(drpF, err = errF, multithread = FALSE, verbose = FALSE)
ddR <- dada(drpR, err = errR, multithread = FALSE, verbose = FALSE)
for (i in 1:2) {
  u <- getUniques(ddF[[i]])
  wdf(data.frame(sequence = names(u), abundance = unname(u)), sprintf("uniques_%d.tsv", i))
  uniquesToFasta(ddF[[i]], f(sprintf("unq_%d.fa", i)))
  uniquesToFasta(ddF[[i]], f(sprintf("unq_%d_w60.fa", i)), width = 60)
  uniquesToFasta(ddF[[i]], f(sprintf("unq_%d_ids.fa", i)), ids = paste0("ASV", seq_along(u)))
}
uniquesToFasta(ddF[[1]], f("unq_append.fa"))
uniquesToFasta(ddF[[2]], f("unq_append.fa"), mode = "a")

bb <- list()
for (i in 1:2) {
  cl <- ddF[[i]]$clustering
  wdf(data.frame(sequence = cl$sequence, birth_ham = cl$birth_ham), sprintf("clust_%d.tsv", i))
  wdf(ddF[[i]]$birth_subs[, c("pos", "clust")], sprintf("bsubs_%d.tsv", i))
  for (omegaB in c(1e-20, 1e-10, 1e-3, 1)) for (minOcc in c(1, 2, 4)) {
    key <- sprintf("%d\t%g\t%d", i, omegaB, minOcc)
    b <- dada2:::getBadBases(cl, ddF[[i]]$birth_subs, omegaB = omegaB, minOccurence = minOcc)
    fp <- dada2:::isBadBaseFP(cl, ddF[[i]]$birth_subs, omegaB = omegaB, minOccurence = minOcc, minFraction = 0.3)
    bb[[length(bb) + 1]] <- sprintf("%s\t%s\t%s", key, paste(b, collapse = ","), paste(as.integer(fp), collapse = ","))
  }
}
wl(unlist(bb), "badbases.txt")

# mergePairs -> sequence tables -> mergeSequenceTables, collapseNoMismatch, isShiftDenovo
mm <- mergePairs(ddF, drpF, ddR, drpR)
st <- makeSequenceTable(mm); rownames(st) <- c("s1", "s2")
stF <- makeSequenceTable(ddF); rownames(stF) <- c("f1", "f2")
stRC <- stF; colnames(stRC) <- rc(colnames(stF)); rownames(stRC) <- c("r1", "r2")
cat("merged table:", dim(st), " forward table:", dim(stF), "\n")
stDup <- st[, seq_len(min(3, ncol(st))), drop = FALSE]; rownames(stDup) <- c("s2", "s3")
wst(st, "st.tsv"); wst(stF, "stF.tsv"); wst(stRC, "stRC.tsv"); wst(stDup, "stDup.tsv")
wst(mergeSequenceTables(st, stF), "merge_plain.tsv")
wst(mergeSequenceTables(st, stF, stRC, tryRC = TRUE), "merge_rc.tsv")
wst(mergeSequenceTables(st, stDup, repeats = "sum", orderBy = "nsamples"), "merge_sum.tsv")
wst(mergeSequenceTables(tables = list(stF, st), orderBy = NULL), "merge_noorder.tsv")
wst(collapseNoMismatch(st), "collapse.tsv")

unqs <- c(colSums(st), colSums(stF))
unqs <- unqs[!duplicated(names(unqs))]
top <- names(sort(unqs, decreasing = TRUE))[1:8]
shifted <- c(substring(top[1:4], 4), paste0("ACG", top[5:6]), substring(top[7:8], 1, nchar(top[7:8]) - 5))
syn <- c(unqs, setNames(seq_along(shifted), shifted))
syn <- setNames(as.integer(syn), names(syn))
wdf(data.frame(sequence = names(syn), abundance = unname(syn)), "shift_input.tsv")
wl(as.integer(isShiftDenovo(syn)), "shift_default.txt")
wl(as.integer(isShiftDenovo(syn, minOverlap = 200)), "shift_200.txt")

# nwextract, strdiff
pairs <- list()
for (k in 1:6) {
  s <- top[k]
  q <- substring(s, 20 + k, 150 + 10 * k)
  chars <- strsplit(s, "")[[1]]
  pos <- sample(length(chars), k)
  chars[pos] <- sapply(chars[pos], function(x) setdiff(c("A", "C", "G", "T"), x)[1])
  m <- paste(chars, collapse = "")
  d <- dada2:::strdiff(s, m)
  pairs[[k]] <- sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s", s, q, m, dada2:::nwextract(q, m),
                        paste(d$pos, collapse = ","), paste(d$nt0, collapse = ","), paste(d$nt1, collapse = ","))
}
wl(unlist(pairs), "nwextract_strdiff.txt")
cat("R done\n")
