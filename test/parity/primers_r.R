# R half of the removePrimers parity check.
#   Rscript primers_r.R <outdir>
# Then: julia --project=. test/parity/primers_julia.jl <outdir>
suppressMessages({library(dada2); library(Biostrings)})
out <- commandArgs(TRUE)[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f <- function(name) file.path(out, name)
set.seed(42)

# 1. The matcher itself, on random patterns/subjects with IUPAC codes.
iupac <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N")
rseq <- function(n, pamb) {
  x <- sample(c("A", "C", "G", "T"), n, replace = TRUE)
  amb <- runif(n) < pamb
  x[amb] <- sample(iupac, sum(amb), replace = TRUE)
  paste(x, collapse = "")
}
cases <- character(0)
for (k in 1:3000) {
  pat <- rseq(sample(2:14, 1), if (k %% 2 == 0) 0.15 else 0)
  subj <- rseq(sample(10:60, 1), if (k %% 3 == 0) 0.05 else 0)
  # plant a (mutated) copy of the pattern half of the time
  if (k %% 2 == 1 && nchar(subj) > nchar(pat)) {
    at <- sample(nchar(subj) - nchar(pat) + 1, 1)
    p <- strsplit(pat, "")[[1]]
    if (runif(1) < 0.5) p[sample(length(p), 1)] <- "A"
    if (runif(1) < 0.3) p <- p[-sample(length(p), 1)]
    substr(subj, at, at + length(p) - 1) <- paste(p, collapse = "")
  }
  mm <- sample(0:3, 1); indels <- k %% 4 == 0
  fixed <- dada2:::C_isACGT(pat)
  if (indels) {
    m <- matchPattern(pat, DNAString(subj), max.mismatch = mm, with.indels = TRUE, fixed = fixed)
  } else {
    m <- vmatchPattern(pat, DNAStringSet(subj), max.mismatch = mm, fixed = fixed)[[1]]
  }
  cases[k] <- sprintf("%s\t%s\t%d\t%d\t%d\t%s", pat, subj, mm, as.integer(indels), as.integer(fixed),
                      paste(sprintf("%d:%d", start(m), end(m)), collapse = ","))
}
writeLines(cases, f("match_cases.txt"))

# 2. removePrimers on the PacBio example (dada2 tutorial primers).
fn <- system.file("extdata", "samPBprimers.fastq.gz", package = "dada2")
writeLines(fn, f("input.txt"))
F27 <- "AGRGTTYGATYMTGGCTCAG"; R1492 <- "RGYTACCTTGTTACGACTT"
runs <- list(
  default  = list(),
  noorient = list(orient = FALSE),
  indels   = list(allow.indels = TRUE),
  notrim   = list(trim.fwd = FALSE, trim.rev = FALSE),
  fwdonly  = list(primer.rev = NULL),
  mm0      = list(max.mismatch = 0),
  mm5      = list(max.mismatch = 5, allow.indels = TRUE)
)
counts <- character(0)
for (nm in names(runs)) {
  args <- modifyList(list(fn = fn, fout = f(paste0("pb_", nm, ".fastq.gz")),
                          primer.fwd = F27, primer.rev = dada2:::rc(R1492)), runs[[nm]])
  if ("primer.rev" %in% names(runs[[nm]])) args$primer.rev <- NULL
  r <- do.call(removePrimers, args)
  counts <- c(counts, sprintf("%s\t%d\t%d", nm, r[1, 1], r[1, 2]))
}
writeLines(counts, f("pb_counts.txt"))
cat("R done\n")
