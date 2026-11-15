# R's versions of the plots in plots.jl, for a side-by-side visual check.
#   Rscript plots_r.R <r_out_dir> <out_dir>
suppressMessages({library(dada2); library(ggplot2)})
a <- commandArgs(TRUE); rdir <- a[1]; out <- a[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
samples <- readLines(file.path(rdir, "samples.txt"))
fF <- file.path(rdir, "filtered", paste0(samples, "_F_filt.fastq.gz"))
raw <- file.path(dirname(rdir), "..", "MiSeq_SOP")
rawF <- file.path(raw, paste0(samples, "_R1_001.fastq"))
sv <- function(p, name) ggsave(file.path(out, name), p, width = 9, height = 6, dpi = 100)
sv(plotQualityProfile(rawF[1:4]), "R_quality_raw_F.png")
sv(plotQualityProfile(fF, aggregate = TRUE), "R_quality_filt_agg.png")
e <- learnErrors(fF[1:3], multithread = FALSE)
sv(plotErrors(e, nominalQ = TRUE), "R_errors_F.png")
sv(plotComplexity(rawF[1:2]), "R_complexity.png")
cat("R plots written to", out, "\n")
