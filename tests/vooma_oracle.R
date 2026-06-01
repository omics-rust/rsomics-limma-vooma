#!/usr/bin/env Rscript
# vooma oracle: read a log-expression matrix TSV (header = sample ids, col 1 =
# gene ids) and a design matrix TSV (header = coefficient names, col 1 = sample
# ids), run limma vooma, and write the precision weights matrix
# (gene, then one column per sample).
#
# Usage: vooma_oracle.R <expr.tsv> <design.tsv> <out.tsv>
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
expr_path <- args[1]
design_path <- args[2]
out_path <- args[3]

E <- as.matrix(read.delim(expr_path, row.names = 1, check.names = FALSE))
design <- as.matrix(read.delim(design_path, row.names = 1, check.names = FALSE))

v <- vooma(E, design)
W <- v$weights

con <- file(out_path, "w")
writeLines(paste(c("gene", colnames(E)), collapse = "\t"), con)
g <- function(x) formatC(x, digits = 10, format = "g", flag = "")
for (i in seq_len(nrow(W))) {
  writeLines(paste(c(rownames(E)[i], g(W[i, ])), collapse = "\t"), con)
}
close(con)
