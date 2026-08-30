# Independent fixed-reference TMM factors; no AMALGKIT implementation is sourced.
args <- commandArgs(trailingOnly=TRUE)
root <- args[[1]]
suppressPackageStartupMessages(library(edgeR))
write(as.character(packageVersion("edgeR")), file=file.path(root, "version.txt"))
for (case in args[-1]) {
    x <- as.matrix(read.delim(file.path(root, paste0("counts-", case, ".tsv")), header=FALSE))
    libs <- scan(file.path(root, paste0("libs-", case, ".txt")), quiet=TRUE)
    factors <- calcNormFactors(x, lib.size=libs, refColumn=1, method="TMM")
    write(sprintf("%.17g", factors), file=file.path(root, paste0("factors-", case, ".txt")), ncolumns=1)
}
