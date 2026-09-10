# Run from this directory with Rscript generate.R. Inputs are fixed fixtures,
# not generated anew; package versions and the RUVr source revision are recorded.
counts <- as.matrix(read.delim("counts.tsv", row.names=1, check.names=FALSE))
metadata <- read.delim("metadata.tsv")
X <- model.matrix(~factor(sample_group), metadata)
Y <- log1p(counts)
sv <- sva::irwsva.build(Y, X, n.sv=1, B=5)$sv
write.table(sv, "sva_factors.tsv", sep="\t", quote=FALSE, col.names=NA)

# Independent R expression of the design-protected removal contract.
residual_sv <- qr.resid(qr(X), sv)
Q <- qr.Q(qr(residual_sv))
cleaned <- Y - t(Q %*% crossprod(Q, t(Y)))
write.table(cleaned, "sva_cleaned.tsv", sep="\t", quote=FALSE, col.names=NA)

combat <- sva::ComBat_seq(counts, batch=metadata$bioproject, group=metadata$sample_group)
write.table(combat, "combat_counts.tsv", sep="\t", quote=FALSE, col.names=NA)

# Evaluate only the upstream RUVr matrix-method definition. The SeqExpressionSet
# overload is not registered and no RUVSeq/EDASeq package installation is needed.
revision <- "f560e84eb493f0c2d72279612f9e4ca99493b434"
source_url <- paste0("https://raw.githubusercontent.com/drisso/RUVSeq/", revision, "/R/RUVr-methods.R")
environment <- new.env()
environment$setMethod <- function(f, signature, definition) {
  if (identical(unname(signature[["x"]]), "matrix")) environment$ruvr <- definition
}
environment$.isWholeNumber <- function(x) abs(x - round(x)) < .Machine$double.eps^0.5
eval(parse(text=readLines(source_url, warn=FALSE)), envir=environment)
residual_url <- paste0("https://raw.githubusercontent.com/drisso/RUVSeq/", revision, "/R/residuals.DGEGLM.R")
environment$negative.binomial <- MASS::negative.binomial
eval(parse(text=readLines(residual_url, warn=FALSE)), envir=environment)
edge <- edgeR::DGEList(counts=counts)
edge <- edgeR::calcNormFactors(edge, method="upperquartile")
edge <- edgeR::estimateDisp(edge, design=X)
fit <- edgeR::glmFit(edge, design=X)
residuals <- environment$residuals.DGEGLM(fit, type="deviance")
ruv <- environment$ruvr(counts, seq_len(nrow(counts)), 1, residuals)
write.table(residuals, "ruv_residuals.tsv", sep="\t", quote=FALSE, col.names=NA)
write.table(ruv$W, "ruv_factors.tsv", sep="\t", quote=FALSE, col.names=NA)
write.table(ruv$normalizedCounts, "ruv_counts.tsv", sep="\t", quote=FALSE, col.names=NA)
writeLines(c(R.version.string, paste("sva", packageVersion("sva")),
             paste("edgeR", packageVersion("edgeR")), paste("RUVr source", source_url)), "versions.txt")
