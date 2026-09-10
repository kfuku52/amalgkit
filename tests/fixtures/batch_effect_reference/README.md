# Batch correction reference fixtures

The input is a fixed synthetic 120-gene, 12-run count matrix. NumPy generator seed
2901 was used with independent Gaussian gene intercept/group/batch loadings;
means are `exp(4.5 + N(0,0.4) + N(0,0.6)*group + N(0,0.8)*batch)` and counts are
negative binomial with size 30. Group is 6+6 and batch alternates across runs.
The checked-in counts, rather than a regenerated random stream, are authoritative.

Run `Rscript generate.R` in this directory to regenerate the reference outputs.
Versions are recorded in `versions.txt`. No installed Python backend is used in
reference generation. The RUVr matrix method and residual function are read from
the pinned upstream revision; the SeqExpressionSet overload is not registered.
Its whole-number predicate is supplied locally and all fixture counts are integers.

- `sva_factors.tsv`: R sva IRW estimate, n.sv=1, B=5 on log1p counts.
- `sva_cleaned.tsv`: independent R calculation subtracting the factor space
  orthogonal to the protected design. This is AMALGKIT's exploratory protection
  contract, not a claim that sva itself exports this corrected matrix.
- `combat_counts.tsv`: R sva::ComBat_seq with group protection, default no shrinkage.
- `ruv_residuals.tsv`: edgeR first-pass GLM deviance residuals using the upstream
  RUVSeq residual function.
- `ruv_factors.tsv` / `ruv_counts.tsv`: upstream RUVr matrix kernel using these
  fixed residuals, all genes as controls, k=1, default rounding and pseudocount.

Numerical projection and RUVr kernel tests use tight tolerances / identical integer
outputs. End-to-end Python SVA and RUV fits are not asserted to equal R: Python's
local-FDR smoothing and RUV dispersion procedure differ. The ComBat comparison
requires the optional inmoose dependency. These small balanced fixtures do not
calibrate automatic k or evaluate unknown biological confounding.
