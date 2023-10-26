
get_singlecopy_bool_index = function(df_gc, spp_filled, percent_singlecopy_threshold=50) {
  is_ge_singlecopy_threshold = function(x, num_sp, percent_singlecopy_threshold) {
    num_singlecopy_species = sum(x==1)
    percent_singlecopy_species = num_singlecopy_species / num_sp * 100
    is_ge_singlecopy = percent_singlecopy_species >= percent_singlecopy_threshold
    return(is_ge_singlecopy)
  }
  num_sp = length(spp_filled)
  is_singlecopy = apply(df_gc[,spp_filled], 1, function(x){is_ge_singlecopy_threshold(x, num_sp, percent_singlecopy_threshold)})
  num_sc = sum(is_singlecopy)
  txt = 'Number of single-copy orthogroups (>=%s percent species) detected for the %s species: %s\n'
  cat(sprintf(txt, percent_singlecopy_threshold, formatC(length(spp_filled), big.mark=','), formatC(num_sc, big.mark=',')))
  return(is_singlecopy)
}

impute_expression = function(dat, num_pc=4) {
  is_all_na_row = apply(dat, 1, function(x){all(is.na(x))})
  tmp = dat[!is_all_na_row,]
  txt = 'Number of removed rows with all NA values in the expression matrix: %s\n'
  cat(sprintf(txt, formatC(sum(is_all_na_row), big.mark=',')))
  num_na = sum(is.na(tmp))
  num_sp = ncol(tmp)
  num_gene = nrow(tmp)
  num_all = num_sp * num_gene
  txt = 'Imputing %s missing values in a total of %s observations (%s genes x %s samples).\n'
  cat(sprintf(txt, formatC(num_na, big.mark=','), formatC(num_all, big.mark=','), formatC(num_gene, big.mark=','), formatC(num_sp, big.mark=',')))
  pc = pcaMethods::pca(tmp, nPcs=num_pc, method='ppca')
  imputed_dat = pcaMethods::completeObs(pc)
  num_negative = sum(imputed_dat < 0)
  txt = 'Number of negative values clipped to zero in the imputed expression matrix: %s\n'
  cat(sprintf(txt, formatC(num_negative, big.mark=',')))
  imputed_dat[imputed_dat < 0] = 0
  return(imputed_dat)
}