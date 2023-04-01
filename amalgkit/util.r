
get_singlecopy_bool_index = function(df_gc, spp_filled) {
  is_singlecopy = TRUE
  for (sp in spp_filled) {
    is_singlecopy = is_singlecopy & (df_gc[,sp]==1)
  }
  num_sc = sum(is_singlecopy)
  txt = 'Number of single-copy orthogroups detected for the %s species: %s\n'
  cat(sprintf(txt, formatC(length(spp_filled), big.mark=','), formatC(num_sc, big.mark=',')))
  return(is_singlecopy)
}

