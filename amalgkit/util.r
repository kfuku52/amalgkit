
get_singlecopy_bool_index = function(df_gc, spp_filled) {
  is_singlecopy = TRUE
  for (sp in spp_filled) {
    is_singlecopy = is_singlecopy & (df_gc[,sp]==1)
  }
  num_sc = sum(is_singlecopy)
  cat(num_sc, 'single-copy orthogroups were detected for the', length(spp_filled), 'species.\n')
  return(is_singlecopy)
}

