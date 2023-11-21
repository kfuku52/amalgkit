
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

write_table_with_index_name = function(df, file_path, index_name='target_id', sort=TRUE) {
    df_index = data.frame(placeholder_name=rownames(df), stringsAsFactors=FALSE)
    colnames(df_index) = index_name
    df = cbind(df_index, df)
    if (sort) {
        df = df[order(df[[index_name]]),]
    }
    write.table(df, file=file_path, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
}

save_exclusion_plot = function(df, out_path, font_size) {
    data_summary = aggregate(cbind(count=exclusion) ~ scientific_name + exclusion, df, length)
    data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN=sum)
    data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
    g = ggplot(data_summary, aes(x=scientific_name, y=count, fill=exclusion))
    g = g + geom_bar(stat = "identity")
    g = g + labs(x = "", y = "Count", fill = "exclusion")
    g = g + theme_bw(base_size=font_size)
    g = g + theme(
        axis.text=element_text(size=font_size, color='black'),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=font_size, color='black'),
        #panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=font_size, color='black'),
        legend.position='bottom',
        rect=element_rect(fill="transparent"),
        plot.margin=unit(rep(0.1, 4), "cm")
    )
    num_spp = length(unique(df[['scientific_name']]))
    plot_width = max(3.6, 0.11 * num_spp)
    ggsave(out_path, plot=g, width=plot_width, height=3.6, units='in')
}