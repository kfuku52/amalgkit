#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(edgeR, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly=TRUE)))

mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')

if (mode=="debug") {
  dir_work = '/Users/s229181/MSN/'
  file_orthogroup_table = file.path(dir_work, "OrthoFinder/Results_Jun22/Orthogroups/Orthogroups.tsv")
  file_genecount = file.path(dir_work, "OrthoFinder/Results_Jun22/Orthogroups/Orthogroups.GeneCount.tsv")
  dir_count = file.path(dir_work, "merge")
  dir_cstmm = file.path(dir_work, "cstmm")
  mode_tmm = 'multi_species'
  setwd(dir_work)
} else if (mode=="batch") {
  args = commandArgs(trailingOnly=TRUE)
  dir_count = args[1]
  file_orthogroup_table = args[2]
  file_genecount = args[3]
  dir_cstmm = args[4]
  mode_tmm = args[5]
  r_util_path = args[6]
}
source(r_util_path)

get_spp_filled = function(dir_count, df_gc=NA) {
  sciname_dirs = list.dirs(dir_count, full.names=FALSE, recursive=FALSE)
  spp_filled = c()
  for (sciname_dir in sciname_dirs) {
    count_files = list.files(path = file.path(dir_count, sciname_dir), pattern = ".*est_counts\\.tsv")
    if (length(count_files)==1) {
      spp_filled = c(spp_filled, count_files)
    } else {
      warning(paste0('Multiple or no est_counts files were detected for ', sciname_dir, ': ', paste(count_files, collapse=', ')))
    }
  }
  spp_filled = sub('_', '|', spp_filled)
  spp_filled = sub('_.*', '', spp_filled)
  spp_filled = sub('\\|', '_', spp_filled)
  if ('data.frame' %in% class(df_gc)) {
    is_missing_in_genecount = (!spp_filled %in% colnames(df_gc))
    if (sum(is_missing_in_genecount)) {
      for (sp in spp_filled[is_missing_in_genecount]) {
        warning(paste0('Species excluded. Not found in the orthogroup table: ', sp))
      }
    }
    spp_filled = spp_filled[!is_missing_in_genecount]
  }
  return(spp_filled)
}

read_est_counts = function(dir_count, sp) {
  sciname_path = file.path(dir_count, sp)
  infile = list.files(path=sciname_path, pattern=".*est_counts\\.tsv")
  if (length(infile)> 1){
    stop(paste0("Multiple *count.tsv files found: ", sp ,"\n"))
  } else if (length(infile)==0) {
    warning(paste0("Skipping. No *est_counts.tsv files found: ", sp ,"\n"))
    return(NULL)
  }
  infile_path = file.path(sciname_path, infile[1])
  cat('Input file found, reading:', infile[1], '\n')
  dat = read.delim(infile_path, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
  dat = dat[,(colnames(dat)!='length'), drop=FALSE]
  colnames(dat) = paste(sp, colnames(dat), sep='_')
  return(dat)
}

get_uncorrected = function(dir_count, file_genecount=NA) {
  if (is.na(file_genecount)) {
    df_gc = NA
  } else {
    df_gc = read.table(file_genecount, header=TRUE, sep='\t', check.names=FALSE, quote='', comment.char='')
  }
  spp_filled = get_spp_filled(dir_count, df_gc)
  uncorrected = list()
  for (sp in spp_filled) {
    dat = read_est_counts(dir_count, sp)
    if (is.null(dat)) {
      next
    }
    uncorrected[[sp]] = dat
  }
  return(uncorrected)
}

get_df_exp_single_copy_ortholog = function(file_genecount, file_orthogroup_table, dir_count, uncorrected) {
  df_gc = read.table(file_genecount, header=TRUE, sep='\t', check.names=FALSE, quote='', comment.char='')
  df_og = read.table(file_orthogroup_table, header=TRUE, sep='\t', row.names=1, check.names=FALSE, quote='', comment.char='')
  spp_filled = get_spp_filled(dir_count, df_gc)
  is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
  df_singleog = df_og[is_singlecopy, spp_filled, drop=FALSE]
  df_sog = df_singleog
  for (sp in spp_filled) {
    if (!sp %in% names(uncorrected)) {
      next
    }
    df_sog = merge(df_sog, uncorrected[[sp]], by.x=sp, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
  }
  df_sog = df_sog[,-(1:length(spp_filled))]
  rownames(df_sog) = rownames(df_singleog)
  return(df_sog)
}

get_df_nonzero = function(df_sog, imputation=TRUE) {
  is_no_count_col = apply(df_sog, 2, function(x){sum(x, na.rm=TRUE)==0})
  txt = 'Removing %s out of %s samples whose read mapping values are all zero.\n'
  cat(sprintf(txt, formatC(sum(is_no_count_col), big.mark=','), formatC(ncol(df_sog), big.mark=',')))
  df_nonzero = df_sog[,!is_no_count_col]
  if (imputation) {
    df_nonzero = impute_expression(df_nonzero)
  } else {
    is_na_containing_row = apply(df_sog, 1, function(x){any(is.na(x))})
    txt = 'Removing %s out of %s orthogroups because missing values are observed in at least one species.\n'
    cat(sprintf(txt, formatC(sum(is_na_containing_row), big.mark=','), formatC(nrow(df_sog), big.mark=',')))
    df_nonzero = df_sog[!is_na_containing_row,]
  }
  return(df_nonzero)
}

create_eff_length_symlink = function(dir_count, dir_cstmm, sp) {
  path_sp = file.path(dir_count, sp)
  eff_length_files = list.files(path = path_sp, pattern = ".*eff_length\\.tsv")
  if (length(eff_length_files)==1) {
    path_target = file.path(path_sp, eff_length_files[1])
    path_link = file.path(dir_cstmm, sp, eff_length_files[1])
    cat('Copying file from', path_target, 'to', path_link, '\n')
    file.copy(from=path_target, to=path_link, overwrite=TRUE)
  } else {
    warning(paste0('No eff_length.tsv file found: ', path_sp))
  }
}

append_tmm_stats_to_metadata = function(df_metadata, cnf_out2) {
  df_nf = cnf_out2[[2]]
  df_nf[['sample']] = rownames(df_nf)
  df_nf[['scientific_name']] = df_nf[['sample']]
  df_nf[['scientific_name']] = sub('_', 'PLACEHOLDER', df_nf[['scientific_name']])
  df_nf[['scientific_name']] = sub('_.*', '', df_nf[['scientific_name']])
  df_nf[['scientific_name']] = sub('PLACEHOLDER', ' ', df_nf[['scientific_name']])
  df_nf[['run']] = sub('.*_', '', df_nf[['sample']])
  df_nf = df_nf[,c('scientific_name', 'run', 'lib.size', 'norm.factors')]
  colnames(df_nf) = c('scientific_name', 'run', 'tmm_library_size', 'tmm_normalization_factor')
  out_cols = c(colnames(df_metadata), colnames(df_nf)[3:ncol(df_nf)])
  df_metadata = merge(df_metadata, df_nf, by=c('scientific_name', 'run'), sort=FALSE, all.x=TRUE, all.y=FALSE)
  df_metadata = df_metadata[,out_cols]
  return(df_metadata)
}

plot_norm_factor_histogram = function(df_metadata, font_size=8) {
  tmp = df_metadata[(!is.na(df_metadata[['tmm_normalization_factor']])),]
  x_limit = max(abs(log2(tmp[['tmm_normalization_factor']])), na.rm=TRUE)
  for (fill_by in c('scientific_name', 'curate_group')) {
    g = ggplot2::ggplot(tmp) +
      geom_histogram(aes(x=log2(tmm_normalization_factor), fill=!!rlang::sym(fill_by)), position="stack", alpha=0.7, bins=40) +
      theme_bw(base_size=font_size) +
      xlim(c(-x_limit, x_limit)) +
      labs(x='log2(TMM normalization factor)', y='Count') +
      guides(fill=guide_legend(ncol=1)) +
      theme(
          axis.text=element_text(size=font_size, color='black'),
          axis.title=element_text(size=font_size, color='black'),
          #panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size=font_size, color='black'),
          rect=element_rect(fill="transparent"),
          plot.margin=unit(rep(0.1, 4), "cm")
      )
    out_path = file.path(dir_cstmm, paste0('tmm_normalization_factor_histogram.', fill_by, '.pdf'))
    ggsave(out_path, plot=g, width=4.8, height=2.4, units='in')
  }
}

plot_norm_factor_scatter = function(df_metadata, font_size=8) {
  tmp = df_metadata[(!is.na(df_metadata[['tmm_normalization_factor']])),]
  x_limit = max(abs(log2(tmp[['tmm_normalization_factor']])), na.rm=TRUE)
  g = ggplot2::ggplot(tmp, aes(x=log10(tmm_library_size), y=log2(tmm_normalization_factor), fill=scientific_name, color=curate_group)) +
    geom_point(shape=21, alpha=0.7) +
    scale_fill_hue(l=65) +
    scale_color_hue(l=45) +
    theme_bw(base_size=font_size) +
    ylim(c(-x_limit, x_limit)) +
    labs(x='log10(Library size)', y='log2(TMM normalization factor)') +
    guides(fill=guide_legend(ncol=1), color=guide_legend(ncol=1)) +
    theme(
        axis.text=element_text(size=font_size, color='black'),
        axis.title=element_text(size=font_size, color='black'),
        #panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.title=element_text(size=font_size, color='black'),
        legend.text=element_text(size=font_size, color='black'),
        rect=element_rect(fill="transparent"),
        plot.margin=unit(rep(0.1, 4), "cm")
    )
  ggsave(file.path(dir_cstmm, 'normalization_factor_scatter.pdf'), plot=g, width=4.8, height=2.0, units='in')
}

if (mode_tmm=='single_species') {
  uncorrected = get_uncorrected(dir_count=dir_count, file_genecount=NA)
  stopifnot(length(names(uncorrected))==1)
  sp = names(uncorrected)[[1]]
  df_sog = uncorrected[[sp]]
  df_nonzero = get_df_nonzero(df_sog)
} else if (mode_tmm=='multi_species') {
  uncorrected = get_uncorrected(dir_count=dir_count, file_genecount=file_genecount)
  df_sog = get_df_exp_single_copy_ortholog(file_genecount, file_orthogroup_table, dir_count, uncorrected)
  df_nonzero = get_df_nonzero(df_sog)
}

cnf_in = edgeR::DGEList(counts = df_nonzero)
cat('Round 1: Performing TMM normalization to determine the appropriate baseline.\n')
cnf_out1 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=NULL)
x = cnf_out1[[2]][['norm.factors']]
cat('Round 1: Median TMM normalization factor =', median(x), '\n')
median_value = sort(x)[ceiling(length(x)/2)]
median_index = (1:length(x))[x==median_value]

cat('Round 2: Performing TMM normalization for output.\n')
cnf_out2 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=median_index)
cat('Round 2: Median TMM normalization factor =', median(cnf_out2[[2]][['norm.factors']]), '\n')

path_metadata = file.path(dir_count, 'metadata.tsv')
df_metadata = read.table(path_metadata, header=TRUE, sep='\t', check.names=FALSE, quote='', comment.char='')
df_metadata = append_tmm_stats_to_metadata(df_metadata, cnf_out2)
out_path = file.path(dir_cstmm, 'metadata.tsv')
write.table(df_metadata, out_path, row.names=FALSE, sep='\t', quote=FALSE)
plot_norm_factor_histogram(df_metadata=df_metadata)
plot_norm_factor_scatter(df_metadata=df_metadata)

for (sp in names(uncorrected)) {
  cat('Applying TMM normalization factors:', sp, '\n')
  dat = uncorrected[[sp]]
  df_nf_sp = cnf_out2[[2]][startsWith(rownames(cnf_out2[[2]]),sp),]
  for (i in 1:length(df_nf_sp[,1])){
    SRR = as.character(row.names(df_nf_sp[i,]))
    tmm_normalization_factor = as.double(df_nf_sp[i,"norm.factors"]) # manually apply normfactor
    dat[,SRR] = dat[,SRR]/tmm_normalization_factor
  }
  dat_out = cbind(target_id=rownames(dat), dat)
  rownames(dat_out) = NULL
  colnames(dat_out) = sub(paste0(sp, '_'), '', colnames(dat_out))
  dir_cstmm_sp = file.path(dir_cstmm, sp)
  if (!file.exists(dir_cstmm_sp)) {
    dir.create(dir_cstmm_sp)
  }
  create_eff_length_symlink(dir_count, dir_cstmm, sp)
  file_path = file.path(dir_cstmm_sp, paste0(sp, "_cstmm_counts.tsv"))
  write.table(dat_out, file_path, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
}
cat('cstmm.r completed!\n')