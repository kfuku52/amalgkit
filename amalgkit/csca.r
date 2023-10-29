#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(amap, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(colorspace, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(dendextend, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(NMF, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(MASS, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pvclust, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(Rtsne, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pcaMethods, quietly=TRUE)))
options(stringsAsFactors = FALSE)

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
font_size = 8

if (debug_mode == "debug") {
  selected_curate_groups = c('root','flower', 'leaf')
  dir_work = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_gfe_pipeline/amalgkit_out'
  dir_csca_input_table = file.path(dir_work,'csca/csca_input_symlinks')
  file_orthogroup = file.path(dir_work, 'csca/multispecies_busco_table.tsv')
  file_genecount = file.path(dir_work, 'csca/multispecies_genecount.tsv')
  r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
  dir_csca = file.path(dir_work, 'csca')
  batch_effect_alg = 'sva'
} else if (debug_mode == "batch") {
  args = commandArgs(trailingOnly = TRUE)
  selected_curate_groups = strsplit(args[1], "\\|")[[1]]
  dir_work = args[2]
  dir_csca_input_table = args[3]
  file_orthogroup = args[4]
  file_genecount = args[5]
  r_util_path = args[6]
  dir_csca = args[7]
  batch_effect_alg = args[8]
}
source(r_util_path)
setwd(dir_csca)
cat('selected_curate_groups:', selected_curate_groups, "\n")
cat('dir_work:', dir_work, "\n")
cat('dir_csca_input_table:', dir_csca_input_table, "\n")
cat('file_orthogroup:', file_orthogroup, "\n")
cat('file_genecount:', file_genecount, "\n")
cat('r_util_path:', r_util_path, "\n")
cat('dir_csca:', dir_csca, "\n")
cat('batch_effect_alg:', batch_effect_alg, "\n")

add_color_to_metadata = function(df, selected_curate_groups) {
  df = df[,(!colnames(df) %in% c('bp_color','sp_color','curate_group_color'))]
  scientific_name = as.character(df[['scientific_name']])
  curate_group = as.character(df[['curate_group']])
  scientific_name_unique = sort(scientific_name[!duplicated(scientific_name)])
  curate_group_unique = sort(curate_group[!duplicated(curate_group)])
  if (length(selected_curate_groups) <= 8) {
    curate_group_color = brewer.pal(length(unique(curate_group)), "Dark2")
    sp_color = rainbow_hcl(length(unique(scientific_name)), c=100)
  } else if (length(selected_curate_groups) <= 12) {
    curate_group_color = brewer.pal(length(unique(curate_group)), "Paired")
    sp_color = rainbow_hcl(length(unique(scientific_name)), c=100)
  } else {
    curate_group_color = rainbow_hcl(length(selected_curate_groups), c=100)
    sp_color = rainbow_hcl(length(unique(scientific_name)), c=150)
  }
  df_curate_group = data.frame(curate_group=sort(curate_group_unique), curate_group_color=curate_group_color[1:length(curate_group_unique)], stringsAsFactors=FALSE)
  df_sp = data.frame(scientific_name=scientific_name_unique, sp_color=sp_color[1:length(scientific_name_unique)], stringsAsFactors=FALSE)
  df = merge(df, df_sp, sort=FALSE, all.y=FALSE)
  df = merge(df, df_curate_group, sort=FALSE, all.y=FALSE)
  if ('bioproject' %in% colnames(df)) {
    bioproject = as.character(df$bioproject)
    if (length(selected_curate_groups) <= 8) {
      bp_color = rainbow_hcl(length(unique(bioproject)), c=50)
    } else if (length(selected_curate_groups) <= 12) {
      bp_color = rainbow_hcl(length(unique(bioproject)), c=50)
    } else {
      bp_color = rainbow_hcl(length(unique(bioproject)), c=50)
    }
    df_bp = data.frame(bioproject=unique(bioproject), bp_color=bp_color[1:length(unique(bioproject))], stringsAsFactors=FALSE)
    df = merge(df, df_bp, sort=FALSE, all.y=FALSE)
  }
  return(df)
}

sort_labels = function(df_label, label_orders) {
  df_tmp = data.frame()
  for (lo in label_orders) {
    splits = strsplit(lo, '_')[[1]]
    scientific_name = paste(splits[1], splits[2])
    curate_group = splits[3]
    df_tmp = rbind(df_tmp, df_label[(df_label[['scientific_name']]==scientific_name)&(df_label[['curate_group']]==curate_group),])
  }
  return(df_tmp)
}

sort_averaged_tc = function(tc) {
  split_colnames = strsplit(colnames(tc), "_")
  genus_names = c()
  specific_names = c()
  curate_group_names = c()
  for (i in 1:length(split_colnames)) {
    genus_names = c(genus_names, split_colnames[[i]][1])
    specific_names = c(specific_names, split_colnames[[i]][2])
    curate_group_names = c(curate_group_names, split_colnames[[i]][3])
  }
  colname_order = order(curate_group_names, genus_names, specific_names)
  tc = tc[, colname_order]
  return(tc)
}

color_children2parent = function(node) {
  if (length(node)!=2) {
    return(node)
  }
  if (!is.null(attributes(node[[1]])) && !is.null(attributes(node[[1]])$edgePar)) {
    child1_color = attributes(node[[1]])$edgePar[['col']]
  } else {
    child1_color = NULL
  }
  if (!is.null(attributes(node[[2]])) && !is.null(attributes(node[[2]])$edgePar)) {
    child2_color = attributes(node[[2]])$edgePar[['col']]
  } else {
    child2_color = NULL
  }
  if (is.null(child1_color) | is.null(child2_color)) {
    return(node)
  }
  if (is.na(child1_color) | is.na(child2_color)) {
    return(node)
  }
  if (child1_color==child2_color) {
    attributes(node)$edgePar[['col']] = child1_color
  }
  return(node)
}

map_color = function(redundant_variables, c) {
  uniq_var = unique(redundant_variables)
  uniq_col = rainbow_hcl(length(uniq_var), c=c)
  df_unique = data.frame(var=uniq_var, col=uniq_col, stringsAsFactors=FALSE)
  df_redundant = data.frame(var=redundant_variables, order=seq(1, length(redundant_variables)), stringsAsFactors=FALSE)
  df_redundant = merge(df_redundant, df_unique, by="var", all.x=TRUE, stringsAsFactors=FALSE)
  df_redundant = df_redundant[order(df_redundant[['order']]),]
  return(df_redundant[['col']])
}

draw_multisp_heatmap = function(tc, df_label) {
  tc_dist_matrix = cor(tc, method='pearson')
  tc_dist_matrix[is.na(tc_dist_matrix)] = 0
  ann_label = df_label[,c('scientific_name','curate_group')]
  colnames(ann_label) = c('species', 'curate_group')
  sp_color = df_label[!duplicated(df_label[['sp_color']]), 'sp_color']
  sp_color = sp_color[order(df_label[!duplicated(df_label[['scientific_name']]), 'scientific_name'])]
  curate_group_color = df_label[!duplicated(df_label[['curate_group_color']]), 'curate_group_color']
  curate_group_color = curate_group_color[order(df_label[!duplicated(df_label[['curate_group']]), 'curate_group'])]
  ann_color = list(species=sp_color, curate_group=curate_group_color)
  breaks = c(0, seq(0.3, 1, 0.01))
  NMF::aheatmap(tc_dist_matrix, color="-RdYlBu2:71", Rowv=NA, Colv=NA, revC=TRUE, legend=TRUE, breaks=breaks,
           annCol=ann_label, annRow=ann_label, annColors=ann_color, annLegend=FALSE, labRow=NA, labCol=NA)
}

draw_multisp_dendrogram = function(tc, df_label, df_metadata, nboot, cex.xlab, cex.yaxis, pvclust_file='pvclust.RData') {
  colnames(tc) = sub("_.*","",sub('_',' ',colnames(tc)))
  dist_fun = function(x){Dist(t(x), method='pearson')}
  if (file.exists(pvclust_file)) {
    if (file.info(pvclust_file)$size) {
      cat('pvclust file found.\n')
      load(pvclust_file)
    }
  } else {
    cat('No pvclust file found. Start bootstrapping.\n')
    result = pvclust(tc, method.dist=dist_fun, method.hclust="average", nboot=nboot, parallel=FALSE) # UPGMA
    save(result, file=pvclust_file)
  }
  dend = as.dendrogram(result)
  dend_colors = df_label[order.dendrogram(dend), 'curate_group_color']
  label_colors = df_label[order.dendrogram(dend), 'sp_color']
  labels_colors(dend) = label_colors
  dend_labels <- df_metadata[order.dendrogram(dend), 'run']
  dend <- color_branches(dend, labels=dend_labels, col=dend_colors)
  dend <- set(dend, "branches_lwd", 2)
  for (i in 1:ncol(tc)) {
    dend = dendrapply(dend, color_children2parent)
  }
  par(cex=cex.xlab)
  plot(dend, las=1, yaxt='n')
  text(result, print.num=FALSE, cex=1, col.pv='black')
  par(cex=cex.yaxis)
  axis(2, las=1)
  mtext(text='Distance', side=2, line=4, cex=cex.yaxis)

  n = ncol(tc)
  f = 100
  curate_group_unique = unique(df_metadata['curate_group'])
  sp_unique = unique(df_metadata[['scientific_name']])
  bp_unique = unique(df_metadata[['bioproject']])
  curate_group_color_unique = unique(df_metadata[['curate_group_color']])
  sp_color_unique = unique(df_metadata[['sp_color']])
  bp_color_unique = unique(df_metadata[['bp_color']])
  legend_text = c(as.character(curate_group_unique), "", as.character(sp_unique), "", as.character(bp_unique))
  legend_bg = c(curate_group_color_unique, "white", sp_color_unique, "white", bp_color_unique)
  legend_fg = c(rep("black", length(curate_group_color_unique)), "white", rep("black", length(sp_color_unique)), "white", rep("black", length(bp_color_unique)))
  #plot.new() ; par(mar=c(0,0,0,0))
  #legend("center", legend=legend_text, cex=1, pch=22, lty=0, lwd=1, pt.bg=legend_bg, col=legend_fg)
}

draw_multisp_pca = function(tc, df_label) {
  tc_dist_matrix = cor(tc, method='pearson')
  tc_dist_matrix[is.na(tc_dist_matrix)] = 0
  set.seed(1)
  pca = prcomp(tc_dist_matrix)
  xlabel = paste0("PC 1 (", round(summary(pca)$importance[2,1]*100, digits=1), "%)")
  ylabel = paste0("PC 2 (", round(summary(pca)$importance[2,2]*100, digits=1), "%)")
  plot(pca$rotation[,1], pca$rotation[,2], pch=21, cex=2, lwd=1, bg=df_label$curate_group_color, col=df_label$sp_color, xlab=xlabel, ylab=ylabel, las=1)
}

draw_multisp_mds = function(tc, df_label) {
  tc_dist_dist = Dist(t(tc), method='pearson') + 0.000000001
  tc_dist_dist[is.na(tc_dist_dist)] = 1
  set.seed(1)
  try_out = tryCatch(
    {isoMDS(tc_dist_dist, k=2, maxit=100)},
    error = function(a){return("MDS failed.")}
  )
  if (mode(try_out)=="character") {
    cat('MDS failed.\n')
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  } else {
    mds <- try_out
    plot(mds$points[,1], mds$points[,2], pch=21, cex=2, lwd=1, bg=df_label$curate_group_color, col=df_label$sp_color, xlab="MDS dimension 1", ylab="MDS dimension 2", las=1)
  }
}

draw_multisp_tsne = function(tc, df_label) {
  perplexity = min(30, floor(ncol(tc)/4), na.rm=TRUE)
  set.seed(1)
  out_tsne = Rtsne(as.matrix(t(tc)), theta=0, check_duplicates=FALSE, verbose=FALSE, perplexity=perplexity, dims=2)
  try_out = tryCatch(
    {
      plot(out_tsne$Y[,1], out_tsne$Y[,2], pch=21, cex=2, lwd=1, bg=df_label$curate_group_color, col=df_label$sp_color,
           xlab="t-SNE dimension 1", ylab="t-SNE dimension 2", las=1)
    },
    error = function(a){return("t-SNE plot failed.")}
  )
  if (mode(try_out)=="character") {
    cat('t-SNE failed.\n')
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  }
}

draw_multisp_legend = function(df_label) {
  cex_axis=0.7
  curate_group_unique = df_label$curate_group[!duplicated(df_label$curate_group)]
  sp_unique = df_label$scientific_name[!duplicated(df_label$scientific_name)]
  curate_group_color_unique = df_label$curate_group_color[!duplicated(df_label$curate_group_color)]
  sp_color_unique = df_label$sp_color[!duplicated(df_label$sp_color)]
  toumei=rgb(1,1,1,0)
  legend_text = c('Tissue', as.character(curate_group_unique), "", 'Species', as.character(sp_unique))
  legend_bg = c(toumei, curate_group_color_unique, toumei, toumei, rep(toumei, length(sp_color_unique)))
  legend_fg = c(toumei, rep(toumei, length(curate_group_color_unique)), toumei, toumei, sp_color_unique)
  legend_pch = c(1, rep(21,length(curate_group_color_unique)), 1, 1, rep(1,length(sp_color_unique)))
  legend_font = c(2, rep(1, length(curate_group_color_unique)), 1, 2, rep(3, length(sp_color_unique)))
  plot.new()
  legend("right", legend=legend_text, pt.cex=1, pch=legend_pch, lty=0, lwd=2, pt.bg=legend_bg, col=legend_fg, cex=cex_axis, text.font=legend_font)
}

prepare_metadata_table = function(dir_csca_input_table, selected_curate_groups, spp) {
  files = list.files(dir_csca_input_table, pattern = ".*metadata.*")
  df_metadata = data.frame()
  for (file in files) {
    metadata_path = file.path(dir_csca_input_table, file)
    tmp_metadata = read.table(metadata_path, header=TRUE, sep='\t', quote='', comment.char='', check.names=FALSE)
    df_metadata = rbind(df_metadata, tmp_metadata)
  }
  df_metadata = df_metadata[(df_metadata[['curate_group']] %in% selected_curate_groups)&(df_metadata[['scientific_name']] %in% spp),]
  df_metadata = df_metadata[,!startsWith(colnames(df_metadata), 'Unnamed')]
  return(df_metadata)
}

get_label_orders = function(df_metadata) {
  order_cg = order(df_metadata[['curate_group']])
  label_orders = unique(paste(df_metadata[order_cg,'scientific_name'], df_metadata[order_cg,'curate_group'], sep='_'))
  label_orders = sub(' ', '_', label_orders)
  return(label_orders)
}

extract_ortholog_mean_expression_table = function(df_singleog, averaged_tcs, label_orders) {
  averaged_orthologs = list()
  averaged_orthologs[['uncorrected']] = df_singleog
  averaged_orthologs[['corrected']] = df_singleog
  row_names = rownames(averaged_orthologs[['corrected']])
  for (d in c('uncorrected', 'corrected')) {
    for (sp_filled in colnames(df_singleog)) {
      tc = averaged_tcs[[d]][[sp_filled]]
      averaged_orthologs[[d]] = merge(averaged_orthologs[[d]], tc, by.x=sp_filled, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
    }
    num_remove_col = ncol(df_singleog)
    averaged_orthologs[[d]] = averaged_orthologs[[d]][,-(1:num_remove_col)]
    rownames(averaged_orthologs[[d]]) = row_names
    averaged_orthologs[[d]] = sort_averaged_tc(averaged_orthologs[[d]])
    available_label_orders = label_orders[label_orders %in% colnames(averaged_orthologs[[d]])]
    averaged_orthologs[[d]] = averaged_orthologs[[d]][, available_label_orders]
  }
  cat(nrow(averaged_orthologs[[d]]), 'orthologs were found before filtering.\n')
  return(averaged_orthologs)
}

load_unaveraged_expression_tables = function(dir_csca_input_table, spp_filled, batch_effect_alg) {
  unaveraged_tcs = list()
  unaveraged_tcs[['uncorrected']] = list()
  unaveraged_tcs[['corrected']] = list()
  all_files = list.files(dir_csca_input_table, pattern="*.tc.tsv")
  uncorrected_files = all_files[grepl("uncorrected", all_files)]
  corrected_files = all_files[((!grepl("uncorrected", all_files))&(grepl(batch_effect_alg, all_files)))]
  for (sp in spp_filled) {
    uncorrected_file = uncorrected_files[startsWith(uncorrected_files, sp)]
    uncorrected_path = file.path(dir_csca_input_table, uncorrected_file)
    corrected_file = corrected_files[startsWith(corrected_files, sp)]
    corrected_path = file.path(dir_csca_input_table, corrected_file)
    if ((length(uncorrected_path)==0)|(length(corrected_path)==0)) {
      cat(paste("Skipping. `amalgkit curate` output(s) not found:", sp, "\n"), file=stderr())
      next
    }
    unaveraged_tcs[['uncorrected']][[sp]] = read.delim(uncorrected_path, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
    unaveraged_tcs[['corrected']][[sp]] = read.delim(corrected_path, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
  }
  return(unaveraged_tcs)
}

extract_ortholog_unaveraged_expression_table = function(df_singleog, unaveraged_tcs) {
  unaveraged_orthologs = list()
  unaveraged_orthologs[['uncorrected']] = df_singleog
  unaveraged_orthologs[['corrected']] = df_singleog
  row_names = rownames(unaveraged_orthologs[['corrected']])
  for (d in c('uncorrected', 'corrected')) {
    for (sp_filled in colnames(df_singleog)) {
      tc = unaveraged_tcs[[d]][[sp_filled]]
      colnames(tc) = paste(sp_filled, colnames(tc), sep='_')
      unaveraged_orthologs[[d]] = merge(unaveraged_orthologs[[d]], tc, by.x=sp_filled, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
    }
    num_remove_col = length(spp)
    unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][,-(1:num_remove_col)]
    rownames(unaveraged_orthologs[[d]]) = row_names
    tc_order = order(sub('.*_','',colnames(unaveraged_orthologs[[d]])))
    unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][,tc_order]
  }
  return(unaveraged_orthologs)
}

get_df_labels_averaged = function(df_metadata, label_orders) {
    metadata_tmp = df_metadata[(df_metadata[['exclusion']]=='no'),]
    df_label = unique(metadata_tmp[,c('scientific_name','curate_group')])
    categories = list(scientific_name=metadata_tmp[['scientific_name']], curate_group=metadata_tmp[['curate_group']])
    df_bp = aggregate(metadata_tmp[['bioproject']], by=categories, function(x){length(unique(x))})
    colnames(df_bp) = c('scientific_name','curate_group','num_bp')
    df_label = merge(df_label, df_bp, all.x=TRUE, all.y=FALSE)
    df_run = aggregate(metadata_tmp[['run']], by=categories, function(x){length(unique(x))})
    colnames(df_run) = c('scientific_name','curate_group','num_run')
    df_label = merge(df_label, df_run, all.x=TRUE, all.y=FALSE)
    df_label = df_label[order(df_label[['curate_group']], df_label[['scientific_name']]),]
    df_label = sort_labels(df_label, label_orders)
    df_label = add_color_to_metadata(df_label, selected_curate_groups)
    df_label = sort_labels(df_label, label_orders)
    rownames(df_label) = NULL
    write.table(df_label, paste0('csca_color_averaged.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
    return(df_label)
}

get_df_labels_unaveraged = function(df_metadata, selected_curate_groups) {
    cols = c('run','bioproject','curate_group','scientific_name','sp_color','curate_group_color','bp_color')
    metadata_tmp = df_metadata[(df_metadata[['exclusion']]=='no'),]
    df_color = add_color_to_metadata(metadata_tmp, selected_curate_groups)
    df_color = df_color[,cols]
    label_order = order(df_color[['run']])
    df_color = df_color[label_order,]
    write.table(df_color, paste0('csca_color_unaveraged.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
    return(df_color)
}

save_averaged_tsne_plot = function(tc, df_label) {
  cat('Generating averaged t-SNE plot.\n')
  perplexity = min(30, floor(ncol(tc)/4), na.rm=TRUE)
  set.seed(1)
  out_tsne = try(Rtsne(as.matrix(t(tc)), theta=0, check_duplicates=FALSE, verbose=FALSE, perplexity=perplexity, dims=2))
  if ("try-error" %in% class(out_tsne)) {
    flag_tsne_success = FALSE
    print(out_tsne)
  } else {
    flag_tsne_success = TRUE
  }
  if (!flag_tsne_success) {
    return()
  }
  try_out = tryCatch(
    {
      plot(out_tsne$Y[,1], out_tsne$Y[,2], pch=21, cex=2, lwd=1, bg=df_label[['curate_group_color']],
           col=df_label[['sp_color']], xlab="t-SNE dimension 1", ylab="t-SNE dimension 2", las=1)
    },
    error = function(a){return("t-SNE plot failed.")}
  )
  if (mode(try_out)=="character") {
    cat('t-SNE failed.\n')
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  }
}

get_pca_coordinates = function(tc, df_label, by='species_curate_group') {
  tc_dist_matrix = cor(tc, method='pearson')
  tc_dist_matrix[is.na(tc_dist_matrix)] = 0
  #set.seed(1)
  pca = prcomp(tc_dist_matrix)
  labels = c()
  for (i in 1:5) {
    labels = c(labels, paste0("Principal component ", i, " (", round(summary(pca)$importance[2,i]*100, digits=1), "%)"))
  }
  PC1 = pca[['x']][,'PC1']
  PC2 = pca[['x']][,'PC2']
  PC3 = pca[['x']][,'PC3']
  PC4 = pca[['x']][,'PC4']
  PC5 = pca[['x']][,'PC5']
  tmp = data.frame(PC1, PC2, PC3, PC4, PC5)
  if (by=='species_curate_group') {
    df_label[by] = paste0(sub(' ', '_', df_label[['scientific_name']]), '_', df_label[['curate_group']])
    tmp[by] = rownames(tmp)
  } else if (by=='run') {
    tmp[by] = sub('.*_','',rownames(tmp))
  } else {
    tmp[by] = rownames(tmp)
  }
  tmp = merge(df_label, tmp, by=by)
  return(list(tmp, labels))
}

save_unaveraged_pca_plot = function(unaveraged_orthologs, df_color_unaveraged, df_metadata) {
  cat('Generating unaveraged PCA plot.\n')
  for (d in c('uncorrected', 'corrected')) {
    out = get_pca_coordinates(tc=unaveraged_orthologs[[d]], df_label=df_color_unaveraged, by='run')
    tmp = out[[1]]
    pc_contributions = out[[2]]
    pc_cols = c('PC1','PC2','PC3','PC4','PC5')
    pc_cols2 = paste(pc_cols, d, sep='_')
    sorted_cols = c(colnames(df_metadata), pc_cols2)
    tmp2 = tmp[,c('run',pc_cols)]
    colnames(tmp2) = c('run', pc_cols2)
    df_metadata = merge(df_metadata, tmp2, all.x=TRUE, by='run', sort=FALSE)
    df_metadata = df_metadata[,sorted_cols]
    for (pcxy in list(c(1,2),c(3,4))) {
      pcx = pcxy[1]
      pcy = pcxy[2]

      colx = paste0('PC', pcx)
      coly = paste0('PC', pcy)
      xmin = min(tmp[[colx]], na.rm=TRUE)
      xmax = max(tmp[[colx]], na.rm=TRUE)
      xunit = (xmax-xmin)*0.01
      xmin = xmin - xunit
      xmax = xmax + xunit

      ymin = min(tmp[[coly]], na.rm=TRUE)
      ymax = max(tmp[[coly]], na.rm=TRUE)
      yunit = (ymax-ymin)*0.01
      ymin = ymin - yunit
      ymax = ymax + yunit

      curate_group_colors = unique(df_color_unaveraged[,c('curate_group','curate_group_color')])[['curate_group_color']]

      g = ggplot(tmp, aes(x=!!rlang::sym(colx), y=!!rlang::sym(coly), color=curate_group))
      g = g + theme_bw()
      g = g + geom_point(size=0.5, alpha=0.3)
      g = g + geom_density_2d(mapping=aes(color=curate_group), bins=12, linewidth=0.25)
      g = g + geom_density_2d(mapping=aes(color=curate_group), bins=12, linewidth=0.25)
      g = g + scale_color_manual(values=curate_group_colors)
      g = g + xlab(pc_contributions[pcx])
      g = g + ylab(pc_contributions[pcy])
      g = g + xlim(xmin, xmax)
      g = g + ylim(ymin, ymax)
      g = g + theme(
        axis.text=element_text(size=font_size),
        axis.title=element_text(size=font_size),
        legend.text=element_text(size=font_size),
        legend.title=element_text(size=font_size)
      )
      filename = paste0('csca_unaveraged_pca_PC', pcx, pcy, '.', d, '.pdf')
      ggsave(file=filename, g, height=2.15, width=4.25)
    }
  }
  return(df_metadata)
}

get_tsne_coordinates = function(tc, df_label, by='run') {
  perplexity = min(30, floor(ncol(tc)/4), na.rm=TRUE)
  set.seed(1)
  out_tsne = Rtsne(as.matrix(t(tc)), theta=0, check_duplicates=FALSE, verbose=FALSE, perplexity=perplexity, dims=2)
  tmp = data.frame(tsne1=out_tsne$Y[,1], tsne2=out_tsne$Y[,2])
  tmp[[by]] = sub('.*_', '', colnames(tc))
  tmp = merge(df_label, tmp, by=by)
  return(tmp)
}

save_unaveraged_tsne_plot = function(unaveraged_orthologs, df_color_unaveraged) {
  cat('Generating unaveraged t-SNE plot.\n')
  for (d in c('uncorrected', 'corrected')) {
    tmp = get_tsne_coordinates(tc=unaveraged_orthologs[[d]], df_label=df_color_unaveraged)
    pcx = 1
    pcy = 2
    colx = paste0('tsne', pcx)
    coly = paste0('tsne', pcy)
    xmin = min(tmp[[colx]], na.rm=TRUE)
    xmax = max(tmp[[colx]], na.rm=TRUE)
    xunit = (xmax-xmin)*0.01
    xmin = xmin - xunit
    xmax = xmax + xunit

    ymin = min(tmp[[coly]], na.rm=TRUE)
    ymax = max(tmp[[coly]], na.rm=TRUE)
    yunit = (ymax-ymin)*0.01
    ymin = ymin - yunit
    ymax = ymax + yunit

    g = ggplot(tmp, aes(x=!!rlang::sym(colx), !!rlang::sym(coly), color=curate_group))
    g = g + theme_bw()
    g = g + geom_point(size=0.5)
    g = g + geom_density_2d(mapping=aes(color=curate_group), bins=12, linewidth=0.25)
    curate_group_colors = unique(df_color_unaveraged[,c('curate_group','curate_group_color')])[['curate_group_color']]
    g = g + scale_color_manual(values=curate_group_colors)
    g = g + xlab('t-SNE dimension 1')
    g = g + ylab('t-SNE dimension 2')
    g = g + xlim(xmin, xmax)
    g = g + ylim(ymin, ymax)
    g = g + theme(
      axis.text=element_text(size=font_size),
      axis.title=element_text(size=font_size),
      legend.text=element_text(size=font_size),
      legend.title=element_text(size=font_size)
    )
    filename = paste0('csca_unaveraged_tsne.', d, '.pdf')
    ggsave(file=filename, g, height=2.15, width=4.25)
  }
}

save_averaged_heatmap_plot = function(averaged_orthologs, df_color_averaged) {
  cat('Generating averaged heatmap.\n')
  file_name='csca_SVA_heatmap.pdf'
  pdf(file_name, height=3.3, width=7.2) # full figure size = 9.7 x 7.2
  layout_matrix=matrix(c(
    1,2,2,2,2,2,2,2,2,2,
    1,3,3,3,3,3,3,3,3,3),
    2,10,byrow=TRUE)
  layout(t(layout_matrix))
  par(mar=c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(0.27,0.5,'Uncorrected', srt=0, cex=2)
  text(0.80,0.5,'Corrected', srt=0, cex=2)
  for (d in c('uncorrected','corrected')) {
    tc = averaged_orthologs[[d]]
    df_label = df_color_averaged
    par(mar=c(0,0,0,0))
    draw_multisp_heatmap(tc=tc, df_label=df_label)
  }
  graphics.off()
}

save_averaged_dendrogram_plot = function(averaged_orthologs, df_color_averaged) {
  cat('Generating averaged dendrogram.\n')
  file_name='csca_SVA_dendrogram.pdf'
  pdf(file_name, height=2.5, width=7.2) # full figure size = 9.7 x 7.2
  layout_matrix=matrix(c(1,2),2,1,byrow=TRUE)
  layout(t(layout_matrix))
  for (d in c('uncorrected','corrected')) {
    tc = averaged_orthologs[[d]]
    df_label = df_color_averaged
    par(cex=0.5, mar=c(10,5.5,0,0), mgp=c(4, 0.7, 0))
    pvclust_file=paste0('csca_pvclust_',d,'.RData')
    draw_multisp_dendrogram(tc=tc, df_label=df_label, df_metadata=df_metadata, pvclust_file=pvclust_file,
                            nboot=1, cex.xlab=0.3, cex.yaxis=0.5)
  }
  graphics.off()
}

save_averaged_dimensionality_reduction_summary = function(averaged_orthologs, df_color_averaged) {
  cat('Generating averaged dimensionality reduction summary.\n')
  par(cex=1)
  file_name='csca_averaged_summary.pdf'
  pdf(file_name, height=7.2, width=7.2) # full figure size = 9.7 x 7.2
  layout_matrix = matrix(c(1,1,1,4,4,4,7,7,2,2,2,5,5,5,7,7,3,3,3,6,6,6,7,7),3,8,byrow=TRUE)
  layout(layout_matrix)
  for (d in c('uncorrected','corrected')) {
    tc = averaged_orthologs[[d]]
    df_label = df_color_averaged
    par(mar=c(4,4,0.1,1)); draw_multisp_pca(tc=tc, df_label=df_label)
    par(mar=c(4,4,0.1,1)); draw_multisp_tsne(tc=tc, df_label=df_label)
    par(mar=c(4,4,0.1,1)); draw_multisp_mds(tc=tc, df_label=df_label)
  }
  par(mar=c(0,0,0,0)); draw_multisp_legend(df_label)
  graphics.off()
}

draw_multisp_boxplot = function(df_metadata, tc_dist_matrix, fontsize=8) {
  is_same_sp = outer(df_metadata[['scientific_name']], df_metadata[['scientific_name']], function(x,y){x==y})
  is_same_curate_group = outer(df_metadata[['curate_group']], df_metadata[['curate_group']], function(x,y){x==y})
  plot(c(0.5, 4.5), c(0, 1), type = 'n', xlab='', ylab="Pearson's correlation\ncoefficient", las=1, xaxt='n')
  boxplot(tc_dist_matrix[(!is_same_sp)&(!is_same_curate_group)], at=1, add=TRUE, col='gray', yaxt='n')
  boxplot(tc_dist_matrix[(is_same_sp)&(!is_same_curate_group)], at=2, add=TRUE, col='gray', yaxt='n')
  boxplot(tc_dist_matrix[(!is_same_sp)&(is_same_curate_group)], at=3, add=TRUE, col='gray', yaxt='n')
  boxplot(tc_dist_matrix[(is_same_sp)&(is_same_curate_group)], at=4, add=TRUE, col='gray', yaxt='n')
  labels = c('bw\nbw', 'bw\nwi', 'wi\nbw', 'wi\nwi')
  axis(side=1, at=c(1,2,3,4), labels=labels, padj=0.5)
  axis(side=1, at=0.35, labels='Group\nSpecies', padj=0.5, hadj=1, tick=FALSE)
}

save_averaged_box_plot = function(averaged_orthologs, df_color_averaged) {
  cat('Generating averaged boxplot.\n')
  file_name='csca_boxplot.pdf'
  pdf(file_name, height=3.6, width=7.2) # full figure size = 9.7 x 7.2
  par(mfrow=c(1,2))
  for (d in c('uncorrected','corrected')) {
    tc = averaged_orthologs[[d]]
    tc[tc<0] = 0
    tc_dist_matrix = cor(tc, method='pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    draw_multisp_boxplot(df_color_averaged, tc_dist_matrix, fontsize=8)
  }
  d = 'corrected'
  tc = averaged_orthologs[[d]]
  df_label = df_color_averaged
  par(mar=c(2.5,2.5,0.3,0.1), cex=1, ps=8, mgp=c(1.5, 0.7, 0)); draw_multisp_tsne(tc=tc, df_label=df_label[['sp_color']])
  graphics.off()
}

calculate_correlation_within_group = function(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_curate_groups, dist_method='pearson') {
    for (d in c('uncorrected', 'corrected')) {
        ortholog_med = data.frame(matrix(NA, nrow(averaged_orthologs[[d]]), length(selected_curate_groups)))
        colnames(ortholog_med) = selected_curate_groups
        rownames(ortholog_med) = rownames(averaged_orthologs[[d]])
        for (curate_group in selected_curate_groups) {
            is_curate_group = endsWith(colnames(averaged_orthologs[[d]]), curate_group)
            ortholog_med[,curate_group] = apply(averaged_orthologs[[d]][,is_curate_group], 1, function(x){median(x, na.rm=TRUE)})
        }
        stopifnot(all(rownames(unaveraged_orthologs[[d]])==rownames(ortholog_med)))
        target_col = paste0('within_group_cor_', d)
        nongroup_col = paste0('max_nongroup_cor_', d)
        df_metadata[,target_col] = NA
        df_metadata[,nongroup_col] = NA
        for (sp_and_run in colnames(unaveraged_orthologs[[d]])) {
            split_string = strsplit(sp_and_run, "_")[[1]]
            sp = paste(split_string[1:2], collapse=" ")
            sra_run = paste(split_string[3:length(split_string)], collapse='_')
            is_sra = (df_metadata[['run']]==sra_run) & (df_metadata[['scientific_name']]==sp)
            sample_cg = df_metadata[is_sra,'curate_group']
            sample_values = unaveraged_orthologs[[d]][,sp_and_run]
            for (curate_group in colnames(ortholog_med)) {
                med_values = ortholog_med[,curate_group]
                is_na = (is.na(sample_values) | is.na(med_values))
                sample_values2 = sample_values[!is_na]
                med_values2 = med_values[!is_na]
                cor_coef = cor(sample_values2, med_values2, method=dist_method)
                if (sample_cg==curate_group) {
                    df_metadata[is_sra,target_col] = cor_coef
                } else {
                    df_metadata[is_sra,nongroup_col] = max(cor_coef, df_metadata[is_sra,nongroup_col], na.rm=TRUE)
                }
            }
        }
    }
    return(df_metadata)
}

save_group_cor_histogram = function(df_metadata, font_size=8) {
  cat('Generating unaveraged group correlation histogram.\n')
  max_count <- 0
  for (col in c('within_group_cor_uncorrected', 'within_group_cor_corrected')) {
    for (fill_by in c('curate_group', 'scientific_name')) {
      tmp = df_metadata[(!is.na(df_metadata[[col]])),]
      bin_counts <- table(cut(tmp[[col]], breaks = seq(0, 1, length.out = 41)))
      max_count <- max(max_count, max(bin_counts, na.rm=TRUE), na.rm=TRUE)  # Update max_count if necessary
    }
  }
  plot_list <- list()
  for (col in c('within_group_cor_uncorrected', 'within_group_cor_corrected')) {
    for (fill_by in c('curate_group', 'scientific_name')) {
      tmp = df_metadata[(!is.na(df_metadata[[col]])),]
      g = ggplot2::ggplot(tmp) +
        geom_histogram(aes(x=!!rlang::sym(col), fill=!!rlang::sym(fill_by)),
                       position="stack", alpha=0.7, bins=40) +
        theme_bw(base_size=font_size) +
        xlim(c(0, 1)) +
        ylim(c(0, max_count)) +
        labs(x=col, y='Count') +
        theme(
          axis.text=element_text(size=font_size, color='black'),
          axis.title=element_text(size=font_size, color='black'),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          rect=element_rect(fill="transparent"),
          plot.margin=unit(rep(0.1, 4), "cm"),
          legend.position=c(0,1),
          legend.justification=c(0, 1),
          legend.key.height=unit(font_size, "pt"),
          legend.key.width=unit(font_size, "pt"),
          legend.text = element_text(size=font_size),
        )
      plot_list[[paste0(col, "_", fill_by)]] <- g
    }
  }
  final_plot <- ( plot_list[['within_group_cor_uncorrected_scientific_name']] +
                plot_list[['within_group_cor_corrected_scientific_name']] ) /
                ( plot_list[['within_group_cor_uncorrected_curate_group']] +
                plot_list[['within_group_cor_corrected_curate_group']] )
  ggsave(filename="csca_within_group_cor.pdf", plot=final_plot, width=7.2, height=6.0)
}

extract_selected_tc_only = function(unaveraged_tcs, df_metadata) {
    selected_runs = df_metadata[(df_metadata[['exclusion']]=='no'),'run']
    for (d in c('uncorrected', 'corrected')) {
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            is_selected = colnames(unaveraged_tcs[[d]][[sci_name]]) %in% selected_runs
            unaveraged_tcs[[d]][[sci_name]] = unaveraged_tcs[[d]][[sci_name]][,is_selected]
        }
    }
    return(unaveraged_tcs)
}

unaveraged2averaged = function(unaveraged_tcs, df_metadata, selected_curate_groups) {
    is_curate_groups = list()
    for (curate_group in selected_curate_groups) {
        is_curate_groups[[curate_group]] = (df_metadata[['curate_group']]==curate_group)
    }
    is_sci_names = list()
    for (sci_name in unique(df_metadata[['scientific_name']])) {
        sci_name_ub = sub(' ', '_', sci_name)
        is_sci_names[[sci_name_ub]] = (df_metadata[['scientific_name']]==sci_name)
    }
    is_not_excluded = (df_metadata[['exclusion']]=='no')
    averaged_tcs = list()
    for (d in c('uncorrected', 'corrected')) {
        averaged_tcs[[d]] = list()
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            n_rows = nrow(unaveraged_tcs[[d]][[sci_name]])
            averaged_tcs[[d]][[sci_name]] = data.frame(matrix(ncol = 0, nrow = n_rows))
            rownames(averaged_tcs[[d]][[sci_name]]) = rownames(unaveraged_tcs[[d]][[sci_name]])
            for (curate_group in selected_curate_groups) {
                is_target = (is_curate_groups[[curate_group]] & is_sci_names[[sci_name]] & is_not_excluded)
                target_runs = df_metadata[is_target,'run']
                target_runs = target_runs[target_runs %in% colnames(unaveraged_tcs[[d]][[sci_name]])]
                if (length(target_runs)==0) {
                    next
                }
                label = paste(sub(' ', '_', sci_name), curate_group, sep='_')
                if (sum(is_target)==1) {
                    averaged_tcs[[d]][[sci_name]][,label] = unaveraged_tcs[[d]][[sci_name]][,target_runs]
                } else {
                    averaged_tcs[[d]][[sci_name]][,label] = apply(unaveraged_tcs[[d]][[sci_name]][,target_runs], 1, mean)
                }
            }
            if (ncol(averaged_tcs[[d]][[sci_name]])==0) {
                averaged_tcs[[d]][[sci_name]] = NULL
            }
        }
    }
    return(averaged_tcs)
}

save_group_cor_scatter = function(df_metadata, font_size=8) {
    cat('Generating unaveraged group correlation scatter plot.\n')
    alpha_value = 0.2
    improvement_xymin = 0.5
    improvement_xymax = 2.0
    df_metadata[['corrected_per_uncorrected_group_cor']] = df_metadata[['within_group_cor_corrected']] / df_metadata[['within_group_cor_uncorrected']]
    df_metadata[['corrected_per_uncorrected_max_nongroup_cor']] = df_metadata[['max_nongroup_cor_corrected']] / df_metadata[['max_nongroup_cor_uncorrected']]
    for (col in c('corrected_per_uncorrected_group_cor', 'corrected_per_uncorrected_max_nongroup_cor')) {
        df_metadata[[col]] = ifelse(df_metadata[[col]] < improvement_xymin, improvement_xymin, df_metadata[[col]])
        df_metadata[[col]] = ifelse(df_metadata[[col]] > improvement_xymax, improvement_xymax, df_metadata[[col]])
    }
    ps = list()
    ps[[1]] = ggplot(df_metadata, aes(x = max_nongroup_cor_uncorrected, y = within_group_cor_uncorrected)) + xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[2]] = ggplot(df_metadata, aes(x = max_nongroup_cor_corrected, y = within_group_cor_corrected)) + xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[3]] = ggplot(df_metadata, aes(x = within_group_cor_uncorrected, y = within_group_cor_corrected)) + xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[4]] = ggplot(df_metadata, aes(x = max_nongroup_cor_uncorrected, y = max_nongroup_cor_corrected)) + xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[5]] = ggplot(df_metadata, aes(x = corrected_per_uncorrected_max_nongroup_cor, y = corrected_per_uncorrected_group_cor)) + xlim(c(improvement_xymin, improvement_xymax)) + ylim(c(improvement_xymin, improvement_xymax))
    for (i in 1:length(ps)) {
        ps[[i]] = ps[[i]] + geom_point(alpha=alpha_value)
        ps[[i]] = ps[[i]] + geom_abline(intercept = 0, slope = 1, linetype='dashed', color='blue')
        ps[[i]] = ps[[i]] + theme_bw()
        ps[[i]] = ps[[i]] + geom_density_2d(bins=12, linewidth=0.25, color='gray')
        ps[[i]] = ps[[i]] + theme(
            text=element_text(size=font_size),
            axis.text=element_text(size=font_size),
            axis.title=element_text(size=font_size),
            legend.text=element_text(size=font_size),
            legend.title=element_text(size=font_size)
        )
    }
    ps[[6]] = ggplot() + theme_void()
    combined_plot = wrap_plots(ps)
    ggsave(filename="csca_group_cor_scatter.pdf", plot=combined_plot, width=7.2, height=4.8)
}

df_og = read.table(file_orthogroup, header=TRUE, sep='\t', row.names=1, quote='', check.names=FALSE)
df_gc = read.table(file_genecount, header=TRUE, sep='\t', quote='', check.names=FALSE)
spp_filled = colnames(df_gc)

is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
df_singleog = df_og[is_singlecopy,spp_filled]
spp = sub('_', ' ', spp_filled)
df_metadata = prepare_metadata_table(dir_csca_input_table, selected_curate_groups, spp)
label_orders = get_label_orders(df_metadata)
df_color_averaged = get_df_labels_averaged(df_metadata, label_orders)
df_color_unaveraged = get_df_labels_unaveraged(df_metadata, selected_curate_groups)
cat('Number of orthologs in input table:', nrow(df_og), '\n')
cat('Number of selected single-copy orthologs:', nrow(df_singleog), '\n')
cat('Number of selected species:', length(spp), '\n')

unaveraged_tcs = load_unaveraged_expression_tables(dir_csca_input_table, spp_filled, batch_effect_alg)
unaveraged_tcs = extract_selected_tc_only(unaveraged_tcs, df_metadata)
unaveraged_orthologs = extract_ortholog_unaveraged_expression_table(df_singleog, unaveraged_tcs)
averaged_tcs = unaveraged2averaged(unaveraged_tcs, df_metadata, selected_curate_groups)
averaged_orthologs = extract_ortholog_mean_expression_table(df_singleog, averaged_tcs, label_orders)

cat('Applying expression level imputation for missing orthologs.\n')
imputed_averaged_orthologs = list()
imputed_unaveraged_orthologs = list()
for (d in c('uncorrected','corrected')) {
  imputed_averaged_orthologs[[d]] = impute_expression(averaged_orthologs[[d]])
  imputed_unaveraged_orthologs[[d]] = impute_expression(unaveraged_orthologs[[d]])
  write.table(averaged_orthologs[[d]], file=paste0('csca_ortholog_averaged.',d,'.tsv'), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(unaveraged_orthologs[[d]], file=paste0('csca_ortholog_unaveraged.',d,'.tsv'), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(imputed_averaged_orthologs[[d]], file=paste0('csca_ortholog_averaged.imputed.',d,'.tsv'), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(imputed_unaveraged_orthologs[[d]], file=paste0('csca_ortholog_unaveraged.imputed.',d,'.tsv'), sep="\t", row.names=FALSE, quote=FALSE)
}
cat(nrow(imputed_unaveraged_orthologs[[d]]), 'orthologs were found after filtering and imputation.\n')
df_metadata = calculate_correlation_within_group(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_curate_groups)
save_group_cor_scatter(df_metadata, font_size=8)
save_group_cor_histogram(df_metadata, font_size=8)
save_averaged_tsne_plot(tc=imputed_unaveraged_orthologs[['corrected']], df_label=df_color_unaveraged)
save_averaged_heatmap_plot(imputed_averaged_orthologs, df_color_averaged)
save_averaged_dendrogram_plot(imputed_averaged_orthologs, df_color_averaged)
save_averaged_dimensionality_reduction_summary(imputed_averaged_orthologs, df_color_averaged)
save_averaged_box_plot(imputed_averaged_orthologs, df_color_averaged)
df_metadata = save_unaveraged_pca_plot(imputed_unaveraged_orthologs, df_color_unaveraged, df_metadata)
save_unaveraged_tsne_plot(imputed_unaveraged_orthologs, df_color_unaveraged)

file_metadata_out = file.path(dir_csca, 'metadata.tsv')
write.table(df_metadata, file_metadata_out, row.names=FALSE, sep='\t', quote=FALSE)

if (file.exists('Rplots.pdf')) {
  file.remove('Rplots.pdf')
}
cat('csca.r completed!\n')