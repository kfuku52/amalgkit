if (FALSE) { 
  library(devtools)
  library(httr)
  set_config(config(ssl_verifypeer = 0L))
  options(repos=structure(c(CRAN="http://cran.rstudio.com/")))
  install.packages(c("phytools","amap","RColorBrewer","colorspace","dendextend","NMF","MASS","caper","pvclust"))
}
suppressPackageStartupMessages(library(amap, quietly=TRUE))
suppressPackageStartupMessages(library(RColorBrewer, quietly=TRUE))
suppressPackageStartupMessages(library(colorspace, quietly=TRUE))
suppressPackageStartupMessages(library(dendextend, quietly=TRUE))
suppressPackageStartupMessages(library(NMF, quietly=TRUE))
suppressPackageStartupMessages(library(MASS, quietly=TRUE))
suppressPackageStartupMessages(library(pvclust, quietly=TRUE))
suppressPackageStartupMessages(library(Rtsne, quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly=TRUE))
options(stringsAsFactors = FALSE)

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
font_size = 8

if (debug_mode == "debug") {
  selected_curate_groups = c('root','flower', 'leaf')
  dir_work = '/Users/s229181/MSN/'
  dir_curated_transcriptome = paste0(dir_work, 'curate/')
  dir_tc = paste0(dir_curated_transcriptome,'tables/')
  dir_uncorrected_curate_group_mean = paste0(dir_curated_transcriptome,'tables/')
  dir_curate_group_mean = paste0(dir_curated_transcriptome, 'tables/')
  dir_sra = paste0(dir_curated_transcriptome,'tables/')
  file_outsra = paste0(dir_work,'/sra_table_amalgamated.tsv')
  file_orthogroup = paste0(dir_work, 'OrthoFinder/Orthogroups.tsv')
} else if (debug_mode == "batch") {
  args = commandArgs(trailingOnly = TRUE)
  print(args)
  selected_curate_groups = strsplit(args[1], "\\|")[[1]]
  dir_work = args[2]
  dir_tc = args[3]
  dir_uncorrected_curate_group_mean = args[4]
  dir_curate_group_mean = args[5]
  dir_sra = args[6]
  file_orthogroup = args[7]
  file_genecount = args[8]
  r_util_path = args[9]
  dir_csca = args[10]
}
source(r_util_path)
setwd(dir_csca)

add_color_to_sra = function(df, selected_curate_groups) {
  df = df[,(!colnames(df) %in% c('bp_color','sp_color','curate_group_color'))]
  scientific_name = as.character(df$scientific_name)
  curate_group = as.character(df$curate_group)
  scientific_name_unique = scientific_name[!duplicated(scientific_name)]
  curate_group_unique = curate_group[!duplicated(curate_group)]
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
    df_tmp = rbind(df_tmp, df_label[(df_label$scientific_name==scientific_name)&(df_label$curate_group==curate_group),])
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
  if (length(node)==2) {
    child1_color = attributes(node[[1]])$edgePar[['col']]
    child2_color = attributes(node[[2]])$edgePar[['col']]
    if ((!is.null(child1_color))&(!is.null(child2_color))) {
      if (child1_color==child2_color) {
        attributes(node)$edgePar[['col']] = child1_color
      }
    }
  }
  return(node)
}

map_color = function(redundant_variables, c) {
  uniq_var = unique(redundant_variables)
  uniq_col = rainbow_hcl(length(uniq_var), c=c)
  df_unique = data.frame(var=uniq_var, col=uniq_col, stringsAsFactors=FALSE)
  df_redundant = data.frame(var=redundant_variables, order=seq(1, length(redundant_variables)), stringsAsFactors=FALSE)
  df_redundant = merge(df_redundant, df_unique, by="var", all.x=TRUE, stringsAsFactors=FALSE)
  df_redundant = df_redundant[order(df_redundant$order),]
  return(df_redundant$col)
}

draw_multisp_heatmap = function(tc, df_label) {
  tc_dist_matrix = cor(tc, method='pearson')
  tc_dist_matrix[is.na(tc_dist_matrix)] = 0
  ann_label = df_label[,c('scientific_name','curate_group')]
  colnames(ann_label) = c('species', 'curate_group')
  sp_color = df_label$sp_color[!duplicated(df_label$sp_color)]
  sp_color = sp_color[order(df_label$scientific_name[!duplicated(df_label$scientific_name)])]
  curate_group_color = df_label$curate_group_color[!duplicated(df_label$curate_group_color)]
  curate_group_color = curate_group_color[order(df_label$curate_group[!duplicated(df_label$curate_group)])]
  ann_color = list(species=sp_color, curate_group=curate_group_color)
  breaks = c(0, seq(0.3, 1, 0.01))
  aheatmap(tc_dist_matrix, color="-RdYlBu2:71", Rowv=NA, Colv=NA, revC=TRUE, legend=TRUE, breaks=breaks, 
           annCol=ann_label, annRow=ann_label, annColors=ann_color, annLegend=FALSE, labRow=NA, labCol=NA)
}

draw_multisp_dendrogram = function(tc, df_label, sra, nboot, cex.xlab, cex.yaxis, pvclust_file='pvclust.RData') {
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
  dend_colors = df_label$curate_group_color[order.dendrogram(dend)]
  label_colors = df_label$sp_color[order.dendrogram(dend)]
  labels_colors(dend) = label_colors
  dend_labels <- sra$run[order.dendrogram(dend)]
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
  curate_group_unique = unique(sra$curate_group)
  sp_unique = unique(sra$scientific_name)
  bp_unique = unique(sra$bioproject)
  curate_group_color_unique = unique(sra$curate_group_color)
  sp_color_unique = unique(sra$sp_color)
  bp_color_unique = unique(sra$bp_color)
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
  perplexity = min(30, floor(ncol(tc)/4))
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

file_outsra = file.path(dir_csca, paste0('sra_table_csca.tsv'))

df_og = read.table(file_orthogroup, header=TRUE, sep='\t', row.names=1, quote='')
df_gc = read.table(file_genecount, header=TRUE, sep='\t', check.names=FALSE, quote='')
spp_filled = colnames(df_gc)
is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
df_singleog = df_og[is_singlecopy,spp_filled]
spp_filled = colnames(df_singleog)
spp = sub('_', ' ', spp_filled)
cat('Number of orthologs in input table:', nrow(df_og), '\n')
cat('Number of single-copy orthologs:', nrow(df_singleog), '\n')
cat('Number of selected species:', length(spp), '\n')
cat('Selected species:', spp, '\n')

files = list.files(dir_sra, pattern = ".*sra.*")
df_sra = data.frame()
for (file in files) {
  sra_path = file.path(dir_sra, file)
  tmp_sra = read.table(sra_path, header=TRUE, sep='\t', quote='', comment.char='')
  df_sra = rbind(df_sra, tmp_sra)
}
write.table(df_sra, file_outsra, row.names=FALSE, sep='\t')
df_sra = df_sra[(df_sra[['curate_group']] %in% selected_curate_groups)&(df_sra[['scientific_name']] %in% spp),]
order_cg = order(df_sra[['curate_group']])
label_orders = unique(paste(df_sra[order_cg,'scientific_name'], df_sra[order_cg,'curate_group'], sep='_'))
label_orders = sub(' ', '_', label_orders)

tcs = list()
tcs[['uncorrected']] = list()
tcs[['corrected']] = list()

all_files = list.files(dir_curate_group_mean, pattern = "*.curate_group.mean.tsv")
uncorrected_files = all_files[grep("uncorrected", all_files)]
corrected_files = all_files[-grep("uncorrected", all_files)]
for (sp in spp_filled) {
  uncorrected_file = uncorrected_files[startsWith(uncorrected_files, sp)]
  corrected_file = corrected_files[startsWith(corrected_files, sp)]
  uncorrected_path = file.path(dir_uncorrected_curate_group_mean, uncorrected_file)
  corrected_path = file.path(dir_curate_group_mean, corrected_file)
  tcs[['uncorrected']][[sp]] = read.delim(uncorrected_path, header=TRUE, row.names=1, sep='\t')
  tcs[['corrected']][[sp]] = read.delim(corrected_path, header=TRUE, row.names=1, sep='\t')
}

# Orthogroup mean expression

ortholog = list()
ortholog[['uncorrected']] = df_singleog
ortholog[['corrected']] = df_singleog
row_names = rownames(ortholog[['corrected']])
selected_species = spp
for (d in c('uncorrected', 'corrected')) {
  cat(d, '\n')
  for (sp in spp) {
    sp_filled = sub(" ", "_", sp)
    tc = tcs[[d]][[sp_filled]]
    colnames(tc) = paste(sp_filled, colnames(tc), sep='_')
    ortholog[[d]] = merge(ortholog[[d]], tc, by.x=sp_filled, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
  }
  num_remove_col = length(selected_species)
  ortholog[[d]] = ortholog[[d]][,-(1:num_remove_col)]
  rownames(ortholog[[d]]) = row_names
  ortholog[[d]] = sort_averaged_tc(ortholog[[d]])
  ortholog[[d]] = ortholog[[d]][, label_orders]
}
cat(nrow(ortholog[[d]]), 'orthologs were found before filtering.\n')

unaveraged_tcs = list()
unaveraged_tcs[['corrected']] = list()
all_files = list.files(dir_curate_group_mean, pattern = "*.tc.*")
uncorrected_files = all_files[grep("uncorrected", all_files)]
corrected_files = all_files[-grep("uncorrected", all_files)]
for (sp in spp_filled) {
  corrected_file = corrected_files[startsWith(corrected_files, sp)]
  corrected_path = file.path(dir_tc, corrected_file)
  unaveraged_tcs[['corrected']][[sp]] = read.delim(corrected_path, header=TRUE, row.names=1, sep='\t')
}
# Unaveraged orthogroup extraction
orthogroup = list()
is_complete = apply(df_og!='', 1, all)
ogtmp = df_og[is_complete,]
orthogroup[['uncorrected']] = ogtmp[,NULL]
orthogroup[['corrected']] = ogtmp[,NULL]
row_names = rownames(orthogroup[['corrected']])
selected_species = spp
for (d in c('corrected')) {
  cat(d, '\n')
  for (sp in spp) {
    cat(sp, '\n')
    sp_filled = sub(" ", "_", sp)
    tc = unaveraged_tcs[[d]][[sp_filled]]
    tc2 = tc[NULL,]
    for (og in row_names) {
      genes = strsplit(ogtmp[og,sp_filled], ', ')[[1]]
      #genes = sub('.*_','',genes)
      if (length(genes)==1) {
        tc2[og,] = tc[genes,]
      } else if (length(genes)>1) {
        tc2[og,] = apply(tc[genes,], 2, median)
      }
    }
    colnames(tc2) = paste(sp_filled, colnames(tc2), sep='_')
    orthogroup[[d]] = merge(orthogroup[[d]], tc2, by="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
    rownames(orthogroup[[d]]) = orthogroup[[d]][['Row.names']]
    orthogroup[[d]][['Row.names']] = NULL
  }
  num_remove_col = length(selected_species)
  orthogroup[[d]] = orthogroup[[d]][,-(1:num_remove_col)]
  rownames(orthogroup[[d]]) = row_names
  orthogroup[[d]] = sort_averaged_tc(orthogroup[[d]])
  #orthogroup[[d]][orthogroup[[d]]<0] = 0
  #orthogroup[[d]] = orthogroup[[d]][, label_orders]
}

# Single-copy unaveraged_ortholog extraction

unaveraged_ortholog = list()
unaveraged_ortholog[['uncorrected']] = df_singleog
unaveraged_ortholog[['corrected']] = df_singleog
row_names = rownames(unaveraged_ortholog[['corrected']])
selected_species = spp
for (d in c('corrected')) {
  cat(d, '\n')
  for (sp in spp) {
    sp_filled = sub(" ", "_", sp)
    tc = unaveraged_tcs[[d]][[sp_filled]]
    colnames(tc) = paste(sp_filled, colnames(tc), sep='_')
    unaveraged_ortholog[[d]] = merge(unaveraged_ortholog[[d]], tc, by.x=sp_filled, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
  }
  num_remove_col = length(selected_species)
  unaveraged_ortholog[[d]] = unaveraged_ortholog[[d]][,-(1:num_remove_col)]
  rownames(unaveraged_ortholog[[d]]) = row_names
  is_not_na = apply(unaveraged_ortholog[[d]], 1, function(x){!any(is.na(x))})
  unaveraged_ortholog[[d]] = unaveraged_ortholog[[d]][is_not_na,]
  #unaveraged_ortholog[[d]] = sort_averaged_tc(unaveraged_ortholog[[d]])
  #unaveraged_ortholog[[d]][unaveraged_ortholog[[d]]<0] = 0
  #unaveraged_ortholog[[d]] = unaveraged_ortholog[[d]][, label_orders]
}
cat(nrow(unaveraged_ortholog[[d]]), 'unaveraged_orthologs were found after filtering.\n')

cols = c('run','bioproject','curate_group','scientific_name','sp_color','curate_group_color','bp_color')
df_label2 = add_color_to_sra(df_sra[(df_sra$exclusion=='no'),], selected_curate_groups)
df_label2 =df_label2[,cols]
label_order = order(df_label2[['run']])
df_label2 = df_label2[label_order,]
tc_order = order(sub('.*_','',colnames(unaveraged_ortholog[['corrected']])))
unaveraged_ortholog[['corrected']] = unaveraged_ortholog[['corrected']][,tc_order]

draw_multisp_tsne2 = function(tc, df_label) {
  perplexity = min(30, floor(ncol(tc)/4))
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

draw_multisp_tsne2(tc=unaveraged_ortholog[['corrected']], df_label=df_label2)

df_labels = list()
for (tpm in c('uncorrected','corrected')) {
  if (tpm=='uncorrected') {
    sra_tmp = df_sra
  } else if (tpm=='corrected') {
    sra_tmp = df_sra[(df_sra$exclusion=='no'),]
  }
  
  df_label = unique(sra_tmp[,c('scientific_name','curate_group')])
  categories = list(scientific_name=sra_tmp$scientific_name, curate_group=sra_tmp$curate_group)
  
  df_bp = aggregate(sra_tmp$bioproject, by=categories, function(x){length(unique(x))})
  colnames(df_bp) = c('scientific_name','curate_group','num_bp')
  df_label = merge(df_label, df_bp, all.x=TRUE, all.y=FALSE)
  
  df_exp = aggregate(sra_tmp$experiment, by=categories, function(x){length(unique(x))})
  colnames(df_exp) = c('scientific_name','curate_group','num_exp')
  df_label = merge(df_label, df_exp, all.x=TRUE, all.y=FALSE)
  df_label = df_label[order(df_label$curate_group, df_label$scientific_name),]
  df_label = sort_labels(df_label, label_orders)
  df_label = add_color_to_sra(df_label, selected_curate_groups)
  df_label = sort_labels(df_label, label_orders)
  df_labels[[tpm]] = df_label
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

max_gene_per_og = 1
cat('Maximum number of genes per orthogroup =', max_gene_per_og, '\n')
preprocessing = 'no_averaging'
tc_in = orthogroup[['corrected']]
tc_in = tc_in[apply(tc_in, 1, function(x){!any(is.na(x))}),]
dfl_in = df_label2

if (preprocessing=='no_averaging') {
  by = 'run'
} else if (preprocessing=='bioproject_mean') {
  dfl_in[['bioproject']] = sub(';.*', '', dfl_in[['bioproject']])
  for (sp in unique(dfl_in[['scientific_name']])) {
    for (curate_group in unique(dfl_in[['curate_group']])) {
      dfl_tmp = dfl_in[(dfl_in[['scientific_name']]==sp)&(dfl_in[['curate_group']]==curate_group),]
      for (bp in unique(dfl_tmp[['bioproject']])) {
        key = paste0(sub(' ','_',sp), '_', curate_group, '_', bp)
        runs = dfl_tmp[(dfl_tmp[['bioproject']]==bp),'run']
        cols = paste0(sub(' ','_',sp), '_', runs)
        if (length(cols)==1) {
          mean_values = tc_in[,cols]
        } else {
          mean_values = apply(tc_in[,cols], 1, mean)
        }
        tc_in[,cols] = NULL
        tc_in[[key]] = mean_values
      }
    }
  }
  dfl_in[['run']] = NULL
  dfl_in = unique(dfl_in)
  dfl_in[['key']] = paste0(sub(' ', '_', dfl_in[['scientific_name']]), '_', dfl_in[['curate_group']], '_', dfl_in[['bioproject']])
  by='key'
} else if (preprocessing=='curate_group_mean') {
  by='species_curate_group'
}
out = get_pca_coordinates(tc=tc_in, df_label=dfl_in, by='run')
tmp = out[[1]]
pc_contributions = out[[2]]

cat('Number of orthogroups:', nrow(tc_in), '\n')

for (pcxy in list(c(1,2),c(3,4))) {    
  pcx = pcxy[1]
  pcy = pcxy[2]
  
  colx = paste0('PC', pcx)
  coly = paste0('PC', pcy)
  xmin = min(tmp[[colx]])
  xmax = max(tmp[[colx]])
  xunit = (xmax-xmin)*0.01
  xmin = xmin - xunit
  xmax = xmax + xunit
  
  ymin = min(tmp[[coly]])
  ymax = max(tmp[[coly]])
  yunit = (ymax-ymin)*0.01
  ymin = ymin - yunit
  ymax = ymax + yunit

  write.table(df_label, 'df_label.tsv', sep='\t', row.names=FALSE)

  g = ggplot(tmp, aes(x=!!rlang::sym(colx), y=!!rlang::sym(coly), color=curate_group))
  g = g + theme_bw()
  g = g + geom_point(size=0.5)
  # g = g + geom_density_2d(mapping=aes(color=curate_group), bins=12, size=0.25) # Temporarily deactivated. Bug fix needed.
  curate_group_colors = unique(df_label[,c('curate_group','curate_group_color')])[['curate_group_color']]
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
  filename = paste0('multisp_pca_maxGenePerOG', max_gene_per_og, '_PC', pcx, pcy, '.pdf')
  ggsave(file=filename, g, height=2.15, width=4.25)
  g
}

preprocessing = 'no_averaging'
tc_in = orthogroup[['corrected']]
tc_in = tc_in[apply(tc_in, 1, function(x){!any(is.na(x))}),]
dfl_in = df_label2

if (preprocessing=='no_averaging') {
  by = 'run'
} else if (preprocessing=='bioproject_mean') {
  dfl_in[['bioproject']] = sub(';.*', '', dfl_in[['bioproject']])
  for (sp in unique(dfl_in[['scientific_name']])) {
    for (curate_group in unique(dfl_in[['curate_group']])) {
      dfl_tmp = dfl_in[(dfl_in[['scientific_name']]==sp)&(dfl_in[['curate_group']]==curate_group),]
      for (bp in unique(dfl_tmp[['bioproject']])) {
        key = paste0(sub(' ','_',sp), '_', curate_group, '_', bp)
        runs = dfl_tmp[(dfl_tmp[['bioproject']]==bp),'run']
        cols = paste0(sub(' ','_',sp), '_', runs)
        if (length(cols)==1) {
          mean_values = tc_in[,cols]
        } else {
          mean_values = apply(tc_in[,cols], 1, mean)
        }
        tc_in[,cols] = NULL
        tc_in[[key]] = mean_values
      }
    }
  }
  dfl_in[['run']] = NULL
  dfl_in = unique(dfl_in)
  dfl_in[['key']] = paste0(sub(' ', '_', dfl_in[['scientific_name']]), '_', dfl_in[['curate_group']], '_', dfl_in[['bioproject']])
  by='key'
} else if (preprocessing=='curate_group_mean') {
  by='species_curate_group'
}

get_tsne_coordinates = function(tc, df_label, by='run') {
  perplexity = min(30, floor(ncol(tc)/4))
  set.seed(1)
  out_tsne = Rtsne(as.matrix(t(tc)), theta=0, check_duplicates=FALSE, verbose=FALSE, perplexity=perplexity, dims=2)
  tmp = data.frame(tsne1=out_tsne$Y[,1], tsne2=out_tsne$Y[,2])
  tmp[[by]] = sub('.*_', '', colnames(tc))
  tmp = merge(df_label, tmp, by=by)
  return(tmp)
}

tmp = get_tsne_coordinates(tc=tc_in, df_label=dfl_in)

pcx = 1
pcy = 2
colx = paste0('tsne', pcx)
coly = paste0('tsne', pcy)
xmin = min(tmp[[colx]])
xmax = max(tmp[[colx]])
xunit = (xmax-xmin)*0.01
xmin = xmin - xunit
xmax = xmax + xunit

ymin = min(tmp[[coly]])
ymax = max(tmp[[coly]])
yunit = (ymax-ymin)*0.01
ymin = ymin - yunit
ymax = ymax + yunit

g = ggplot(tmp, aes(x=!!rlang::sym(colx), !!rlang::sym(coly), color=curate_group))
g = g + theme_bw()
g = g + geom_point(size=0.5)
#g = g+ geom_density_2d(mapping=aes(color=curate_group), bins=12, size=0.25) # Temporarily deactivated. Bug fix needed.
curate_group_colors = unique(df_label[,c('curate_group','curate_group_color')])[['curate_group_color']]
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
filename = paste0('multisp_tsne.pdf')
ggsave(file=filename, g, height=2.15, width=4.25)
g 

do_heatmap=TRUE
do_dendrogram=TRUE
do_pca_mds=TRUE

if (do_heatmap) {
  file_name='Multisp.SVA.heatmap.pdf'
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
  for (tpm in c('uncorrected','corrected')) {
    tc = ortholog[[tpm]]
    df_label = df_labels[[tpm]]
    par(mar=c(0,0,0,0))
    draw_multisp_heatmap(tc=tc, df_label=df_label)
  }
  graphics.off()
}

if (do_dendrogram) {
  file_name='Multisp.SVA.dendrogram.pdf'
  pdf(file_name, height=2.5, width=7.2) # full figure size = 9.7 x 7.2
  layout_matrix=matrix(c(1,2),2,1,byrow=TRUE)
  layout(t(layout_matrix))
  for (tpm in c('uncorrected','corrected')) {
    tc = ortholog[[tpm]]
    df_label = df_labels[[tpm]]
    par(cex=0.5, mar=c(10,5.5,0,0), mgp=c(4, 0.7, 0))
    pvclust_file=paste0('Multisp.pvclust.',tpm,'.RData')
    draw_multisp_dendrogram(tc=tc, df_label=df_label, sra=df_sra, pvclust_file=pvclust_file, 
                            nboot=1, cex.xlab=0.3, cex.yaxis=0.5)
  }
  graphics.off()
}

if (do_pca_mds) {
  par(cex=1)
  cex_axis=0.7
  cex_plot=1
  file_name='Multisp.SVA.PCA.MDS.pdf'
  pdf(file_name, height=4.5, width=7.2) # full figure size = 9.7 x 7.2
  layout_matrix = matrix(c(1,1,1,3,3,3,5,5,2,2,2,4,4,4,5,5),2,8,byrow=TRUE)
  layout(layout_matrix)
  for (tpm in c('uncorrected','corrected')) {
    tc = ortholog[[tpm]]
    df_label = df_labels[[tpm]]
    par(mar=c(4,4,0.1,1)); draw_multisp_pca(tc=tc, df_label=df_label)
    par(mar=c(4,4,0.1,1)); draw_multisp_tsne(tc=tc, df_label=df_label)
    #par(mar=c(4,4,0.1,1)); draw_multisp_mds(tc=tc, df_label=df_label)
  }
  par(mar=c(0,0,0,0)); draw_multisp_legend(df_label)
  graphics.off()
}
savefig=TRUE
do_combined=TRUE
panel_cex=1
chars = c('A','B','C','D','E','G','F','H')

if (do_combined) {
  if (savefig) {
    file_name='Multisp.SVA.combined.pdf'
    pdf(file_name, height=9.7, width=7.2) # full figure size = 9.7 x 7.2        
  }　else {
    options(repr.plot.height=9.7, repr.plot.width=7.2)
  }
  layout_matrix=matrix(c(
    1,1,1,1,2,2,2,2,
    3,3,3,3,4,4,4,4,
    3,3,3,3,4,4,4,4,
    3,3,3,3,4,4,4,4,
    3,3,3,3,4,4,4,4,
    3,3,3,3,4,4,4,4,
    3,3,3,3,4,4,4,4,
    3,3,3,3,4,4,4,4,
    3,3,3,3,4,4,4,4,
    5,5,5,5,6,6,6,6,
    5,5,5,5,6,6,6,6,
    5,5,5,5,6,6,6,6,
    5,5,5,5,6,6,6,6,
    5,5,5,5,6,6,6,6,
    5,5,5,5,6,6,6,6,
    7,7,7,9,9,9,11,11,
    7,7,7,9,9,9,11,11,
    7,7,7,9,9,9,11,11,
    7,7,7,9,9,9,11,11,
    7,7,7,9,9,9,11,11,
    7,7,7,9,9,9,11,11,
    8,8,8,10,10,10,11,11,
    8,8,8,10,10,10,11,11,
    8,8,8,10,10,10,11,11,
    8,8,8,10,10,10,11,11,
    8,8,8,10,10,10,11,11,
    8,8,8,10,10,10,11,11),
    27,8,byrow=TRUE)
  layout(layout_matrix)
  par(cex=1, xpd=T, mar=c(0,0,0,0)); plot(c(0, 1), c(0, 1), ann=F, bty='n',type='n',xaxt='n', yaxt='n')
  text(0.50,0.31,'Uncorrected', srt=0, cex=1, font=1, adj=c(0.5,1))
  text(-0.04, 0.31, chars[1], srt=0, font=2, cex=panel_cex, adj=c(0,1)) ; chars = chars[2:length(chars)]
  par(cex=1, xpd=T, mar=c(0,0,0,0)); plot(c(0, 1), c(0, 1), ann=F, bty='n',type='n',xaxt='n', yaxt='n')
  text(0.50,0.31,'Corrected', srt=0, cex=1, font=1, adj=c(0.5,1))
  text(-0.04, 0.31, chars[1], srt=0, font=2, cex=panel_cex, adj=c(0,1)) ; chars = chars[2:length(chars)]
  for (tpm in c('uncorrected','corrected')) {
    tc = ortholog[[tpm]]
    df_label = df_labels[[tpm]]
    par(mar=c(0,0,0,0))
    draw_multisp_heatmap(tc=tc, df_label=df_label)
    text(0.01, 0.93, 'Sp', srt=90, font=1, cex=0.7)
    text(0.05, 0.93, 'Tis', srt=90, font=1, cex=0.7)
    text(0.805, 0.97, 'Sp', srt=0, font=1, cex=0.7)
    text(0.805, 0.93, 'Tis', srt=0, font=1, cex=0.7)
  }
  for (tpm in c('uncorrected','corrected')) {
    tc = ortholog[[tpm]]
    df_label = df_labels[[tpm]]
    par(cex=0.5, mar=c(11,5.5,0,0), mgp=c(4, 0.7, 0))
    pvclust_file=paste0('Multisp.pvclust.',tpm,'.RData')
    draw_multisp_dendrogram(tc=tc, df_label=df_label, sra=df_sra, pvclust_file=pvclust_file, nboot=1000, cex.xlab=0.3, cex.yaxis=0.5)
    par(new=T, cex=1, xpd=T, mar=c(0,0,0,0)); plot(c(0, 1), c(0, 1), ann=F, bty='n',type='n',xaxt='n', yaxt='n')
    text(-0.04, 1.03, chars[1], srt=0, font=2, cex=panel_cex, adj=c(0,1)) ; chars = chars[2:length(chars)]
  }
  par(mar=c(0,0,0,0))
  for (tpm in c('uncorrected','corrected')) {
    tc = ortholog[[tpm]]
    df_label = df_labels[[tpm]]
    par(mar=c(4,4,0.1,1)); draw_multisp_pca(tc=tc, df_label=df_label)
    par(new=T, cex=1, xpd=T, mar=c(0,0,0,0)); plot(c(0, 1), c(0, 1), ann=F, bty='n',type='n',xaxt='n', yaxt='n')
    text(-0.04, 1.03, chars[1], srt=0, font=2, cex=panel_cex, adj=c(0,1)) ; chars = chars[2:length(chars)]
    par(mar=c(4,4,0.1,1)); draw_multisp_tsne(tc=tc, df_label=df_label)
    par(new=T, cex=1, xpd=T, mar=c(0,0,0,0)); plot(c(0, 1), c(0, 1), ann=F, bty='n',type='n',xaxt='n', yaxt='n')
    text(-0.04, 1.03, chars[1], srt=0, font=2, cex=panel_cex, adj=c(0,1)) ; chars = chars[2:length(chars)]
  }
  par(mar=c(0,0,0,0), xpd=T); draw_multisp_legend(df_label)
  if (savefig) {
    graphics.off()
  }
}

draw_multisp_boxplot = function(sra, tc_dist_matrix, fontsize=7) {
  is_same_sp = outer(sra$scientific_name, sra$scientific_name, function(x,y){x==y})
  is_same_curate_group = outer(sra$curate_group, sra$curate_group, function(x,y){x==y})
  plot(c(0.5, 4.5), c(0, 1), type = 'n', xlab='', ylab="Pearson's correlation\ncoefficient", las=1, xaxt='n')
  boxplot(tc_dist_matrix[(!is_same_sp)&(!is_same_curate_group)], at=1, add=TRUE, col='gray', yaxt='n')
  boxplot(tc_dist_matrix[(is_same_sp)&(!is_same_curate_group)], at=2, add=TRUE, col='gray', yaxt='n')
  boxplot(tc_dist_matrix[(!is_same_sp)&(is_same_curate_group)], at=3, add=TRUE, col='gray', yaxt='n')
  boxplot(tc_dist_matrix[(is_same_sp)&(is_same_curate_group)], at=4, add=TRUE, col='gray', yaxt='n')
  labels = c('bw\nbw', 'bw\nwi', 'wi\nbw', 'wi\nwi')
  axis(side=1, at=c(1,2,3,4), labels=labels, padj=0.5)
  axis(side=1, at=0.35, labels='Group\nSpecies', padj=0.5, hadj=1, tick=FALSE)
}

file_name='Multisp.boxplot.pdf'
pdf(file_name, height=3.6, width=7.2) # full figure size = 9.7 x 7.2
par(mfrow=c(1,2))
for (tpm in c('uncorrected','corrected')) {
  tc = ortholog[[tpm]]
  tc[tc<0] = 0
  tc_dist_matrix = cor(tc, method='pearson')
  tc_dist_matrix[is.na(tc_dist_matrix)] = 0
  draw_multisp_boxplot(df_labels[[tpm]], tc_dist_matrix, fontsize=7)
}
tpm = 'corrected'
tc = ortholog[[tpm]]
df_label = df_labels[[tpm]]
par(mar=c(2.5,2.5,0.3,0.1), cex=1, ps=8, mgp=c(1.5, 0.7, 0)); draw_multisp_tsne(tc=tc, df_label=df_label$sp_color)
graphics.off()
if (file.exists('Rplots.pdf')) {
  file.remove('Rplots.pdf')
}