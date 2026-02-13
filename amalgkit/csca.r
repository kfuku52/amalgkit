#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(amap, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(colorspace, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(dendextend, quietly = TRUE)))
detected_cores = tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
if (is.na(detected_cores) && (is.null(getOption("cores")) || is.na(getOption("cores")))) {
    options(cores = 1L)
}
suppressWarnings(suppressPackageStartupMessages(library(NMF, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(MASS, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pvclust, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(Rtsne, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pcaMethods, quietly = TRUE)))
options(stringsAsFactors = FALSE)

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
font_size = 8

if (debug_mode == "debug") {
    developer = 'kf'
    if (developer == 'mf') {
        selected_sample_groups = c('root', 'flower', 'leaf')
        dir_work = '/home/s229181/projects/amalgkit_paper/Plant_set'
        #selected_sample_groups = c('Amicoumacin_low','Amicoumacin_high','anaerobic','Azithromycin_low','Azithromycin_high','control','ciprofloxacin_low','ciprofloxacin_high','colistin_low','colistin_high','H2O2','heat','meropenem_low','meropenem_high','NaCl','Anaerobic','IPTG','biofilm_medium','acid','H202','butyrate','aerobic','Ampicillin','Vancomycin','ciprofloxacin','colistin','glucose','Glucose','rifampicin','probiotic','cleaner','ampicillin','tetracycline','pediocin','glycerol','Pyruvate','Ca2+','Glycerol','H2o2','anhydrotetracycline','TB47','treated','iron','Lactate','rion','phage','Ag','biofilm','MIC','AZM','citrate','NaNO2','Acetate','sucrose','coumermycin','copper','mitomycin','arabinose','Cefotaxime','Cellulose','vancomycin','mupirocin','galactose','macrophages','tobramycin')
        sample_group_colors = 'DEFAULT'
        #dir_work = '/home/s229181/projects/amalgkit_paper/Prokaryote_set'
        dir_csca_input_table = file.path(dir_work, 'csca/csca_input_symlinks')
        file_orthogroup = file.path(dir_work, 'csca/multispecies_busco_table.tsv')
        file_genecount = file.path(dir_work, 'csca/multispecies_genecount.tsv')
        r_util_path = '/home/s229181/projects/amalgkit_paper/amalgkit/amalgkit/util.r'
        dir_csca = file.path(dir_work, 'csca')
        batch_effect_alg = 'sva'
    } else if (developer == 'kf') {
        selected_sample_groups = c('root', 'flower', 'leaf')
        sample_group_colors = c('#d95f02ff', '#1b9e77ff', '#7570b3ff')
        dir_work = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
        dir_csca_input_table = file.path(dir_work, 'csca/csca_input_symlinks')
        file_orthogroup = file.path(dir_work, 'csca/multispecies_busco_table.tsv')
        file_genecount = file.path(dir_work, 'csca/multispecies_genecount.tsv')
        r_util_path = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/repos/amalgkit/amalgkit/util.r'
        dir_csca = file.path(dir_work, 'csca')
        batch_effect_alg = 'sva'
    }
} else if (debug_mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    selected_sample_groups = strsplit(args[1], "\\|")[[1]]
    sample_group_colors = strsplit(args[2], ",")[[1]]
    dir_work = args[3]
    dir_csca_input_table = args[4]
    file_orthogroup = args[5]
    file_genecount = args[6]
    r_util_path = args[7]
    dir_csca = args[8]
    batch_effect_alg = args[9]
}
source(r_util_path)
setwd(dir_csca)
cat('selected_sample_groups:', selected_sample_groups, "\n")
cat('selected_sample_group_colors:', sample_group_colors, "\n")
cat('dir_work:', dir_work, "\n")
cat('dir_csca_input_table:', dir_csca_input_table, "\n")
cat('file_orthogroup:', file_orthogroup, "\n")
cat('file_genecount:', file_genecount, "\n")
cat('r_util_path:', r_util_path, "\n")
cat('dir_csca:', dir_csca, "\n")
cat('batch_effect_alg:', batch_effect_alg, "\n")

sort_labels = function(df_label, label_orders) {
    df_tmp = data.frame()
    for (lo in label_orders) {
        splits = strsplit(lo, '_')[[1]]
        scientific_name = paste(splits[1], splits[2])
        sample_group = paste(splits[3:length(splits)], collapse = '_')
        df_tmp = rbind(df_tmp, df_label[(df_label[['scientific_name']] == scientific_name) & (df_label[['sample_group']] == sample_group),])
    }
    return(df_tmp)
}

sort_averaged_tc = function(tc) {
    split_colnames = strsplit(colnames(tc), "_")
    genus_names = c()
    specific_names = c()
    sample_group_names = c()
    for (i in 1:length(split_colnames)) {
        genus_names = c(genus_names, split_colnames[[i]][1])
        specific_names = c(specific_names, split_colnames[[i]][2])
        sample_group_names = c(sample_group_names, paste0(split_colnames[[i]][3:length(split_colnames[[i]])], collapse = '_'))
    }
    colname_order = order(sample_group_names, genus_names, specific_names)
    tc = tc[, colname_order]
    return(tc)
}

color_children2parent = function(node) {
    if (length(node) != 2) {
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
    if (child1_color == child2_color) {
        attributes(node)$edgePar[['col']] = child1_color
    }
    return(node)
}

map_color = function(redundant_variables, c) {
    uniq_var = unique(redundant_variables)
    uniq_col = rainbow_hcl(length(uniq_var), c = c)
    df_unique = data.frame(var = uniq_var, col = uniq_col, stringsAsFactors = FALSE)
    df_redundant = data.frame(var = redundant_variables, order = seq(1, length(redundant_variables)), stringsAsFactors = FALSE)
    df_redundant = merge(df_redundant, df_unique, by = "var", all.x = TRUE, stringsAsFactors = FALSE)
    df_redundant = df_redundant[order(df_redundant[['order']]),]
    return(df_redundant[['col']])
}

draw_multisp_heatmap = function(tc, df_label) {
    tc_dist_matrix = cor(tc, method = 'pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    ann_label = df_label[, c('scientific_name', 'sample_group')]
    colnames(ann_label) = c('species', 'sample_group')
    is_sp_first_appearance = (!duplicated(df_label[['sp_color']]))
    sp_color = df_label[is_sp_first_appearance, 'sp_color']
    sp_color = sp_color[order(df_label[is_sp_first_appearance, 'scientific_name'])]
    is_cg_first_appearance = (!duplicated(df_label[['sample_group_color']]))
    sample_group_color = df_label[is_cg_first_appearance, 'sample_group_color']
    sample_group_color = sample_group_color[order(df_label[is_cg_first_appearance, 'sample_group'])]
    ann_color = list(species = sp_color, sample_group = sample_group_color)
    breaks = c(0, seq(0.3, 1, 0.01))
    NMF::aheatmap(tc_dist_matrix, color = "-RdYlBu2:71", Rowv = NA, Colv = NA, revC = TRUE, legend = TRUE, breaks = breaks,
                  annCol = ann_label, annRow = ann_label, annColors = ann_color, annLegend = FALSE, labRow = NA, labCol = NA)
}

draw_multisp_dendrogram = function(tc, df_label, df_metadata, nboot, cex.xlab, cex.yaxis, pvclust_file = 'pvclust.RData') {
    colnames(tc) = sub("_.*", "", sub('_', ' ', colnames(tc)))
    dist_fun = function(x) { Dist(t(x), method = 'pearson') }
    if (file.exists(pvclust_file)) {
        if (file.info(pvclust_file)$size) {
            cat('pvclust file found.\n')
            load(pvclust_file)
        }
    } else {
        cat('No pvclust file found. Start bootstrapping.\n')
        result = pvclust(tc, method.dist = dist_fun, method.hclust = "average", nboot = nboot, parallel = FALSE) # UPGMA
        save(result, file = pvclust_file)
    }
    dend = as.dendrogram(result)
    dend_colors = df_label[order.dendrogram(dend), 'sample_group_color']
    label_colors = df_label[order.dendrogram(dend), 'sp_color']
    labels_colors(dend) = label_colors
    dend_labels <- df_metadata[order.dendrogram(dend), 'run']
    dend <- color_branches(dend, labels = dend_labels, col = dend_colors)
    dend <- set(dend, "branches_lwd", 2)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    par(cex = cex.xlab)
    plot(dend, las = 1, yaxt = 'n')
    text(result, print.num = FALSE, cex = 1, col.pv = 'black')
    par(cex = cex.yaxis)
    axis(2, las = 1)
    mtext(text = 'Distance', side = 2, line = 4, cex = cex.yaxis)

    n = ncol(tc)
    f = 100
    sample_group_unique = unique(df_metadata['sample_group'])
    sp_unique = unique(df_metadata[['scientific_name']])
    bp_unique = unique(df_metadata[['bioproject']])
    sample_group_color_unique = unique(df_metadata[['sample_group_color']])
    sp_color_unique = unique(df_metadata[['sp_color']])
    bp_color_unique = unique(df_metadata[['bp_color']])
    legend_text = c(as.character(sample_group_unique), "", as.character(sp_unique), "", as.character(bp_unique))
    legend_bg = c(sample_group_color_unique, "white", sp_color_unique, "white", bp_color_unique)
    legend_fg = c(rep("black", length(sample_group_color_unique)), "white", rep("black", length(sp_color_unique)), "white", rep("black", length(bp_color_unique)))
    #plot.new() ; par(mar=c(0,0,0,0))
    #legend("center", legend=legend_text, cex=1, pch=22, lty=0, lwd=1, pt.bg=legend_bg, col=legend_fg)
}

draw_multisp_pca = function(tc, df_label) {
    tc_dist_matrix = cor(tc, method = 'pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    set.seed(1)
    pca = prcomp(tc_dist_matrix)
    xlabel = paste0("PC 1 (", round(summary(pca)$importance[2, 1] * 100, digits = 1), "%)")
    ylabel = paste0("PC 2 (", round(summary(pca)$importance[2, 2] * 100, digits = 1), "%)")
    plot(
        pca[['x']][, 1],
        pca[['x']][, 2],
        pch = 21,
        cex = 2,
        lwd = 1,
        bg = df_label[['sample_group_color']],
        col = df_label[['sp_color']],
        xlab = xlabel,
        ylab = ylabel,
        las = 1
    )
}

draw_multisp_mds = function(tc, df_label) {
    tc_dist_dist = Dist(t(tc), method = 'pearson') + 0.000000001
    tc_dist_dist[is.na(tc_dist_dist)] = 1
    set.seed(1)
    try_out = tryCatch(
    { isoMDS(tc_dist_dist, k = 2, maxit = 100) },
    error = function(a) { return("MDS failed.") }
    )
    if (mode(try_out) == "character") {
        cat('MDS failed.\n')
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    } else {
        mds <- try_out
        plot(mds$points[, 1], mds$points[, 2], pch = 21, cex = 2, lwd = 1, bg = df_label$sample_group_color, col = df_label$sp_color, xlab = "MDS dimension 1", ylab = "MDS dimension 2", las = 1)
    }
}

draw_multisp_tsne = function(tc, df_label) {
    num_samples = ncol(tc)
    perplexity = min(30, floor(num_samples / 4), na.rm = TRUE)
    if ((num_samples < 4) || (!is.finite(perplexity)) || (perplexity <= 0)) {
        cat('t-SNE skipped: insufficient number of samples.\n')
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        return(invisible(NULL))
    }
    set.seed(1)
    out_tsne = Rtsne(as.matrix(t(tc)), theta = 0, check_duplicates = FALSE, verbose = FALSE, perplexity = perplexity, dims = 2)
    try_out = tryCatch(
    {
        plot(out_tsne$Y[, 1], out_tsne$Y[, 2], pch = 21, cex = 2, lwd = 1, bg = df_label$sample_group_color, col = df_label$sp_color,
             xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", las = 1)
    },
        error = function(a) { return("t-SNE plot failed.") }
    )
    if (mode(try_out) == "character") {
        cat('t-SNE failed.\n')
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    }
}

draw_multisp_legend = function(df_label) {
    cex_axis = 0.7
    sample_group_unique = df_label$sample_group[!duplicated(df_label$sample_group)]
    sp_unique = df_label$scientific_name[!duplicated(df_label$scientific_name)]
    sample_group_color_unique = df_label$sample_group_color[!duplicated(df_label$sample_group_color)]
    sp_color_unique = df_label$sp_color[!duplicated(df_label$sp_color)]
    toumei = rgb(1, 1, 1, 0)
    legend_text = c('Sample group', as.character(sample_group_unique), "", 'Species', as.character(sp_unique))
    legend_bg = c(toumei, sample_group_color_unique, toumei, toumei, rep(toumei, length(sp_color_unique)))
    legend_fg = c(toumei, rep(toumei, length(sample_group_color_unique)), toumei, toumei, sp_color_unique)
    legend_pch = c(1, rep(21, length(sample_group_color_unique)), 1, 1, rep(1, length(sp_color_unique)))
    legend_font = c(2, rep(1, length(sample_group_color_unique)), 1, 2, rep(3, length(sp_color_unique)))
    plot.new()
    legend("right", legend = legend_text, pt.cex = 1, pch = legend_pch, lty = 0, lwd = 2, pt.bg = legend_bg, col = legend_fg, cex = cex_axis, text.font = legend_font)
}

prepare_metadata_table = function(dir_csca_input_table, selected_sample_groups, spp) {
    files = list.files(dir_csca_input_table, pattern = ".*metadata.*")
    df_metadata = data.frame()
    for (file in files) {
        metadata_path = file.path(dir_csca_input_table, file)
        tmp_metadata = read.table(metadata_path, header = TRUE, sep = '\t', quote = '', comment.char = '', check.names = FALSE)
        df_metadata = rbind(df_metadata, tmp_metadata)
    }
    df_metadata = df_metadata[(df_metadata[['sample_group']] %in% selected_sample_groups) & (df_metadata[['scientific_name']] %in% spp),]
    df_metadata = df_metadata[, !startsWith(colnames(df_metadata), 'Unnamed')]
    return(df_metadata)
}

get_label_orders = function(df_metadata) {
    order_cg = order(df_metadata[['sample_group']])
    label_orders = unique(paste(df_metadata[order_cg, 'scientific_name'], df_metadata[order_cg, 'sample_group'], sep = '_'))
    label_orders = sub(' ', '_', label_orders)
    return(label_orders)
}

extract_ortholog_mean_expression_table = function(df_singleog, averaged_tcs, label_orders) {
    averaged_orthologs = list()
    rowname_order = rownames(df_singleog)
    for (d in c('uncorrected', 'corrected')) {
        averaged_orthologs[[d]] = df_singleog
        averaged_orthologs[[d]][, 'busco_id'] = rownames(averaged_orthologs[[d]])
        for (sp_filled in colnames(df_singleog)) {
            tc = averaged_tcs[[d]][[sp_filled]]
            averaged_orthologs[[d]] = merge(averaged_orthologs[[d]], tc, by.x = sp_filled, by.y = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE)
        }
        num_remove_col = ncol(df_singleog) + 1
        rownames(averaged_orthologs[[d]]) = averaged_orthologs[[d]][['busco_id']]
        averaged_orthologs[[d]] = averaged_orthologs[[d]][rowname_order,]
        averaged_orthologs[[d]] = averaged_orthologs[[d]][, -(1:num_remove_col)]
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
    all_files = list.files(dir_csca_input_table, pattern = "*.tc.tsv")
    uncorrected_files = all_files[grepl("uncorrected", all_files)]
    corrected_files = all_files[((!grepl("uncorrected", all_files)) & (grepl(batch_effect_alg, all_files)))]
    for (sp in spp_filled) {
        uncorrected_file = uncorrected_files[startsWith(uncorrected_files, sp)]
        uncorrected_path = file.path(dir_csca_input_table, uncorrected_file)
        corrected_file = corrected_files[startsWith(corrected_files, sp)]
        corrected_path = file.path(dir_csca_input_table, corrected_file)
        if ((length(uncorrected_path) == 0) | (length(corrected_path) == 0)) {
            cat(paste("Skipping. `amalgkit curate` output(s) not found:", sp, "\n"), file = stderr())
            next
        }
        unaveraged_tcs[['uncorrected']][[sp]] = read.delim(uncorrected_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
        unaveraged_tcs[['corrected']][[sp]] = read.delim(corrected_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
    }
    return(unaveraged_tcs)
}

extract_ortholog_unaveraged_expression_table = function(df_singleog, unaveraged_tcs) {
    unaveraged_orthologs = list()
    rowname_order = rownames(df_singleog)
    for (d in c('uncorrected', 'corrected')) {
        unaveraged_orthologs[[d]] = df_singleog
        unaveraged_orthologs[[d]][, 'busco_id'] = rownames(unaveraged_orthologs[[d]])
        for (sp_filled in colnames(df_singleog)) {
            tc = unaveraged_tcs[[d]][[sp_filled]]
            colnames(tc) = paste(sp_filled, colnames(tc), sep = '_')
            unaveraged_orthologs[[d]] = merge(unaveraged_orthologs[[d]], tc, by.x = sp_filled, by.y = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE)
        }
        num_remove_col = length(spp) + 1
        rownames(unaveraged_orthologs[[d]]) = unaveraged_orthologs[[d]][['busco_id']]
        unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][rowname_order,]
        unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][, -(1:num_remove_col)]
        tc_order = order(sub('.*_', '', colnames(unaveraged_orthologs[[d]])))
        unaveraged_orthologs[[d]] = unaveraged_orthologs[[d]][, tc_order]
    }
    return(unaveraged_orthologs)
}

get_df_labels_averaged = function(df_metadata, label_orders, selected_sample_groups, sample_group_colors) {
    metadata_tmp = df_metadata[(df_metadata[['exclusion']] == 'no'),]
    df_label = unique(metadata_tmp[, c('scientific_name', 'sample_group')])
    categories = list(scientific_name = metadata_tmp[['scientific_name']], sample_group = metadata_tmp[['sample_group']])
    df_bp = aggregate(metadata_tmp[['bioproject']], by = categories, function(x) { length(unique(x)) })
    colnames(df_bp) = c('scientific_name', 'sample_group', 'num_bp')
    df_label = merge(df_label, df_bp, all.x = TRUE, all.y = FALSE)
    df_run = aggregate(metadata_tmp[['run']], by = categories, function(x) { length(unique(x)) })
    colnames(df_run) = c('scientific_name', 'sample_group', 'num_run')
    df_label = merge(df_label, df_run, all.x = TRUE, all.y = FALSE)
    df_label = df_label[order(df_label[['sample_group']], df_label[['scientific_name']]),]
    df_label = sort_labels(df_label, label_orders)
    df_label = add_color_to_metadata(df_label, selected_sample_groups, sample_group_colors)
    df_label = sort_labels(df_label, label_orders)
    rownames(df_label) = NULL
    write.table(df_label, paste0('csca_color_averaged.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
    return(df_label)
}

get_df_labels_unaveraged = function(df_metadata, selected_sample_groups, sample_group_colors) {
    cols = c('run', 'bioproject', 'sample_group', 'scientific_name', 'sp_color', 'sample_group_color', 'bp_color')
    metadata_tmp = df_metadata[(df_metadata[['exclusion']] == 'no'),]
    df_color = add_color_to_metadata(metadata_tmp, selected_sample_groups, sample_group_colors)
    df_color = df_color[, cols]
    label_order = order(df_color[['run']])
    df_color = df_color[label_order,]
    write.table(df_color, paste0('csca_color_unaveraged.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
    return(df_color)
}

add_color_to_metadata = function(df, selected_sample_groups, sample_group_colors) {
    df = df[, (!colnames(df) %in% c('bp_color', 'sp_color', 'sample_group_color'))]
    scientific_name = as.character(df[['scientific_name']])
    sample_group = as.character(df[['sample_group']])
    scientific_name_unique = sort(scientific_name[!duplicated(scientific_name)])
    if (length(sample_group_colors) == 1 && sample_group_colors == 'DEFAULT') {
        if (length(selected_sample_groups) <= 8) {
            sample_group_color = brewer.pal(length(unique(sample_group)), "Dark2")
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
        } else if (length(selected_sample_groups) <= 12) {
            sample_group_color = brewer.pal(length(unique(sample_group)), "Paired")
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
        } else {
            sample_group_color = rainbow_hcl(length(selected_sample_groups), c = 100)
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 150)
        }
    } else {
        if (length(sample_group_colors) != length(selected_sample_groups)) {
            stop("Length of sample_group_colors must match length of selected_sample_groups")
        }
        sample_group_color = sample_group_colors
        sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
    }
    sample_group_unique = selected_sample_groups
    df_sample_group = data.frame(sample_group = sample_group_unique, sample_group_color = sample_group_color[1:length(sample_group_unique)], stringsAsFactors = FALSE)
    df_sp = data.frame(scientific_name = scientific_name_unique, sp_color = sp_color[1:length(scientific_name_unique)], stringsAsFactors = FALSE)
    df = merge(df, df_sp, sort = FALSE, all.y = FALSE)
    df = merge(df, df_sample_group, sort = FALSE, all.y = FALSE)
    if ('bioproject' %in% colnames(df)) {
        bioproject = as.character(df$bioproject)
        bp_color = rainbow_hcl(length(unique(bioproject)), c = 50)
        df_bp = data.frame(bioproject = unique(bioproject), bp_color = bp_color[1:length(unique(bioproject))], stringsAsFactors = FALSE)
        df = merge(df, df_bp, sort = FALSE, all.y = FALSE)
    }
    return(df)
}

save_averaged_tsne_plot = function(tc, df_label) {
    cat('Generating averaged t-SNE plot.\n')
    num_samples = ncol(tc)
    perplexity = min(30, floor(num_samples / 4), na.rm = TRUE)
    if ((num_samples < 4) || (!is.finite(perplexity)) || (perplexity <= 0)) {
        cat('Skipping averaged t-SNE: insufficient number of samples.\n')
        return()
    }
    set.seed(1)
    out_tsne = try(Rtsne(as.matrix(t(tc)), theta = 0, check_duplicates = FALSE, verbose = FALSE, perplexity = perplexity, dims = 2))
    if ("try-error" %in% class(out_tsne)) {
        flag_tsne_success = FALSE
        cat("t-SNE failed and will be skipped.\n")
    } else {
        flag_tsne_success = TRUE
    }
    if (!flag_tsne_success) {
        return()
    }
    try_out = tryCatch(
    {
        plot(out_tsne$Y[, 1], out_tsne$Y[, 2], pch = 21, cex = 2, lwd = 1, bg = df_label[['sample_group_color']],
             col = df_label[['sp_color']], xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", las = 1)
    },
        error = function(a) { return("t-SNE plot failed.") }
    )
    if (mode(try_out) == "character") {
        cat('t-SNE failed.\n')
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    }
}

get_pca_coordinates = function(tc, df_label, by = 'species_sample_group') {
    tc_dist_matrix = cor(tc, method = 'pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    #set.seed(1)
    pca = prcomp(tc_dist_matrix)
    pca_summary = summary(pca)
    num_pc_available = ncol(pca[['x']])
    labels = c()
    for (i in 1:5) {
        if (i <= num_pc_available) {
            labels = c(labels, paste0("Principal component ", i, " (", round(pca_summary$importance[2, i] * 100, digits = 1), "%)"))
        } else {
            labels = c(labels, paste0("Principal component ", i, " (NA%)"))
        }
    }
    pca_x = pca[['x']]
    empty_pc = rep(NA_real_, nrow(pca_x))
    PC1 = if ('PC1' %in% colnames(pca_x)) pca_x[, 'PC1'] else empty_pc
    PC2 = if ('PC2' %in% colnames(pca_x)) pca_x[, 'PC2'] else empty_pc
    PC3 = if ('PC3' %in% colnames(pca_x)) pca_x[, 'PC3'] else empty_pc
    PC4 = if ('PC4' %in% colnames(pca_x)) pca_x[, 'PC4'] else empty_pc
    PC5 = if ('PC5' %in% colnames(pca_x)) pca_x[, 'PC5'] else empty_pc
    tmp = data.frame(PC1, PC2, PC3, PC4, PC5)
    if (by == 'species_sample_group') {
        df_label[by] = paste0(sub(' ', '_', df_label[['scientific_name']]), '_', df_label[['sample_group']])
        tmp[by] = rownames(tmp)
    } else if (by == 'run') {
        tmp[by] = sub('.*_', '', rownames(tmp))
    } else {
        tmp[by] = rownames(tmp)
    }
    tmp = merge(df_label, tmp, by = by)
    return(list(tmp, labels))
}

save_unaveraged_pca_plot = function(unaveraged_orthologs, df_color_unaveraged, df_metadata) {
    cat('Generating unaveraged PCA plot.\n')
    for (d in c('uncorrected', 'corrected')) {
        out = get_pca_coordinates(tc = unaveraged_orthologs[[d]], df_label = df_color_unaveraged, by = 'run')
        tmp = out[[1]]
        pc_contributions = out[[2]]
        pc_cols = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')
        pc_cols2 = paste(pc_cols, d, sep = '_')
        sorted_cols = c(colnames(df_metadata), pc_cols2)
        tmp2 = tmp[, c('run', pc_cols)]
        colnames(tmp2) = c('run', pc_cols2)
        df_metadata = merge(df_metadata, tmp2, all.x = TRUE, by = 'run', sort = FALSE)
        df_metadata = df_metadata[, sorted_cols]

        df_color_uniq = unique(df_color_unaveraged[, c('sample_group', 'sample_group_color')])
        sample_group_colors = df_color_uniq[['sample_group_color']]
        names(sample_group_colors) = df_color_uniq[['sample_group']]

        for (pcxy in list(c(1, 2), c(3, 4))) {
            pcx = pcxy[1]
            pcy = pcxy[2]
            colx = paste0('PC', pcx)
            coly = paste0('PC', pcy)

            # Check if data is finite
            xvals = tmp[[colx]]
            yvals = tmp[[coly]]
            finite_x = xvals[is.finite(xvals)]
            finite_y = yvals[is.finite(yvals)]

            if (length(finite_x) == 0 || length(finite_y) == 0) {
                message(sprintf("Skipping PC%d vs PC%d plot for '%s' - no finite data.", pcx, pcy, d))
                next
            }

            xmin = min(finite_x)
            xmax = max(finite_x)
            ymin = min(finite_y)
            ymax = max(finite_y)

            # Check if there's any range to plot
            if (xmin == xmax || ymin == ymax) {
                message(sprintf("Skipping PC%d vs PC%d plot for '%s' - no variation in data.", pcx, pcy, d))
                next
            }

            # Add a small buffer around the min/max
            xunit = (xmax - xmin) * 0.01
            yunit = (ymax - ymin) * 0.01
            xmin = xmin - xunit
            xmax = xmax + xunit
            ymin = ymin - yunit
            ymax = ymax + yunit

            g = ggplot(tmp, aes(x = !!rlang::sym(colx), y = !!rlang::sym(coly), color = sample_group)) +
                theme_bw() +
                geom_point(size = 0.5, alpha = 0.3) +
                # Add density lines only if there is more than one unique point
                # This reduces the chance of zero contour issues.
                geom_density_2d(mapping = aes(color = sample_group), bins = 12, linewidth = 0.25) +
                scale_color_manual(values = sample_group_colors) +
                xlab(pc_contributions[pcx]) +
                ylab(pc_contributions[pcy]) +
                xlim(xmin, xmax) +
                ylim(ymin, ymax) +
                theme(
                    axis.text = element_text(size = font_size),
                    axis.title = element_text(size = font_size),
                    legend.text = element_text(size = font_size),
                    legend.title = element_text(size = font_size)
                )

            filename = paste0('csca_unaveraged_pca_PC', pcx, pcy, '.', d, '.pdf')
            tryCatch(
                {
                    ggsave(file = filename, g, height = 2.15, width = 4.25)
                },
                error = function(cond) {
                    message(paste("PCA could not be computed for file", filename))
                }
            )
        }
    }
    return(df_metadata)
}

get_tsne_coordinates = function(tc, df_label, by = 'run') {
    num_samples = ncol(tc)
    perplexity = min(30, floor(num_samples / 4), na.rm = TRUE)
    if ((num_samples < 4) || (!is.finite(perplexity)) || (perplexity <= 0)) {
        return(NULL)
    }
    set.seed(1)
    out_tsne = Rtsne(as.matrix(t(tc)), theta = 0, check_duplicates = FALSE, verbose = FALSE, perplexity = perplexity, dims = 2)
    tmp = data.frame(tsne1 = out_tsne$Y[, 1], tsne2 = out_tsne$Y[, 2])
    tmp[[by]] = sub('.*_', '', colnames(tc))
    tmp = merge(df_label, tmp, by = by)
    return(tmp)
}

save_unaveraged_tsne_plot = function(unaveraged_orthologs, df_color_unaveraged) {
    cat('Generating unaveraged t-SNE plot.\n')
    df_color_uniq = unique(df_color_unaveraged[, c('sample_group', 'sample_group_color')])
    sample_group_colors = df_color_uniq[['sample_group_color']]
    names(sample_group_colors) = df_color_uniq[['sample_group']]
    for (d in c('uncorrected', 'corrected')) {
        tmp = get_tsne_coordinates(tc = unaveraged_orthologs[[d]], df_label = df_color_unaveraged)
        if (is.null(tmp)) {
            cat(sprintf('Skipping unaveraged t-SNE (%s): insufficient number of samples.\n', d))
            next
        }
        pcx = 1
        pcy = 2
        colx = paste0('tsne', pcx)
        coly = paste0('tsne', pcy)
        xmin = min(tmp[[colx]], na.rm = TRUE)
        xmax = max(tmp[[colx]], na.rm = TRUE)
        xunit = (xmax - xmin) * 0.01
        xmin = xmin - xunit
        xmax = xmax + xunit

        ymin = min(tmp[[coly]], na.rm = TRUE)
        ymax = max(tmp[[coly]], na.rm = TRUE)
        yunit = (ymax - ymin) * 0.01
        ymin = ymin - yunit
        ymax = ymax + yunit

        g = ggplot(tmp, aes(x = !!rlang::sym(colx), !!rlang::sym(coly), color = sample_group))
        g = g + theme_bw()
        g = g + geom_point(size = 0.5)
        g = g + geom_density_2d(mapping = aes(color = sample_group), bins = 12, linewidth = 0.25)
        g = g + scale_color_manual(values = sample_group_colors)
        g = g + xlab('t-SNE dimension 1')
        g = g + ylab('t-SNE dimension 2')
        g = g + xlim(xmin, xmax)
        g = g + ylim(ymin, ymax)
        g = g + theme(
            axis.text = element_text(size = font_size),
            axis.title = element_text(size = font_size),
            legend.text = element_text(size = font_size),
            legend.title = element_text(size = font_size)
        )
        filename = paste0('csca_unaveraged_tsne.', d, '.pdf')
        tryCatch(
        {
            ggsave(file = filename, g, height = 2.15, width = 4.25)
        },
            error = function(cond) {
                message(paste("t-SNE could not be computed for file", filename))
                #message("original error message:")
                #message(conditionMessage(cond))
            }
        )
    }
}

save_averaged_heatmap_plot = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged heatmap.\n')
    file_name = 'csca_SVA_heatmap.pdf'
    pdf(file_name, height = 3.3, width = 7.2) # full figure size = 9.7 x 7.2
    layout_matrix = matrix(c(
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        1, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                           2, 10, byrow = TRUE)
    layout(t(layout_matrix))
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(0.27, 0.5, 'Uncorrected', srt = 0, cex = 2)
    text(0.80, 0.5, 'Corrected', srt = 0, cex = 2)
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        df_label = df_color_averaged
        par(mar = c(0, 0, 0, 0))
        draw_multisp_heatmap(tc = tc, df_label = df_label)
    }
    graphics.off()
}

save_averaged_dendrogram_plot = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged dendrogram.\n')
    file_name = 'csca_SVA_dendrogram.pdf'
    pdf(file_name, height = 2.5, width = 7.2) # full figure size = 9.7 x 7.2
    layout_matrix = matrix(c(1, 2), 2, 1, byrow = TRUE)
    layout(t(layout_matrix))
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        df_label = df_color_averaged
        par(cex = 0.5, mar = c(10, 5.5, 0, 0), mgp = c(4, 0.7, 0))
        pvclust_file = paste0('csca_pvclust_', d, '.RData')
        draw_multisp_dendrogram(tc = tc, df_label = df_label, df_metadata = df_metadata, pvclust_file = pvclust_file,
                                nboot = 1, cex.xlab = 0.3, cex.yaxis = 0.5)
    }
    graphics.off()
}

save_averaged_dimensionality_reduction_summary = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged dimensionality reduction summary.\n')
    par(cex = 1)
    file_name = 'csca_averaged_summary.pdf'
    pdf(file_name, height = 7.2, width = 7.2) # full figure size = 9.7 x 7.2
    layout_matrix = matrix(c(1, 1, 1, 4, 4, 4, 7, 7, 2, 2, 2, 5, 5, 5, 7, 7, 3, 3, 3, 6, 6, 6, 7, 7), 3, 8, byrow = TRUE)
    layout(layout_matrix)
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        df_label = df_color_averaged
        par(mar = c(4, 4, 0.1, 1)); draw_multisp_pca(tc = tc, df_label = df_label)
        par(mar = c(4, 4, 0.1, 1)); draw_multisp_tsne(tc = tc, df_label = df_label)
        par(mar = c(4, 4, 0.1, 1)); draw_multisp_mds(tc = tc, df_label = df_label)
    }
    par(mar = c(0, 0, 0, 0)); draw_multisp_legend(df_label)
    graphics.off()
}

draw_multisp_boxplot = function(df_metadata, tc_dist_matrix, fontsize = 8) {
    is_same_sp = outer(df_metadata[['scientific_name']], df_metadata[['scientific_name']], function(x, y) { x == y })
    is_same_sample_group = outer(df_metadata[['sample_group']], df_metadata[['sample_group']], function(x, y) { x == y })
    plot(c(0.5, 4.5), c(0, 1), type = 'n', xlab = '', ylab = "Pearson's correlation\ncoefficient", las = 1, xaxt = 'n')
    boxplot(tc_dist_matrix[(!is_same_sp) & (!is_same_sample_group)], at = 1, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(tc_dist_matrix[(is_same_sp) & (!is_same_sample_group)], at = 2, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(tc_dist_matrix[(!is_same_sp) & (is_same_sample_group)], at = 3, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(tc_dist_matrix[(is_same_sp) & (is_same_sample_group)], at = 4, add = TRUE, col = 'gray', yaxt = 'n')
    labels = c('bw\nbw', 'bw\nwi', 'wi\nbw', 'wi\nwi')
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = 'Sample group\nSpecies', padj = 0.5, hadj = 1, tick = FALSE)
}

save_averaged_box_plot = function(averaged_orthologs, df_color_averaged) {
    cat('Generating averaged boxplot.\n')
    file_name = 'csca_boxplot.pdf'
    pdf(file_name, height = 3.6, width = 7.2) # full figure size = 9.7 x 7.2
    par(mfrow = c(1, 2))
    for (d in c('uncorrected', 'corrected')) {
        tc = averaged_orthologs[[d]]
        tc[tc < 0] = 0
        tc_dist_matrix = cor(tc, method = 'pearson')
        tc_dist_matrix[is.na(tc_dist_matrix)] = 0
        draw_multisp_boxplot(df_color_averaged, tc_dist_matrix, fontsize = 8)
    }
    graphics.off()
}

calculate_correlation_within_group = function(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_sample_groups, dist_method = 'pearson') {
    my_fun1 = function(x) { median(x, na.rm = TRUE) }
    for (d in c('uncorrected', 'corrected')) {
        ortholog_med = data.frame(matrix(NA, nrow(averaged_orthologs[[d]]), length(selected_sample_groups)))
        colnames(ortholog_med) = selected_sample_groups
        rownames(ortholog_med) = rownames(averaged_orthologs[[d]])
        for (sample_group in selected_sample_groups) {
            is_sample_group = endsWith(colnames(averaged_orthologs[[d]]), sample_group)
            ortholog_med[, sample_group] = apply(as.data.frame(averaged_orthologs[[d]][, is_sample_group]), 1, my_fun1)
        }
        stopifnot(all(rownames(unaveraged_orthologs[[d]]) == rownames(ortholog_med)))
        target_col = paste0('within_group_cor_', d)
        nongroup_col = paste0('max_nongroup_cor_', d)
        df_metadata[, target_col] = NA
        df_metadata[, nongroup_col] = NA
        for (sp_and_run in colnames(unaveraged_orthologs[[d]])) {
            split_string = strsplit(sp_and_run, "_")[[1]]
            sp = paste(split_string[1:2], collapse = " ")
            sra_run = paste(split_string[3:length(split_string)], collapse = '_')
            is_sra = (df_metadata[['run']] == sra_run) & (df_metadata[['scientific_name']] == sp)
            if (sum(is_sra) == 0) {
                warning(paste('Sample skipped:', sp_and_run))
                next
            }
            sample_cg = df_metadata[is_sra, 'sample_group']
            sample_values = unaveraged_orthologs[[d]][, sp_and_run]
            for (sample_group in colnames(ortholog_med)) {
                med_values = ortholog_med[, sample_group]
                is_na = (is.na(sample_values) | is.na(med_values))
                sample_values2 = sample_values[!is_na]
                med_values2 = med_values[!is_na]
                cor_coef = cor(sample_values2, med_values2, method = dist_method)
                if (sample_cg == sample_group) {
                    df_metadata[is_sra, target_col] = cor_coef
                } else {
                    max_value = max(cor_coef, df_metadata[is_sra, nongroup_col], na.rm = TRUE)
                    df_metadata[is_sra, nongroup_col] = max_value
                }
            }
        }
    }
    return(df_metadata)
}

save_group_cor_histogram = function(df_metadata, font_size = 8) {
    cat('Generating unaveraged group correlation histogram.\n')

    # Columns we will iterate over
    cor_cols = c('within_group_cor_uncorrected', 'within_group_cor_corrected')
    fill_by_vars = c('sample_group', 'scientific_name')

    # Pre-filter the data to remove rows with NA/Non-finite values in these columns
    df_clean = df_metadata
    for (col in cor_cols) {
        df_clean = df_clean[is.finite(df_clean[[col]]), ]
    }

    # Limit values to avoid out-of-range warnings if needed, e.g. to [0, 1]
    # Adjust or remove this if not required.
    for (col in cor_cols) {
        df_clean = df_clean[!is.na(df_clean[[col]]) & df_clean[[col]] >= 0 & df_clean[[col]] <= 1, ]
    }

    max_count <- 0
    plot_list <- list()

    for (col in cor_cols) {
        for (fill_by in fill_by_vars) {
            tmp = df_clean[!is.na(df_clean[[col]]), ]
            g = ggplot2::ggplot(tmp) +
                geom_histogram(aes(x = !!rlang::sym(col), fill = !!rlang::sym(fill_by)),
                               position = "stack", alpha = 0.7, bins = 40, na.rm = TRUE) +
                theme_bw(base_size = font_size) +
                xlim(c(0, 1)) +
                labs(x = col, y = 'Count') +
                theme(
                    axis.text = element_text(size = font_size, color = 'black'),
                    axis.title = element_text(size = font_size, color = 'black'),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    rect = element_rect(fill = "transparent"),
                    plot.margin = unit(rep(0.1, 4), "cm")
                )
            plot_list[[paste0(col, "_", fill_by)]] <- g
        }
    }

    final_plot <- (plot_list[['within_group_cor_uncorrected_scientific_name']] +
        plot_list[['within_group_cor_corrected_scientific_name']]) /
        (plot_list[['within_group_cor_uncorrected_sample_group']] +
            plot_list[['within_group_cor_corrected_sample_group']])

    ggsave(filename = "csca_within_group_cor.pdf", plot = final_plot, width = 7.2, height = 6.0)
}

extract_selected_tc_only = function(unaveraged_tcs, df_metadata) {
    selected_runs = df_metadata[(df_metadata[['exclusion']] == 'no'), 'run']
    for (d in c('uncorrected', 'corrected')) {
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            is_selected = colnames(unaveraged_tcs[[d]][[sci_name]]) %in% selected_runs
            if (sum(is_selected) == 0) {
                warning(paste('No', d, 'samples were selected:', sci_name))
                next
            }
            unaveraged_tcs[[d]][[sci_name]] = unaveraged_tcs[[d]][[sci_name]][, is_selected]
        }
    }
    return(unaveraged_tcs)
}

unaveraged2averaged = function(unaveraged_tcs, df_metadata, selected_sample_groups) {
    is_sample_groups = list()
    for (sample_group in selected_sample_groups) {
        is_sample_groups[[sample_group]] = (df_metadata[['sample_group']] == sample_group)
    }
    is_sci_names = list()
    for (sci_name in unique(df_metadata[['scientific_name']])) {
        sci_name_ub = sub(' ', '_', sci_name)
        is_sci_names[[sci_name_ub]] = (df_metadata[['scientific_name']] == sci_name)
    }
    is_not_excluded = (df_metadata[['exclusion']] == 'no')
    averaged_tcs = list()
    for (d in c('uncorrected', 'corrected')) {
        averaged_tcs[[d]] = list()
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            n_rows = nrow(unaveraged_tcs[[d]][[sci_name]])
            averaged_tcs[[d]][[sci_name]] = data.frame(matrix(ncol = 0, nrow = n_rows))
            rownames(averaged_tcs[[d]][[sci_name]]) = rownames(unaveraged_tcs[[d]][[sci_name]])
            for (sample_group in selected_sample_groups) {
                is_target = (is_sample_groups[[sample_group]] &
                    is_sci_names[[sci_name]] &
                    is_not_excluded)
                target_runs = df_metadata[is_target, 'run']
                target_runs = target_runs[target_runs %in% colnames(unaveraged_tcs[[d]][[sci_name]])]
                if (length(target_runs) == 0) {
                    next
                }
                label = paste(sub(' ', '_', sci_name), sample_group, sep = '_')
                if (sum(is_target) == 1) {
                    averaged_tcs[[d]][[sci_name]][, label] = unaveraged_tcs[[d]][[sci_name]][, target_runs]
                } else {
                    averaged_tcs[[d]][[sci_name]][, label] = apply(unaveraged_tcs[[d]][[sci_name]][, target_runs], 1, mean)
                }
            }
            if (ncol(averaged_tcs[[d]][[sci_name]]) == 0) {
                averaged_tcs[[d]][[sci_name]] = NULL
            }
        }
    }
    return(averaged_tcs)
}

save_group_cor_scatter = function(df_metadata, font_size = 8) {
    cat('Generating unaveraged group correlation scatter plot.\n')
    alpha_value = 0.2
    improvement_xymin = 0.5
    improvement_xymax = 2.0

    # Compute new columns
    df_metadata[['corrected_per_uncorrected_group_cor']] = df_metadata[['within_group_cor_corrected']] / df_metadata[['within_group_cor_uncorrected']]
    df_metadata[['corrected_per_uncorrected_max_nongroup_cor']] = df_metadata[['max_nongroup_cor_corrected']] / df_metadata[['max_nongroup_cor_uncorrected']]

    # Limit values to the specified range
    df_metadata[['corrected_per_uncorrected_group_cor']] = pmax(pmin(df_metadata[['corrected_per_uncorrected_group_cor']], improvement_xymax), improvement_xymin)
    df_metadata[['corrected_per_uncorrected_max_nongroup_cor']] = pmax(pmin(df_metadata[['corrected_per_uncorrected_max_nongroup_cor']], improvement_xymax), improvement_xymin)

    # Remove rows with NA, NaN, or Inf
    df_clean <- df_metadata[is.finite(df_metadata$within_group_cor_uncorrected) &
                            is.finite(df_metadata$within_group_cor_corrected) &
                            is.finite(df_metadata$max_nongroup_cor_uncorrected) &
                            is.finite(df_metadata$max_nongroup_cor_corrected) &
                            is.finite(df_metadata$corrected_per_uncorrected_group_cor) &
                            is.finite(df_metadata$corrected_per_uncorrected_max_nongroup_cor), ]

    # Plotting
    ps = list()
    ps[[1]] = ggplot(df_clean, aes(x = max_nongroup_cor_uncorrected, y = within_group_cor_uncorrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[2]] = ggplot(df_clean, aes(x = max_nongroup_cor_corrected, y = within_group_cor_corrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[3]] = ggplot(df_clean, aes(x = within_group_cor_uncorrected, y = within_group_cor_corrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[4]] = ggplot(df_clean, aes(x = max_nongroup_cor_uncorrected, y = max_nongroup_cor_corrected)) +
        xlim(c(0, 1)) + ylim(c(0, 1))
    ps[[5]] = ggplot(df_clean, aes(x = corrected_per_uncorrected_max_nongroup_cor, y = corrected_per_uncorrected_group_cor)) +
        xlim(c(improvement_xymin, improvement_xymax)) + ylim(c(improvement_xymin, improvement_xymax))

    base_layers = list(
        geom_point(alpha = alpha_value, na.rm = TRUE),
        geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'blue'),
        theme_bw(),
        theme(
            text = element_text(size = font_size),
            axis.text = element_text(size = font_size),
            axis.title = element_text(size = font_size),
            legend.text = element_text(size = font_size),
            legend.title = element_text(size = font_size)
        )
    )
    density_layer = stat_density_2d(bins = 12, linewidth = 0.25, color = 'gray', na.rm = TRUE)

    # Try with density contours; fall back to scatter-only if density estimation fails
    for (i in seq_along(ps)) {
        ps[[i]] = Reduce('+', base_layers, ps[[i]])
    }
    ps_with_density = lapply(ps, function(p) p + density_layer)
    ps_with_density[[6]] = ggplot() + theme_void()
    ps[[6]] = ggplot() + theme_void()

    combined_plot = wrap_plots(ps_with_density)
    tryCatch({
        ggsave(filename = "csca_group_cor_scatter.pdf", plot = combined_plot, width = 7.2, height = 4.8)
    }, error = function(e) {
        cat("Density estimation failed (insufficient data points), plotting scatter only.\n")
        combined_plot = wrap_plots(ps)
        ggsave(filename = "csca_group_cor_scatter.pdf", plot = combined_plot, width = 7.2, height = 4.8)
    })
}

write_pivot_table = function(df_metadata, unaveraged_tcs, selected_sample_groups) {
    d = 'corrected'
    is_selected_cg = (df_metadata[['sample_group']] %in% selected_sample_groups)
    is_loaded_run = FALSE
    for (sci_name_ub in names(unaveraged_tcs[[d]])) {
        sci_name = sub('_', ' ', sci_name_ub)
        is_loaded_run_sp = ((df_metadata[['run']] %in% colnames(unaveraged_tcs[[d]][[sci_name_ub]])) & (df_metadata[['scientific_name']] == sci_name))
        is_loaded_run = (is_loaded_run | is_loaded_run_sp)
    }
    is_selected = (is_selected_cg & is_loaded_run)
    tmp = df_metadata[is_selected, c('scientific_name', 'sample_group')]
    pivot_table = as.data.frame.matrix(table(tmp[['scientific_name']], tmp[['sample_group']]))
    pivot_table = cbind(data.frame(scientific_name = rownames(pivot_table)), pivot_table)
    write.table(pivot_table, 'csca_pivot_selected_samples.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
}

save_delta_pcc_plot = function(directory, plot_title) {

    first_row_means <- list()
    last_row_means <- list()
    all_data_points <- list()

    file_list <- list.files(directory, pattern = "\\.correlation_statistics\\.tsv$")
    for (filename in file_list) {
        filepath <- file.path(directory, filename)

        df <- read.table(filepath, sep = '\t', header = TRUE)

        first_row_means[[filename]] <- df[1, c('bwbw_mean', 'bwwi_mean', 'wiwi_mean', 'wibw_mean')]
        last_row_means[[filename]] <- df[nrow(df), c('bwbw_mean', 'bwwi_mean', 'wiwi_mean', 'wibw_mean')]

        all_data_points[[filename]] <- as.numeric(unlist(df[c('bwbw_mean', 'bwwi_mean', 'wiwi_mean', 'wibw_mean')]))
    }

    # Concatenate the mean values for the first and last rows into DataFrames
    first_row_means_df <- do.call(rbind, first_row_means)
    colnames(first_row_means_df) <- c('bwbw_uncorrected', 'bwwi_uncorrected', 'wiwi_uncorrected', 'wibw_uncorrected')

    last_row_means_df <- do.call(rbind, last_row_means)
    colnames(last_row_means_df) <- c('bwbw_corrected', 'bwwi_corrected', 'wiwi_corrected', 'wibw_corrected')

    combined_means_df <- cbind(first_row_means_df, last_row_means_df)

    #cat("Combined Means DataFrame:\n")
    #print(combined_means_df)

    # Calculate deltas
    combined_means_df$delta_bwbw_wibw_uncorrected <- combined_means_df$bwbw_uncorrected - combined_means_df$wibw_uncorrected
    combined_means_df$delta_bwbw_wibw_corrected <- combined_means_df$bwbw_corrected - combined_means_df$wibw_corrected
    combined_means_df$delta_wiwi_bwwi_uncorrected <- combined_means_df$bwwi_uncorrected - combined_means_df$wiwi_uncorrected
    combined_means_df$delta_wiwi_bwwi_corrected <- combined_means_df$bwwi_corrected - combined_means_df$wiwi_corrected

    # Select only the delta columns and take the absolute values
    delta_means_df <- combined_means_df[, c('delta_bwbw_wibw_uncorrected', 'delta_bwbw_wibw_corrected',
                                            'delta_wiwi_bwwi_uncorrected', 'delta_wiwi_bwwi_corrected')]
    delta_means_df <- abs(delta_means_df)
    delta_means_df <- na.omit(delta_means_df)

    #cat("\nDelta Means DataFrame:\n")
    #print(delta_means_df)

    if (nrow(delta_means_df) > 2) {
        cat("\nShapiro-Wilk Test Results:\n")
        cat("delta_bwbw_wibw_uncorrected: ", shapiro.test(delta_means_df$delta_bwbw_wibw_uncorrected)$p.value, "\n")
        cat("delta_bwbw_wibw_corrected: ", shapiro.test(delta_means_df$delta_bwbw_wibw_corrected)$p.value, "\n")
        cat("delta_wiwi_bwwi_uncorrected: ", shapiro.test(delta_means_df$delta_wiwi_bwwi_uncorrected)$p.value, "\n")
        cat("delta_wiwi_bwwi_corrected: ", shapiro.test(delta_means_df$delta_wiwi_bwwi_corrected)$p.value, "\n")

        cat("\nDependent T-Test Results:\n")
        t_test_result_bw = t.test(delta_means_df$delta_bwbw_wibw_uncorrected, delta_means_df$delta_bwbw_wibw_corrected, paired = TRUE)$p.value
        t_test_result_wi = t.test(delta_means_df$delta_wiwi_bwwi_uncorrected, delta_means_df$delta_wiwi_bwwi_corrected, paired = TRUE)$p.value
        cat("delta_bwbw_wibw_uncorrected vs delta_bwbw_wibw_corrected: ", t_test_result_bw, "\n")
        cat("delta_wiwi_bwwi_uncorrected vs delta_wiwi_bwwi_corrected: ", t_test_result_wi, "\n")

        p_label1 <- paste("p =", signif(t_test_result_bw, 3))
        p_label2 <- paste("p =", signif(t_test_result_wi, 3))
    } else {
        cat("Not enough data points to perform statistical tests. P values will not be shown.\n")
    }


    pdf(plot_title, height = 4.5, width = 4.5, fonts = "Helvetica", pointsize = 7)

    plot(c(0.5, 4.5), c(0, 0.45), type = 'n', xlab = '', ylab = expression(Delta ~ "mean PCC"), las = 1, xaxt = 'n')
    boxplot(delta_means_df$delta_bwbw_wibw_uncorrected, at = 1, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(delta_means_df$delta_bwbw_wibw_corrected, at = 2, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(delta_means_df$delta_wiwi_bwwi_uncorrected, at = 3, add = TRUE, col = 'gray', yaxt = 'n')
    boxplot(delta_means_df$delta_wiwi_bwwi_corrected, at = 4, add = TRUE, col = 'gray', yaxt = 'n')

    if (nrow(delta_means_df) > 2) {
        segments(x0 = 1, y0 = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.2, x1 = 2, y1 = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.2, col = "black", lwd = 0.5)
        segments(x0 = 3, y0 = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.2, x1 = 4, y1 = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.2, col = "black", lwd = 0.5)
        text(x = 1.5, y = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.22, labels = p_label1, xpd = TRUE, srt = 0)
        text(x = 3.5, y = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.22, labels = p_label2, xpd = TRUE, srt = 0)
    }

    axis(side = 1, at = c(1, 2, 3, 4), labels = c("uncorr.\n", "corr.\n", "uncorr.\n", "corr.\n"), padj = 0.5, tick = FALSE)
    axis(side = 1, at = 0.35, labels = 'Correction\nSample group', padj = 0.5, hadj = 1, tick = FALSE)
    axis(side = 1, at = c(1.5, 3.5), labels = c("\nbetween sample group", "\nwithin sample group"), padj = 0.5, tick = FALSE)
    graphics.off()
}

save_sample_number_heatmap <- function(df_metadata, font_size = 8, dpi = 300) {
    sampled_data <- df_metadata[df_metadata$exclusion == 'no', c('scientific_name', 'sample_group')]
    freq_table <- table(sampled_data$scientific_name, sampled_data$sample_group)
    df_sample_count <- as.data.frame(freq_table)
    names(df_sample_count) <- c('scientific_name', 'sample_group', 'Freq')
    all_scientific_names <- unique(df_metadata$scientific_name)
    all_sample_groups <- unique(df_metadata$sample_group)
    complete_grid <- expand.grid(scientific_name = all_scientific_names, sample_group = all_sample_groups, stringsAsFactors = FALSE)
    df_sample_count <- merge(complete_grid, df_sample_count, by = c("scientific_name", "sample_group"), all.x = TRUE)
    df_sample_count$Freq[is.na(df_sample_count$Freq)] <- 0
    df_sample_count$log2_Freq <- ifelse(df_sample_count$Freq > 0, log2(df_sample_count$Freq + 1), 0)
    fill_max <- max(df_sample_count$log2_Freq)
    freq_breaks <- c(0, 1, 2, 4, 8, 16, 32, 64, 128)
    freq_breaks <- freq_breaks[freq_breaks <= max(df_sample_count$Freq)]
    if (!0 %in% freq_breaks) { freq_breaks <- c(0, freq_breaks) }
    if (!1 %in% freq_breaks) { freq_breaks <- c(freq_breaks, 1) }
    freq_breaks <- sort(unique(freq_breaks))
    log2_breaks <- log2(freq_breaks + 1)
    log2_breaks <- log2_breaks[log2_breaks <= fill_max]
    n_sample_groups <- length(all_sample_groups)
    n_scientific_names <- length(all_scientific_names)
    base_width <- 5
    base_height <- 5
    per_sample_group_width <- 0.1
    per_scientific_name_height <- 0.1
    width <- base_width + (per_sample_group_width * n_sample_groups)
    height <- base_height + (per_scientific_name_height * n_scientific_names)
    p <- ggplot(data = df_sample_count, aes(x = sample_group, y = scientific_name, fill = log2_Freq)) +
        geom_tile(color = "grey80") +
        geom_text(aes(label = Freq), size = 3, color = "black") +
        scale_fill_gradientn(
            colors = c("white", "#440154", "#482878", "#3E4989", "#31688E", "#35B779", "#4DCD63", "#76D730", "#B8DE29", "#FDE725"),
            values = sort(unique(c(0, log2_breaks / fill_max, 1))),
            limits = c(0, fill_max),
            name = "# of samples",
            breaks = log2_breaks,
            labels = freq_breaks
        ) +
        xlab('') +
        ylab('') +
        scale_y_discrete(
            limits = rev(unique(df_sample_count$scientific_name))
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_text(size = font_size),
            axis.text.x = element_text(size = font_size, angle = 90, hjust = 1),
            legend.title = element_text(size = font_size),
            legend.text = element_text(size = font_size)
        )
    ggsave('csca_sample_number_heatmap.pdf', plot = p, width = width, height = height, dpi = dpi)
}

df_og = read.table(file_orthogroup, header = TRUE, sep = '\t', row.names = 1, quote = '', check.names = FALSE)
df_gc = read.table(file_genecount, header = TRUE, sep = '\t', quote = '', check.names = FALSE)
rownames(df_gc) = df_gc[['orthogroup_id']]
df_gc[, 'orthogroup_id'] = NULL
spp_filled = colnames(df_gc)

is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
df_singleog = df_og[is_singlecopy, spp_filled, drop = FALSE]
spp = sub('_', ' ', spp_filled)
df_metadata = prepare_metadata_table(dir_csca_input_table, selected_sample_groups, spp)
label_orders = get_label_orders(df_metadata)
df_color_averaged = get_df_labels_averaged(df_metadata, label_orders, selected_sample_groups, sample_group_colors)
df_color_unaveraged = get_df_labels_unaveraged(df_metadata, selected_sample_groups, sample_group_colors)
cat('Number of orthologs in input table:', nrow(df_og), '\n')
cat('Number of selected single-copy orthologs:', nrow(df_singleog), '\n')
cat('Number of selected species:', length(spp), '\n')
save_sample_number_heatmap(df_metadata, font_size = font_size, dpi = 300)

unaveraged_tcs = load_unaveraged_expression_tables(dir_csca_input_table, spp_filled, batch_effect_alg)
unaveraged_tcs = extract_selected_tc_only(unaveraged_tcs, df_metadata)

# if a species was skipped during load_unaveraged_expression_tables(), it will cause indexing issues down the line, due to df_singleog having more species entries than the tcs
if (length(unaveraged_tcs[['corrected']]) < length(colnames(df_singleog))) {
    df_singleog = df_singleog[, names(unaveraged_tcs[['corrected']]), drop = FALSE]
}

unaveraged_orthologs = extract_ortholog_unaveraged_expression_table(df_singleog, unaveraged_tcs)
averaged_tcs = unaveraged2averaged(unaveraged_tcs, df_metadata, selected_sample_groups)
averaged_orthologs = extract_ortholog_mean_expression_table(df_singleog, averaged_tcs, label_orders)
write_pivot_table(df_metadata, unaveraged_tcs, selected_sample_groups)

cat('Applying expression level imputation for missing orthologs.\n')
imputed_averaged_orthologs = list()
imputed_unaveraged_orthologs = list()
for (d in c('uncorrected', 'corrected')) {
    imputed_averaged_orthologs[[d]] = impute_expression(averaged_orthologs[[d]])
    imputed_unaveraged_orthologs[[d]] = impute_expression(unaveraged_orthologs[[d]])
    write_table_with_index_name(
        df = averaged_orthologs[[d]],
        file_path = paste0('csca_ortholog_averaged.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
    write_table_with_index_name(
        df = unaveraged_orthologs[[d]],
        file_path = paste0('csca_ortholog_unaveraged.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
    write_table_with_index_name(
        df = imputed_averaged_orthologs[[d]],
        file_path = paste0('csca_ortholog_averaged.imputed.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
    write_table_with_index_name(
        df = imputed_unaveraged_orthologs[[d]],
        file_path = paste0('csca_ortholog_unaveraged.imputed.', d, '.tsv'),
        index_name = 'target_id',
        sort = FALSE
    )
}
cat(nrow(imputed_unaveraged_orthologs[[d]]), 'orthologs were found after filtering and imputation.\n')
df_metadata = calculate_correlation_within_group(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_sample_groups)
tryCatch(save_group_cor_scatter(df_metadata, font_size = 8), error = function(e) cat("Warning: save_group_cor_scatter skipped:", conditionMessage(e), "\n"))
tryCatch(save_group_cor_histogram(df_metadata, font_size = 8), error = function(e) cat("Warning: save_group_cor_histogram skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_tsne_plot(tc = imputed_unaveraged_orthologs[['corrected']], df_label = df_color_unaveraged), error = function(e) cat("Warning: save_averaged_tsne_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_heatmap_plot(imputed_averaged_orthologs, df_color_averaged), error = function(e) cat("Warning: save_averaged_heatmap_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_dendrogram_plot(imputed_averaged_orthologs, df_color_averaged), error = function(e) cat("Warning: save_averaged_dendrogram_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_dimensionality_reduction_summary(imputed_averaged_orthologs, df_color_averaged), error = function(e) cat("Warning: save_averaged_dimensionality_reduction_summary skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_box_plot(imputed_averaged_orthologs, df_color_averaged), error = function(e) cat("Warning: save_averaged_box_plot skipped:", conditionMessage(e), "\n"))
tryCatch({df_metadata = save_unaveraged_pca_plot(imputed_unaveraged_orthologs, df_color_unaveraged, df_metadata)}, error = function(e) cat("Warning: save_unaveraged_pca_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_unaveraged_tsne_plot(imputed_unaveraged_orthologs, df_color_unaveraged), error = function(e) cat("Warning: save_unaveraged_tsne_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_delta_pcc_plot(directory = dir_csca_input_table, plot_title = 'csca_delta_pcc_boxplot.pdf'), error = function(e) cat("Warning: save_delta_pcc_plot skipped:", conditionMessage(e), "\n"))

file_metadata_out = file.path(dir_csca, 'metadata.tsv')
write.table(df_metadata, file_metadata_out, row.names = FALSE, sep = '\t', quote = FALSE)

cat(sprintf('Number of SRA samples for exclusion potting: %s\n', formatC(nrow(df_metadata), format = 'd', big.mark = ',')))
out_path = file.path(dir_csca, 'csca_exclusion.pdf')
save_exclusion_plot(df = df_metadata, out_path = out_path, font_size = 8)

if (file.exists('Rplots.pdf')) {
    file.remove('Rplots.pdf')
}
cat('csca.r completed!\n')
