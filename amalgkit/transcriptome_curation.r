#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(pcaMethods, quietly = TRUE))
suppressPackageStartupMessages(library(colorspace, quietly = TRUE))
suppressPackageStartupMessages(library(RColorBrewer, quietly = TRUE))
suppressPackageStartupMessages(library(sva, quietly = TRUE))
suppressPackageStartupMessages(library(MASS, quietly = TRUE))
suppressPackageStartupMessages(library(NMF, quietly = TRUE))
suppressPackageStartupMessages(library(dendextend, quietly = TRUE))
suppressPackageStartupMessages(library(amap, quietly = TRUE))
suppressPackageStartupMessages(library(pvclust, quietly = TRUE))
suppressPackageStartupMessages(library(Rtsne, quietly = TRUE))

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
log_prefix = "transcriptome_curation.r:"
cat(log_prefix, "mode =", debug_mode, "\n")
if (debug_mode == "debug") {
    infile = '/Users/s229181/MSN/cstmm/cross_species_tmm_normalized_counts/Helianthus_annuus_cstmm_counts.tsv'
    eff_file = '/Users/s229181/MSN/eff_length/Helianthus_annuus_eff_length.tsv'
   # infile = "/Users/kf/Dropbox (Personal)/collaborators/Ken Naito/20210509_Vigna/gfe_data/cross_species_tmm_normalized_counts/Vigna_angularis_cstmm_counts.tsv"
    #eff_file = "/Users/kf/Dropbox (Personal)/collaborators/Ken Naito/20210509_Vigna/gfe_data/merge/Vigna_angularis_eff_length.tsv"
    dist_method = "pearson"
    mapping_rate_cutoff = 0
    min_dif = 0
    plot_intermediate = 0
    selected_tissues = c("root", "flower", "leaf")
    dir_count = "counts/"
    dir_eff_length = "eff_length/"
    mode = "msm"
    transform_method = "fpkm"
    # tmm norm debug
    correlation_threshold = 0.4
    tmm_norm = "yes"
    dir_work = '/Users/s229181/MSN/'
    srafile = '/Users/s229181/MSN/metadata/metadata_may2021_updated.tsv'
   # dir_work = "/Users/kf/Dropbox (Personal)/collaborators/Ken Naito/20210509_Vigna/gfe_data"
   # srafile = "/Users/kf/Dropbox (Personal)/collaborators/Ken Naito/20210509_Vigna/gfe_data/metadata/metadata/metadata_manual.tsv"

    dist_method = "pearson"
    min_dif = 0
    plot_intermediate = 1
    selected_tissues = c("root", "flower", "leaf")
    stop_after_tmm = TRUE

    # selected_tissues = strsplit('root|flower|leaf', '\\|')[[1]]
} else if (debug_mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    print(args)
    infile = args[1]
    srafile = args[2]
    dir_work = args[3]
    eff_file = args[4]
    dist_method = args[5]
    mapping_rate_cutoff = as.numeric(args[6])
    min_dif = as.numeric(args[7])
    plot_intermediate = as.integer(args[8])
    selected_tissues = strsplit(args[9], "\\|")[[1]]
    transform_method = args[10]
    one_outlier_per_iteration = as.integer(args[11])
    correlation_threshold = as.numeric(args[12])
    

}
if (!endsWith(dir_work, "/")) {
    dir_work = paste0(dir_work, "/")
}

# set directory
cat(log_prefix, "dir_work =", dir_work, "\n")
dir.create(dir_work, showWarnings = FALSE)

dir_curate = file.path(dir_work,'curate')
if (!file.exists(dir_curate)) {
  dir.create(dir_curate)
}
setwd(dir_curate)

# create output directories
dir_pdf = file.path(dir_curate,'plots')
if (!file.exists(dir_pdf)) {
  dir.create(dir_pdf)
}
dir_rdata = file.path(dir_curate,'rdata')
if (!file.exists(dir_rdata)) {
  dir.create(dir_rdata)
}
dir_tsv = file.path(dir_curate,'tables')
if (!file.exists(dir_tsv)) {
  dir.create(dir_tsv)
}

tc_sra_intersect = function(tc, sra) {
    sra_run = sra[['run']]
    tc = tc[, colnames(tc) %in% sra_run]
    sra = sra[sra[['run']] %in% colnames(tc), ]
    return(list(tc = tc, sra = sra))
}

remove_nonexpressed_gene = function(tc) {
    gene_sum = apply(tc, 1, sum)
    tc_ex = tc[gene_sum > 0, ]
    tc_ne = tc[gene_sum == 0, ]
    return(list(tc_ex = tc_ex, tc_ne = tc_ne))
}

add_color_to_sra = function(sra, selected_tissues) {
    sra = sra[, (!names(sra) %in% c("bp_color", "sp_color", "tissue_color"))]
    bioproject = as.character(sra[['bioproject']])
    bioproject_u = sort(unique(bioproject))
    scientific_name = as.character(sra[['scientific_name']])
    scientific_name_u = sort(unique(scientific_name))
    tissue = as.character(sra[['tissue']])
    tissue_u = sort(unique(tissue))
    if (length(selected_tissues) <= 8) {
        tissue_color = brewer.pal(length(tissue_u), "Dark2")
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
    } else if (length(selected_tissues) <= 12) {
        tissue_color = brewer.pal(length(unique(tissue)), "Paired")
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
    } else {
        tissue_color = rainbow_hcl(length(selected_tissues), c = 100)
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 150)
    }
    df_tissue = data.frame(tissue=tissue_u, tissue_color=tissue_color[1:length(tissue_u)], stringsAsFactors = FALSE)
    df_bp = data.frame(bioproject=bioproject_u, bp_color=bp_color[1:length(bioproject_u)], stringsAsFactors = FALSE)
    df_sp = data.frame(scientific_name=scientific_name_u, sp_color=sp_color[1:length(scientific_name_u)], stringsAsFactors = FALSE)
    sra = merge(sra, df_bp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_sp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_tissue, sort = FALSE, all.y = FALSE)
    return(sra)
}

sort_tc_and_sra = function(tc, sra, sort_columns = c("tissue", "scientific_name", "bioproject")) {
    for (column in rev(sort_columns)) {
        sra = sra[order(sra[column]), ]
    }
    sra_intersection = sra[(sra[['run']] %in% colnames(tc)),'run']
    tc = tc[, sra_intersection]
    return(list(tc = tc, sra = sra))
}

sort_averaged_tc = function(tc) {
    split_colnames = strsplit(colnames(tc), "_")
    genus_names = c()
    specific_names = c()
    tissue_names = c()
    for (i in 1:length(split_colnames)) {
        genus_names = c(genus_names, split_colnames[[i]][1])
        specific_names = c(specific_names, split_colnames[[i]][2])
        tissue_names = c(tissue_names, split_colnames[[i]][3])
    }
    colname_order = order(tissue_names, genus_names, specific_names)
    tc = tc[, colname_order]
    return(tc)
}

cleanY = function(y, mod, svs) {
    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    P = ncol(mod)
    return(y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P), ]))
}

tissue_mean = function(tc, sra, selected_tissues = NA, balance.bp = FALSE) {
    # if the data are SVA-corrected, balance.bp would not be necessary because project-specfic effects
    # were removed in SVA already.
    if (all(is.na(selected_tissues))) {
        sp_tissues = unique(sra[['tissue']])
    } else {
        sp_tissues = selected_tissues[selected_tissues %in% unique(sra[['tissue']])]
    }
    tc_ave = data.frame(matrix(rep(NA, length(sp_tissues) * nrow(tc)), nrow = nrow(tc)))
    colnames(tc_ave) = sp_tissues
    rownames(tc_ave) = rownames(tc)
    for (tissue in sp_tissues) {
        run_tissue = sra[(sra[['tissue']] == tissue),'run']
        run_tissue = run_tissue[run_tissue %in% colnames(tc)]
        if (length(run_tissue) == 1) {
            exp_tissue = tc[, run_tissue]
        } else {
            if (balance.bp) {
                is_run = (sra[['run']] %in% colnames(tc))
                is_tissue = (sra[['tissue']] == tissue)
                is_no_exclusion = (sra[['exclusion']] == "no")
                bps = unique(sra[is_run & is_tissue & is_no_exclusion, "bioproject"])
                df_tmp = data.frame(matrix(rep(NA, nrow(tc) * length(bps)), nrow = nrow(tc), ncol = length(bps)))
                colnames(df_tmp) = bps
                for (bp in bps) {
                    is_bp = (sra[['bioproject']] == bp)
                    is_tissue = (sra[['tissue']] == tissue)
                    is_no_exclusion = (sra[['exclusion']] == "no")
                    sra_ids = sra[is_bp & is_tissue & is_no_exclusion, "run"]
                    tc_bp = tc[, sra_ids]
                    if (class(tc_bp) == "numeric") {
                        df_tmp[bp] = tc_bp
                    } else {
                        df_tmp[bp] = rowMeans(tc_bp)
                    }
                }
                exp_tissue = rowMeans(df_tmp)
            } else {
                exp_tissue = rowMeans(tc[, run_tissue])
            }
        }
        tc_ave[, tissue] = exp_tissue
    }
    return(tc_ave)
}

tissue2tau = function(tc_tissue, rich.annotation = TRUE, unlog = FALSE) {
    if (rich.annotation) {
        cols = c("tau", "highest", "order")
    } else {
        cols = c("tau")
    }
    df_tau = data.frame(matrix(rep(NA, length(cols) * nrow(tc_tissue)), nrow = nrow(tc_tissue)))
    colnames(df_tau) = cols
    rownames(df_tau) = rownames(tc_tissue)
    if (unlog) {
        tc_tissue = exp(tc_tissue) - 1
        tc_tissue[tc_tissue < 0] = 0
    }
    xmax = apply(tc_tissue, 1, max)
    df_tau[,'tau'] = apply((1 - (tc_tissue/xmax))/(ncol(tc_tissue) - 1), 1, sum)
    if (rich.annotation) {
        tc_tissue[is.na(tc_tissue)] = 0
        for (i in 1:nrow(tc_tissue)) {
            is_nonzero = tc_tissue[i, ] > 0
            if (sum(is_nonzero) > 0) {
                exp_order = order(tc_tissue[i, is_nonzero], decreasing = TRUE)
                tissue_ordered = colnames(tc_tissue)[is_nonzero][exp_order]
                df_tau[i, "highest"] = tissue_ordered[1]
                df_tau[i, "order"] = paste(tissue_ordered, collapse = "|")
            }
        }
    }
    return(df_tau)
}

check_mapping_rate = function(tc, sra, mapping_rate_cutoff) {
    cat(paste0('Mapping rate cutoff: ', mapping_rate_cutoff*100, '%\n'))
    is_mapping_good = (sra[['mapping_rate']] >= mapping_rate_cutoff*100)
    is_mapping_good[is.na(is_mapping_good)] = TRUE
    if (any(!is_mapping_good)) {
        cat("Removed due to low mapping rate:\n")
        df_tmp = sra[!is_mapping_good,]
        for (i in rownames(df_tmp)) {
            sra_id = df_tmp[i,'run']
            mapping_rate = df_tmp[i,'mapping_rate']
            cat(paste0(sra_id, ': mapping rate = ', mapping_rate, '%\n'))
        }
        tc = tc[, colnames(tc) %in% sra[is_mapping_good, "run"]]
    } else {
        cat("No entry removed due to low mapping rate.\n")
    }
    sra[!is_mapping_good, "exclusion"] = "low_mapping_rate"
    return(list(tc = tc, sra = sra))
}

check_within_tissue_correlation = function(tc, sra, dist_method, min_dif, selected_tissues, one_out_per_iter = TRUE, correlation_threshold) {
    out = tc_sra_intersect(tc, sra)
    tc = out[["tc"]]
    sra2 = out[["sra"]]
    
    sra2[,'num_other_run_same_bp_tissue'] = 0
    selected_tissues = selected_tissues[selected_tissues %in% unique(sra2[['tissue']])]
    num_tissue = length(selected_tissues)
    exclude_runs = c()
    for (sra_run in colnames(tc)) {
        is_sra = (sra2[['run']] == sra_run)
        my_tissue = sra2[is_sra, "tissue"]
        my_bioproject = sra2[is_sra, "bioproject"]
        is_not_my_bp = (sra2[['bioproject']] != my_bioproject)
        is_my_tissue = (sra2[['tissue']] == my_tissue)
        run_other_bp = sra2[(is_not_my_bp | !is_my_tissue), "run"]
        run_other_bp = run_other_bp[run_other_bp %in% colnames(tc)]
        tc_other_bp = tc[, run_other_bp]
        num_other_run_same_bp_tissue = length(unique(sra2[(is_not_my_bp & is_my_tissue), "bioproject"]))
        sra2[is_sra, "num_other_run_same_bp_tissue"] = num_other_run_same_bp_tissue
        sra2_other_bp <- sra2[sra2[['run']] %in% run_other_bp,]
        
        # If one tissue is completely sourced from the same bioproject, we can't remove the whole bioproject for tc_ave_other_bp

        num_other_bp_same_tissue = sum(sra2_other_bp[['tissue']] == my_tissue, na.rm=TRUE)
        if (num_other_bp_same_tissue == 0) {
            tc_ave_other_bp = tissue_mean(tc, sra2, selected_tissues)
        } else {
            tc_ave_other_bp = tissue_mean(tc_other_bp, sra2, selected_tissues)
        }
        tc_ave = tissue_mean(tc, sra2, selected_tissues)
        coef = c()
        coef_other_bp = c()
        for (tissue in selected_tissues) {
            tmp_coef = cor(tc[, sra_run], tc_ave[, tissue], method = dist_method)

            if ((num_other_bp_same_tissue == 0) & (tissue == my_tissue)) {
              tmp_coef_other_bp = cor(tc[, sra_run], tc_ave[, tissue], method = dist_method)
            } else {
              tmp_coef_other_bp = cor(tc[, sra_run], tc_ave_other_bp[, tissue], method = dist_method)
            }
            
            if (tissue == my_tissue) {
                tmp_coef = tmp_coef - min_dif
                tmp_coef_other_bp = tmp_coef_other_bp - min_dif
            }
            coef = c(coef, tmp_coef)
            coef_other_bp = c(coef_other_bp, tmp_coef_other_bp)
        }
        names(coef) = selected_tissues
        names(coef_other_bp) = selected_tissues
        if (max(coef) != coef[my_tissue]) {
            cat('Registered as a candidate for exclusion. Better correlation to other categories:', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
        if (coef_other_bp[my_tissue] < correlation_threshold) {
            cat('Registered as a candidate for exclusion. Low within-category correlation:', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
    }
    if (length(exclude_runs)) {
      if(one_out_per_iter == TRUE){
          cat("Excluding only one outlier per bioproject or same tissue. \n")
          exclude_run_bps_and_tissue = sra2[(sra2$run %in% exclude_runs), c("bioproject", "run", "tissue")]
          first_bp_hit = exclude_run_bps_and_tissue[match(unique(exclude_run_bps_and_tissue$bioproject), exclude_run_bps_and_tissue$bioproject),]
          first_same_tissue_hit = exclude_run_bps_and_tissue[match(unique(exclude_run_bps_and_tissue$tissue), exclude_run_bps_and_tissue$tissue),]
          # if a first_same_tissue_hit is part of the same bioproject as the other removal candidates, ommit the same tissue candidates
          if(any(first_same_tissue_hit$bioproject %in% first_bp_hit$bioproject)){
              exclude_runs_tmp = c(first_bp_hit$run,first_same_tissue_hit[!first_same_tissue_hit$bioproject %in% first_bp_hit$bioproject]$run)
          }else{
              exclude_runs_tmp = c(first_bp_hit$run,first_same_tissue_hit$run)
          }
          exclude_runs = unique(exclude_runs_tmp)
          # TODO This boolean vector should be all TRUE by definition. ???: exclude_run_bps_and_tissue$run %in% exclude_runs
          exclude_bps = exclude_run_bps_and_tissue[exclude_run_bps_and_tissue$run %in% exclude_runs, "bioproject"]
      } else {
        exclude_run_bps = sra2[(sra2[['run']] %in% exclude_runs), c("bioproject", "run", "num_other_run_same_bp_tissue")]
        exclude_bp_counts = data.frame(table(exclude_run_bps[['bioproject']]))
        exclude_run_bps = merge(exclude_run_bps, exclude_bp_counts, by.x = "bioproject", by.y = "Var1")
        exclude_run_bps = exclude_run_bps[order(exclude_run_bps[['num_other_run_same_bp_tissue']], exclude_run_bps[['Freq']]),]
        rownames(exclude_run_bps) = 1:nrow(exclude_run_bps)
        min_other_run_same_bp_tissue = exclude_run_bps[1, "num_other_run_same_bp_tissue"]
        semimin_bp_count = exclude_run_bps[1, "Freq"]
        cat("minimum number of other BioProjects within tissue:", min_other_run_same_bp_tissue, "\n")
        cat("semi-minimum count of exclusion-candidate BioProjects:", semimin_bp_count, "\n")
        conditions = (exclude_run_bps[['Freq']] == semimin_bp_count)
        conditions = conditions & (exclude_run_bps[['num_other_run_same_bp_tissue']] == min_other_run_same_bp_tissue)
        exclude_bps = unique(exclude_run_bps[conditions, "bioproject"])
        exclude_runs = exclude_run_bps[(exclude_run_bps[['bioproject']] %in% exclude_bps), "run"]
      }
        
    }
    if (length(exclude_runs)) {
        cat("Partially removed BioProjects due to low within-tissue correlation:\n")
        print(exclude_bps)
        cat("Removed Runs due to low within-tissue correlation:\n")
        print(exclude_runs)
    }
    tc = tc[, !colnames(tc) %in% exclude_runs]
    sra[(sra[['run']] %in% exclude_runs), "exclusion"] = "low_within_tissue_correlation"
    return(list(tc = tc, sra = sra))
}

sva_subtraction = function(tc, sra) {
    out = tc_sra_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    out = sort_tc_and_sra(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    out = remove_nonexpressed_gene(tc)
    tc = out[["tc_ex"]]
    tc_ne = out[["tc_ne"]]
    mod = model.matrix(~tissue, data = sra)
    mod0 = model.matrix(~1, data = sra)
    set.seed(1)
    sva1 = try(sva(dat = as.matrix(tc), mod = mod, mod0 = mod0, B = 10))
    if (class(sva1) != "try-error") {
        cat("SVA correction was correctly performed.\n")
        tc = cleanY(y = tc, mod = mod, svs = sva1[['sv']])
        tc = rbind(tc, tc_ne)
    } else {
        cat("SVA correction failed. Returning the original transcriptome data.\n")
        tc = rbind(tc, tc_ne)
    }
    return(list(tc = tc, sva = sva1))
}

map_color = function(redundant_variables, c) {
    uniq_var = unique(redundant_variables)
    uniq_col = rainbow_hcl(length(uniq_var), c = c)
    df_unique = data.frame(var = uniq_var, col = uniq_col, stringsAsFactors = FALSE)
    df_redundant = data.frame(var = redundant_variables, order = seq(1, length(redundant_variables)),
                              stringsAsFactors = FALSE)
    df_redundant = merge(df_redundant, df_unique, by = "var", all.x = TRUE, stringsAsFactors = FALSE)
    df_redundant = df_redundant[order(df_redundant[['order']]), ]
    return(df_redundant[['col']])
}

draw_heatmap = function(sra, tc_dist_matrix, legend = TRUE, fontsize = 7) {
    bp_fac = factor(sub(";.*", "", sra[, c("bioproject")]))
    tissue_fac = factor(sra[, c("tissue")])
    ann_label = data.frame(bioproject = bp_fac, tissue = tissue_fac)
    bp_col_uniq = unique(sra[order(sra[['bioproject']]), 'bp_color'])
    tissue_col_uniq = unique(sra[order(sra[['tissue']]), 'tissue_color'])
    ann_color = list(bioproject = bp_col_uniq, tissue = tissue_col_uniq)
    breaks = c(0, seq(0.3, 1, 0.01))
    aheatmap(tc_dist_matrix, color = "-RdYlBu2:71", Rowv = NA, Colv = NA, revC = TRUE, legend = TRUE,
             breaks = breaks, annCol = ann_label, annRow = ann_label, annColors = ann_color, annLegend = legend,
             fontsize = fontsize)
}

color_children2parent = function(node) {
    if (length(node) == 2) {
        child1_color = attributes(node[[1]])[['edgePar']][["col"]]
        child2_color = attributes(node[[2]])[['edgePar']][["col"]]
        if ((!is.null(child1_color)) & (!is.null(child2_color))) {
            if (child1_color == child2_color) {
                attributes(node)[['edgePar']][["col"]] = child1_color
            }
        }
    }
    return(node)
}

draw_dendrogram = function(sra, tc_dist_dist, fontsize = 7) {
    dend <- as.dendrogram(hclust(tc_dist_dist))
    dend_colors = sra[order.dendrogram(dend),'tissue_color']
    labels_colors(dend) <- dend_colors
    dend_labels <- sra[order.dendrogram(dend), 'run']
    dend <- color_branches(dend, labels = dend_labels, col = dend_colors)
    dend <- set(dend, "branches_lwd", 1)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    cex.xlab = min(fontsize, max(0.2, 0.5/log10(nrow(sra))))
    par(cex = cex.xlab)
    plot(dend, las = 1, axes = FALSE)
    par(cex = 1)
    axis(side = 2, line = 0, las = 1)
    mtext("Distance", side = 2, line = 8.5, outer = FALSE)
    n = nrow(sra)
    symbols(
      1:n,
      rep(0, n),
      circles = rep(1, n),
      add = TRUE,
      inches = 0.02,
      xpd = TRUE,
      lwd = 1,
      bg = sra[order.dendrogram(dend), 'tissue_color'],
      fg = sra[order.dendrogram(dend), 'bp_color']
    )
}

draw_dendrogram_pvclust = function(sra, tc, nboot, pvclust_file, fontsize = 7) {
    dist_fun = function(x) {
        Dist(t(x), method = "pearson")
    }
    sp = sub(" ", "_", sra[['scientific_name']][1])
    if (file.exists(pvclust_file)) {
        if (file.info(pvclust_file)[['size']]) {
            print("pvclust file found.")
            load(pvclust_file)
        }
    } else {
        print("no pvclust file found. Start bootstrapping.")
        result = pvclust(tc, method.dist = dist_fun, method.hclust = "average", nboot = nboot, parallel = FALSE)  # UPGMA
        save(result, file = pvclust_file)
    }
    dend = as.dendrogram(result)
    dend_colors = sra[order.dendrogram(dend), 'tissue_color']
    labels_colors(dend) = dend_colors
    dend_labels = sra[order.dendrogram(dend), 'run']
    dend = color_branches(dend, labels = dend_labels, col = dend_colors)
    dend = set(dend, "branches_lwd", 2)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    cex.xlab = min(0.2 + 1/log10(fontsize))
    par(cex = cex.xlab)
    plot(dend, las = 1, ylab = "Distance", cex.axis = 1/cex.xlab, cex.lab = 1/cex.xlab)
    par(cex = 1)
    n = nrow(sra)
    symbols(
      1:n,
      rep(0, n),
      circles = rep(1, n),
      add = TRUE,
      inches = 0.04,
      xpd = TRUE,
      lwd = 2,
      bg = sra[order.dendrogram(dend), 'tissue_color'],
      fg = sra[order.dendrogram(dend), 'bp_color']
    )
    text(result, print.num = FALSE, cex = 1, col.pv = "black")
}

draw_pca = function(sra, tc_dist_matrix, fontsize = 7) {
    set.seed(1)
    pca = prcomp(tc_dist_matrix)
    xlabel = paste0("PC 1 (", round(summary(pca)[['importance']][2, 1] * 100, digits = 1), "%)")
    ylabel = paste0("PC 2 (", round(summary(pca)[['importance']][2, 2] * 100, digits = 1), "%)")
    plot(
      pca[['rotation']][, 1],
      pca[['rotation']][, 2],
      pch = 21,
      cex = 2,
      lwd = 1,
      bg = sra[['tissue_color']],
      col = sra[['bp_color']],
      xlab = xlabel,
      ylab = ylabel,
      las = 1
    )
    # plot(pca$x[,1], pca$x[,2], pch=21, cex=2, lwd=2, bg=sra$tissue_color, col=sra$bp_color, main=title,
    # xlab=xlabel, ylab=ylabel, las=1)
}

draw_mds = function(sra, tc_dist_dist, fontsize = 7) {
    set.seed(1)
    try_out = tryCatch({
        isoMDS(tc_dist_dist, k = 2, maxit = 100)
    }, error = function(a) {
        return("MDS failed.")
    })
    if (mode(try_out) == "character") {
        cat("MDS failed.\n")
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    } else {
        mds <- try_out
        plot(
          mds[['points']][, 1],
          mds[['points']][, 2],
          pch = 21,
          cex = 2,
          lwd = 1,
          bg = sra[['tissue_color']],
          col = sra[['bp_color']],
          xlab = "MDS dimension 1",
          ylab = "MDS dimension 2",
          las = 1
        )
    }
}

draw_tsne = function(sra, tc, fontsize = 7) {
    perplexity = min(30, floor(nrow(sra)/4))
    set.seed(1)
    out_tsne = Rtsne(as.matrix(t(tc)), theta = 0, check_duplicates = FALSE, verbose = FALSE, perplexity = perplexity,
                     dims = 2)
    try_out = tryCatch({
        plot(
          out_tsne[['Y']][, 1],
          out_tsne[['Y']][, 2],
          pch = 21,
          cex = 2,
          lwd = 1,
          bg = sra[['tissue_color']],
          col = sra[['bp_color']],
          xlab = "t-SNE dimension 1",
          ylab = "t-SNE dimension 2",
          las = 1
        )
    }, error = function(a) {
        return("t-SNE plot failed.")
    })
    if (mode(try_out) == "character") {
        cat("t-SNE failed.\n")
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    }
}

draw_sva_summary = function(sva_out, tc, sra, fontsize) {
    if ((is.null(sva_out)) | (class(sva_out) == "try-error")) {
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        df = NA
    } else {
        out = tc_sra_intersect(tc, sra)
        tc = out[["tc"]]
        sra = out[["sra"]]
        out = sort_tc_and_sra(tc, sra)
        tc = out[["tc"]]
        sra = out[["sra"]]
        cols = c("tissue","bioproject","lib_selection","instrument","mapping_rate")
        label_cols = c("tissue","BioProject","library selection","instrument","mapping rate")

        num_sv = sva_out[['n.sv']]
        df = data.frame(matrix(NA, num_sv, length(cols)))
        colnames(df) = cols
        rownames(df) = paste0("SV", 1:nrow(df))
        for (i in 1:length(cols)) {
            for (j in 1:num_sv) {
                if (length(unique(sra[, cols[i]])) == 1) {
                    df[j, i] = NA
                } else {
                    lm_summary = summary(lm(sva_out[['sv']][, j] ~ sra[, cols[i]]))
                    df[j, i] = lm_summary[['adj.r.squared']]
                }
            }
        }
        colnames(df) = label_cols
        breaks = seq(0, 1, 0.02)
        colors = colorRampPalette(c("blue", "yellow", "red"))(length(breaks))
        df2 = t(df)
        df2[df2 < 0] = 0
        aheatmap(df2, color = colors, Rowv = NA, Colv = NA, revC = TRUE, breaks = breaks, fontsize = fontsize)
    }
    return(df)
}

draw_boxplot = function(sra, tc_dist_matrix, fontsize = 7) {
    is_same_bp = outer(sra[['bioproject']], sra[['bioproject']], function(x, y) {
        x == y
    })
    is_same_tissue = outer(sra[['tissue']], sra[['tissue']], function(x, y) {
        x == y
    })
    plot(c(0.5, 4.5), c(0, 1), type = "n", xlab = "", ylab = "Pearson's correlation\ncoefficient", las = 1,
         xaxt = "n")
    boxplot(tc_dist_matrix[(!is_same_bp) & (!is_same_tissue)], at = 1, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(is_same_bp) & (!is_same_tissue)], at = 2, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(!is_same_bp) & (is_same_tissue)], at = 3, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(is_same_bp) & (is_same_tissue)], at = 4, add = TRUE, col = "gray", yaxt = "n")
    labels = c("bw\nbw", "bw\nwi", "wi\nbw", "wi\nwi")
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = "Tissue\nBioProject", padj = 0.5, hadj = 1, tick = FALSE)

}

draw_tau_histogram = function(tc, sra, selected_tissues, fontsize = 7) {
    df_tau = tissue2tau(tissue_mean(tc, sra, selected_tissues), rich.annotation = FALSE, unlog = TRUE)
    hist_out = hist(df_tau[['tau']], breaks = seq(0, 1, 0.05), las = 1, xlab = "Tau (expression specificity)",
                    ylab = "Gene count", main = "", col = "gray")
    num_noexp = sum(is.na(df_tau[['tau']]))
    num_all = nrow(df_tau)
    # num_exp = nrow(df_tau) - num_noexp text_noexp = paste('Expressed genes:', num_exp,
    # '\nNon-expressed genes:', num_noexp)
    text_noexp = paste0("Excluded due to\nno expression:\n", num_noexp, "/", num_all, " genes")
    text(0, max(hist_out[['counts']]) * 0.85, text_noexp, pos = 4)
}

draw_exp_level_histogram = function(tc, sra, selected_tissues, fontsize = 7, transform_method) {
    tc_tissue = tissue_mean(tc, sra, selected_tissues)
    xmax = apply(tc_tissue, 1, max)
    xmax[xmax < 0] = 0
    xmax[xmax > 15] = 15
    breaks = seq(0, 15, 1)
    hist_out = hist(xmax, breaks = breaks, las = 1, xlab = paste0("Max expression (log ",transform_method,"+1)"), ylab = "Gene count",
                    main = "", col = "gray")
}

draw_legend = function(sra, new = TRUE, pos = "center", fontsize = 7, nlabel.in.col) {
    if (new) {
        plot.new()
    }
    tissue_unique = unique(sra[['tissue']])
    bp_unique = unique(sub(";.*", "", sra[['bioproject']]))
    tissue_color_unique = unique(sra[['tissue_color']])
    bp_color_unique = unique(sra[['bp_color']])
    ncol = ceiling((length(tissue_unique) + length(bp_unique) + 2)/nlabel.in.col)
    legend_text = c("Tissue", as.character(tissue_unique), "", "BioProject", as.character(bp_unique))
    legend_color = c(rgb(1, 1, 1, 0), rep(rgb(1, 1, 1, 0), length(tissue_color_unique)), rgb(1, 1, 1,
                                                                                             0), rgb(1, 1, 1, 0), bp_color_unique)
    legend_bg = c(rgb(1, 1, 1, 0), tissue_color_unique, rgb(1, 1, 1, 0), rgb(1, 1, 1, 0), rep(rgb(1,
                                                                                                  1, 1, 0), length(bp_color_unique)))
    legend_font = c(2, rep(1, length(tissue_color_unique)), 1, 2, rep(1, length(bp_color_unique)))
    legend(pos, legend = legend_text, pch = 21, lwd = 1, lty = 0, col = legend_color, pt.bg = legend_bg,
           text.font = legend_font, ncol = ncol, bty = "n")
}

save_plot = function(tc, sra, sva_out, dist_method, file, selected_tissues, fontsize = 7, transform_method) {
    out = tc_sra_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    sra = add_color_to_sra(sra, selected_tissues)
    out = sort_tc_and_sra(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    pdf(paste0(dir_pdf,'/',file, ".pdf"), height = 8, width = 7.2, fonts = "Helvetica", pointsize = fontsize)
    layout_matrix = matrix(c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1,
                             2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1,
                             1, 1, 1, 1, 1, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3,
                             3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9,
                             9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10,
                             10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), 14, 12,
                           byrow = TRUE)
    layout(layout_matrix)
    tc_dist_matrix = cor(tc, method = dist_method)
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    tc_dist_dist = Dist(t(tc), method = dist_method) + 1e-09
    tc_dist_dist[is.na(tc_dist_dist)] = 1
    par(mar = c(6, 6, 1, 0))
    draw_dendrogram(sra, tc_dist_dist, fontsize)
    par(mar = c(0, 0, 0, 0))
    draw_heatmap(sra, tc_dist_matrix, legend = FALSE)
    # draw_dendrogram(sra, tc, nboot=1000, cex.xlab=0.6, pvclust_file=paste0(file, '.pvclust.RData'))
    par(mar = c(4, 4, 0.1, 1))
    draw_pca(sra, tc_dist_matrix, fontsize)
    par(mar = c(4, 4, 0.1, 1))
    draw_mds(sra, tc_dist_dist, fontsize)
    par(mar = c(4, 4, 0.1, 1))
    draw_tsne(sra, tc, fontsize)
    par(mar = c(4, 5, 0.1, 1))
    draw_boxplot(sra, tc_dist_matrix, fontsize)
    par(mar = c(4, 4, 1, 1))
    draw_exp_level_histogram(tc, sra, selected_tissues, fontsize, transform_method)
    par(mar = c(4, 4, 1, 1))
    draw_tau_histogram(tc, sra, selected_tissues, fontsize)
    par(mar = rep(0.1, 4))
    df_r2 = draw_sva_summary(sva_out, tc, sra, fontsize)
    if (!all(is.na(df_r2))) {
        write.table(df_r2, paste0(dir_tsv,'/',file, ".r2.tsv"), sep = "\t", row.names = FALSE)
    }
    par(mar = rep(0.1, 4))
    draw_legend(sra, new = TRUE, pos = "center", fontsize = fontsize, nlabel.in.col = 8)
    graphics.off()
}

transform_raw_to_fpkm = function(counts, effective_lengths) {
    res = counts / effective_lengths / sum(counts) * 1e+09
    return(as.data.frame(res))
}

transform_raw_to_tpm = function(counts, len) {
    x <- counts/len
    return(t(t(x) * 1e+06/colSums(x)))
}

get_tmm_scaled_fpkm = function(dat, df_nf, efflen) {
    dat_fpkm = dat
    dat_fpkm[, ] = NA
    for (sample in colnames(dat)) {
        if (sample %in% rownames(df_nf)) {
            dat_sample = data.frame(dat[, sample])
            colnames(dat_sample) = sample
            nf_sample = df_nf[sample, ]
            dl_sample = list(counts = dat_sample, samples = nf_sample)
            class(dl_sample) = "DGEList"
            col_fpkm = edgeR::rpkm(dl_sample, gene.length = efflen[, sample], normalized.lib.sizes = TRUE)
            dat_fpkm[, sample] = col_fpkm
        } else {
            dat_fpkm[, sample] = dat[, sample]
        }
    }
    return(dat_fpkm)
}


####################END OF FUNCTION DECLARATION####################################################

################################################# START OF SVA CORRECTION ########################################################


fontsize = 7

# read SRA table
sra_all = read.table(srafile, sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "",
                     stringsAsFactors = FALSE, check.names=FALSE)
for (col in c('instrument','bioproject')) {
    is_missing = (sra_all[,col] == "")|(is.na(sra_all[,col]))
    sra_all[is_missing, col] = "not_provided"
}

tc <- read.table(infile, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, row.names = 1, check.names=FALSE)

tc_eff_length <- read.table(eff_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, check.names=FALSE)
# row.names(tc)<-tc[,1] tc <- tc[,-1]

# read transcriptome
scientific_name = unique(sra_all[sra_all[['run']] %in% colnames(tc), "scientific_name"])
is_sp = (sra_all[,'scientific_name'] == scientific_name)
is_tissue = (sra_all[,'tissue'] %in% selected_tissues)
cat('Number of SRA runs for this species:', sum(is_sp), '\n')
cat('Number of SRA runs for selected tisssues:', sum(is_tissue), '\n')
sra = sra_all[(is_sp & is_tissue),]
conditions = (sra[['exclusion']] == "no") & (!sra[['run']] %in% colnames(tc))
if (any(conditions)) {
    cat("Failed quantification:", sra[conditions, "run"], "\n")
    sra[conditions, "exclusion"] = "failed_quantification"
}

is_not_excluded = (sra[['exclusion']]=='no')
cat('Number of non-excluded SRA runs (exclusion=="no"):', sum(is_not_excluded), '\n')
tc = tc[,sra[is_not_excluded,'run']]
out = sort_tc_and_sra(tc, sra) ; tc = out[["tc"]] ; sra = out[["sra"]]

# log transform AFTER mappingrate
row.names(tc_eff_length) <- tc_eff_length[, 1]
tc_eff_length <- tc_eff_length[, colnames(tc)]
if (transform_method == "fpkm") {
    tc <- transform_raw_to_fpkm(tc, tc_eff_length)
    tc <- log(tc + 1)
}
if (transform_method == "tpm") {
    tc <- transform_raw_to_tpm(tc, tc_eff_length)
    tc <- log(tc + 1)
}


file_name = paste0(dir_tsv,'/',sub(" ", "_", scientific_name), ".uncorrected.tc.tsv")
write.table(tc, file = file_name,
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
tc_uncorrected = tc
tc_tissue_uncorrected = tissue_mean(tc, sra, selected_tissues)
file_name = paste0(dir_tsv,'/',sub(" ", "_", scientific_name), ".uncorrected.tissue.mean.tsv")
write.table(tc_tissue_uncorrected, file = file_name,
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


round = 0
sva_out = NULL
tc_sva = NULL
save_plot(tc, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".original"),
          selected_tissues, fontsize, transform_method)
out = sva_subtraction(tc, sra)
tc_sva = out[["tc"]]
sva_out = out[["sva"]]
save(sva_out, file = paste0(dir_rdata,'/', sub(" ", "_", scientific_name), ".sva.", round, ".RData"))
save_plot(tc_sva, sra, sva_out, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".original.sva"),
          selected_tissues, fontsize, transform_method)

round = 1
sva_out = NULL
tc_sva = NULL
out = check_mapping_rate(tc, sra, mapping_rate_cutoff)
tc = out[["tc"]]
sra = out[["sra"]]
tc = tc[, sra[sra[['exclusion']] == "no", "run"]]
save_plot(tc, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".mapping_cutoff"),
          selected_tissues, fontsize, transform_method)
out = sva_subtraction(tc, sra)
tc_sva = out[["tc"]]
sva_out = out[["sva"]]
save(sva_out, file = paste0(dir_rdata,'/',sub(" ", "_", scientific_name), ".sva.", round, ".RData"))
save_plot(tc_sva, sra, sva_out, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".mapping_cutoff.sva"),
          selected_tissues, fontsize, transform_method)

round = 2
end_flag = 0
while (end_flag == 0) {
    cat("iteratively checking within-tissue correlation, round:", round, "\n")
    tc_cwtc = NULL
    num_run_before = sum(sra[['exclusion']] == "no")
    out = check_within_tissue_correlation(tc, sra, dist_method, min_dif, selected_tissues, one_outlier_per_iteration, correlation_threshold)
    tc_cwtc = out[["tc"]]
    sra = out[["sra"]]
    num_run_after = sum(sra[['exclusion']] == "no")
    if ((num_run_before == num_run_after) | (plot_intermediate)) {
        sva_out = NULL
        tc_sva = NULL
        out = sva_subtraction(tc_cwtc, sra)
        tc_sva = out[["tc"]]
        sva_out = out[["sva"]]
        save(sva_out, file = paste0(dir_rdata,'/',sub(" ", "_", scientific_name), ".sva.", round, ".RData"))
        save_plot(tc_cwtc, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round,
                                                          ".correlation_cutoff"), selected_tissues, fontsize, transform_method)
        save_plot(tc_sva, sra, sva_out, dist_method, paste0(sub(" ", "_", scientific_name), ".", round,
                                                            ".correlation_cutoff.sva"), selected_tissues, fontsize, transform_method)
    }
    # out = check_within_tissue_correlation(tc_sva, sra, dist_method, min_dif, selected_tissues) ; tc_sva= out[['tc']] ; sra = out[['sra']]
    
    cat("round:", round, ": # before =", num_run_before, ": # after =", num_run_after, "\n\n")
    # tc = tc[,colnames(tc_sva)]
    if (num_run_before == num_run_after) {
        end_flag = 1
    }
    tc = tc_cwtc
    round = round + 1
}
cat("finished checking within-tissue correlation.\n")
file = paste0(dir_tsv,'/',sub(" ", "_", scientific_name), ".sra.tsv")
write.table(sra, file = file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
file = paste0(dir_tsv,'/',sub(" ", "_", scientific_name), ".tc.tsv")
write.table(tc_sva, file = file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
tc_tissue = tissue_mean(tc_sva, sra, selected_tissues)
file = paste0(dir_tsv,'/',sub(" ", "_", scientific_name), ".tissue.mean.tsv")
write.table(tc_tissue, file = file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
tc_tau = tissue2tau(tc_tissue, rich.annotation = TRUE, unlog = TRUE)
file = paste0(dir_tsv,'/',sub(" ", "_", scientific_name), ".tau.tsv")
write.table(tc_tau, file = file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
cat("Done!\n")


#
#
# par_sva = function(i, species_names, dir_count, dir_eff_length) {
#   species = species_names[i]
#   cat("Files found for", species, ": ")
#   species = sub(' ', '_', species)
#   infile = list.files(dir_count, pattern = species)[1]
#   infile = paste0(dir_count, infile)
#   if(exists("dir_eff_length")){
#     eff_file = list.files(dir_eff_length, pattern = species)[1]
#     eff_file = paste0(dir_eff_length, eff_file)
#     cat(infile, eff_file, "\n")
#   }
#   else{
#     cat(infile, "\n")
#     eff_file = NA
#   }
#
#   tc <- read.table(infile, sep="\t", stringsAsFactors=FALSE, header = T, quote="", fill =F)
#   tc_eff_length <- read.table(eff_file, sep="\t", stringsAsFactors=FALSE, header = T, quote="", fill =F)
#   row.names(tc)<-tc[,1]
#   tc<-tc[,-1]
#   row.names(tc_eff_length)<-tc_eff_length[,1]
#   tc_eff_length<-tc_eff_length[,-1]
#   tc_eff_length = tc_eff_length[,colnames(tc)]
#   # calculate tpm and log transform transcriptome (adds pseudocount +1), important to do AFTER get_mapping_rate()
#   tc = transform_raw_to_fpkm(tc, tc_eff_length)
#   main_sva()
# }
#


#}
#
# if(mode == 'msm'){
#   #########################
#   # read SRA table
#   sra_all = read.table(srafile, sep="\t", header=TRUE, quote="", fill=TRUE, comment.char="", stringsAsFactors=FALSE)
#   sra_all$instrument[sra_all$instrument==""] = "not_provided"
#   sra_all$bioproject[sra_all$bioproject==""] = "not_provided"
#
#   #setup parallel backend to use many processors
#   cores=detectCores()
#   #cl <- makeCluster(cores[1]-1) #not to overload your computer
#   cl <- parallel::makeCluster(cores[1]-1, setup_strategy = "sequential")
#
#   registerDoParallel(cl)
#   req_packages = c("Biobase",
#                    "pcaMethods",
#                    "colorspace",
#                    "RColorBrewer",
#                    "sva",
#                    "MASS",
#                    "NMF",
#                    "dendextend",
#                    "amap",
#                    "pvclust",
#                    "Rtsne")
#
#   species_names = unique(sra_all[,'scientific_name'])
#   cat("Species in the provided metadata file: ", species_names, "\n")
#
#
#
#
#
#   out = list()
#   foreach(i=1:length(species_names), .packages = req_packages) %dopar%{
#     out[[species_names[i]]] = par_sva(i, species_names, dir_count, dir_eff_length)
#   }
#
#
#
#  stopCluster()
#   #
# foreach(i=1:length(species_names), .packages = req_packages) %dopar%{
#   species = species_names[i]
#   cat("Files found for", species, ": ")
#   species = sub(' ', '_', species)
#   infile = list.files(dir_count, pattern = species)[1]
#   infile = paste0(dir_count, infile)
#   if(exists("dir_eff_length")){
#   eff_file = list.files(dir_eff_length, pattern = species)[1]
#   eff_file = paste0(dir_eff_length, eff_file)
#   cat(infile, eff_file, "\n")
#   }
#   else{
#   cat(infile, "\n")
#   eff_file = NA
#   }
#
#   tc <- read.table(infile, sep="\t", stringsAsFactors=FALSE, header = T, quote="", fill =F)
#   tc_eff_length <- read.table(eff_file, sep="\t", stringsAsFactors=FALSE, header = T, quote="", fill =F)
#   row.names(tc)<-tc[,1]
#   tc<-tc[,-1]
#   row.names(tc_eff_length)<-tc_eff_length[,1]
#   tc_eff_length<-tc_eff_length[,-1]
#   tc_eff_length = tc_eff_length[,colnames(tc)]
#   # calculate tpm and log transform transcriptome (adds pseudocount +1), important to do AFTER get_mapping_rate()
#   tc = transform_raw_to_fpkm(tc, tc_eff_length)
#   main_sva()
# }
#
# stopCluster()
#

#}
