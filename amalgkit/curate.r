#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(pcaMethods, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(colorspace, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(MASS, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(NMF, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(dendextend, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(amap, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(pvclust, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(Rtsne, quietly = TRUE)))

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
log_prefix = "transcriptome_curation.r:"
cat(log_prefix, "mode =", debug_mode, "\n")
if (debug_mode == "debug") {
    out_dir = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
    metadata_path = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out/cstmm/metadata.tsv'
    est_counts_path = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out/cstmm/Actinidia_chinensis/Actinidia_chinensis_cstmm_counts.tsv'
    eff_length_path = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out/cstmm/Actinidia_chinensis/Actinidia_chinensis_eff_length.tsv'
    dist_method = "pearson"
    mapping_rate_cutoff = .20
    min_dif = 0
    plot_intermediate = as.logical(0)
    selected_curate_groups = c("root", "flower", "leaf")
    transform_method = "log2p1-fpkm"
    one_outlier_per_iteration = as.logical(0)
    correlation_threshold = 0.3
    batch_effect_alg = 'sva'
    dist_method = "pearson"
    clip_negative = as.logical(1)
    maintain_zero = as.logical(1)
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
    setwd(file.path(out_dir, 'curate'))
} else if (debug_mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    est_counts_path = args[1]
    metadata_path = args[2]
    out_dir = args[3]
    eff_length_path = args[4]
    dist_method = args[5]
    mapping_rate_cutoff = as.numeric(args[6])
    min_dif = as.numeric(args[7])
    plot_intermediate = as.logical(as.integer(args[8]))
    selected_curate_groups = strsplit(args[9], "\\|")[[1]]
    transform_method = args[10]
    one_outlier_per_iteration = as.integer(args[11])
    correlation_threshold = as.numeric(args[12])
    batch_effect_alg = args[13]
    clip_negative = as.logical(as.integer(args[14]))
    maintain_zero = as.logical(as.integer(args[15]))
    r_util_path = file.path(args[16])
    skip_curation_flag = as.logical(as.integer(args[17]))
}
cat('est_counts_path:', est_counts_path, "\n")
cat('metadata_path:', metadata_path, "\n")
cat('out_dir:', out_dir, "\n")
cat('eff_length_path:', eff_length_path, "\n")
cat('dist_method:', dist_method, "\n")
cat('mapping_rate_cutoff:', mapping_rate_cutoff, "\n")
cat('min_dif:', min_dif, "\n")
cat('plot_intermediate:', plot_intermediate, "\n")
cat('selected_curate_groups:', selected_curate_groups, "\n")
cat('transform_method:', transform_method, "\n")
cat('one_outlier_per_iteration:', one_outlier_per_iteration, "\n")
cat('correlation_threshold:', correlation_threshold, "\n")
cat('batch_effect_alg:', batch_effect_alg, "\n")
cat('clip_negative:', clip_negative, "\n")
cat('maintain_zero:', maintain_zero, "\n")
cat('r_util_path:', r_util_path, "\n")
cat('skip_curation_flag:', skip_curation_flag, "\n")

source(r_util_path)

if (batch_effect_alg == "ruvseq") {
    suppressWarnings(suppressPackageStartupMessages(library(RUVSeq, quietly = TRUE)))
} else {
    suppressWarnings(suppressPackageStartupMessages(library(sva, quietly = TRUE)))
}

tc_metadata_intersect = function(tc, sra) {
    sra_run = sra[['run']]
    tc = tc[, colnames(tc) %in% sra_run, drop=FALSE]
    sra = sra[sra[['run']] %in% colnames(tc), ]
    return(list(tc = tc, sra = sra))
}

remove_nonexpressed_gene = function(tc) {
    gene_sum = apply(tc, 1, sum)
    tc_ex = tc[gene_sum > 0, ]
    tc_ne = tc[gene_sum == 0, ]
    return(list(tc_ex = tc_ex, tc_ne = tc_ne))
}

add_color_to_metadata = function(sra, selected_curate_groups) {
    sra = sra[, (!colnames(sra) %in% c("bp_color", "sp_color", "curate_group_color"))]
    bioproject = as.character(sra[['bioproject']])
    bioproject_u = sort(unique(bioproject))
    scientific_name = as.character(sra[['scientific_name']])
    scientific_name_u = sort(unique(scientific_name))
    curate_group = as.character(sra[['curate_group']])
    curate_group_u = sort(unique(curate_group))
    if (length(selected_curate_groups) <= 8) {
        curate_group_color = brewer.pal(8, "Dark2")
        curate_group_color = curate_group_color[1:length(selected_curate_groups)] # To avoid the warning "minimal value for n is 3, returning requested palette with 3 different levels"
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
    } else if (length(selected_curate_groups) <= 12) {
        curate_group_color = brewer.pal(length(unique(curate_group)), "Paired")
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
    } else {
        curate_group_color = rainbow_hcl(length(selected_curate_groups), c = 100)
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 150)
    }
    df_curate_group = data.frame(curate_group=curate_group_u, curate_group_color=curate_group_color[1:length(curate_group_u)], stringsAsFactors = FALSE)
    df_bp = data.frame(bioproject=bioproject_u, bp_color=bp_color[1:length(bioproject_u)], stringsAsFactors = FALSE)
    df_sp = data.frame(scientific_name=scientific_name_u, sp_color=sp_color[1:length(scientific_name_u)], stringsAsFactors = FALSE)
    sra = merge(sra, df_bp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_sp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_curate_group, sort = FALSE, all.y = FALSE)
    return(sra)
}

sort_tc_and_metadata = function(tc, sra, sort_columns = c("curate_group", "scientific_name", "bioproject")) {
    for (column in rev(sort_columns)) {
        sra = sra[order(sra[[column]]), ]
    }
    sra_intersection = sra[(sra[['run']] %in% colnames(tc)),'run']
    tc = tc[, sra_intersection, drop=FALSE]
    return(list(tc = tc, sra = sra))
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
    tc = tc[, colname_order, drop=FALSE]
    return(tc)
}

cleanY = function(y, mod, svs) {
    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    P = ncol(mod)
    return(y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P), ]))
}

curate_group_mean = function(tc, sra, selected_curate_groups = NA, balance_bp = FALSE) {
    # if the data are SVA-corrected, balance_bp would not be necessary because project-specfic effects
    # were removed in SVA already.
    if (all(is.na(selected_curate_groups))) {
        sp_curate_groups = unique(sra[['curate_group']])
    } else {
        sp_curate_groups = selected_curate_groups[selected_curate_groups %in% unique(sra[['curate_group']])]
    }
    tc_ave = data.frame(matrix(rep(NA, length(sp_curate_groups) * nrow(tc)), nrow = nrow(tc)))
    colnames(tc_ave) = sp_curate_groups
    rownames(tc_ave) = rownames(tc)
    for (curate_group in sp_curate_groups) {
        exclusion_curate_group = sra[(sra[['curate_group']] == curate_group),'exclusion']
        run_curate_group = sra[(sra[['curate_group']] == curate_group),'run']
        run_curate_group = run_curate_group[run_curate_group %in% colnames(tc)]
        if (all(exclusion_curate_group!= "no")) {
          warning_message = paste0('All samples of curate_group ', curate_group, ' are marked for exclusion. This curate_group will be omitted from further analysis.')
          selected_curate_groups = selected_curate_groups[!selected_curate_groups==curate_group]
          tc_ave = tc_ave[, !names(tc_ave) %in% c(curate_group)]
          warning(warning_message)
          next
        }
        if (length(run_curate_group) == 1) {
            exp_curate_group = tc[,run_curate_group]
        } else {
            if (balance_bp) {
                is_run = (sra[['run']] %in% colnames(tc))
                is_curate_group = (sra[['curate_group']] == curate_group)
                is_no_exclusion = (sra[['exclusion']] == "no")
                bps = unique(sra[is_run & is_curate_group & is_no_exclusion, "bioproject"])
                df_tmp = data.frame(matrix(rep(NA, nrow(tc) * length(bps)), nrow = nrow(tc), ncol = length(bps)))
                colnames(df_tmp) = bps
                for (bp in bps) {
                    is_bp = (sra[['bioproject']] == bp)
                    is_curate_group = (sra[['curate_group']] == curate_group)
                    is_no_exclusion = (sra[['exclusion']] == "no")
                    sra_ids = sra[is_bp & is_curate_group & is_no_exclusion, "run"]
                    tc_bp = tc[, sra_ids]
                    if (class(tc_bp) == "numeric") {
                        df_tmp[bp] = tc_bp
                    } else {
                        df_tmp[bp] = rowMeans(tc_bp)
                    }
                }
                exp_curate_group = rowMeans(df_tmp)
            } else {
                exp_curate_group = rowMeans(tc[, run_curate_group])
            }
        }
        tc_ave[, curate_group] = exp_curate_group
    }
    return(list(tc_ave=tc_ave, selected_curate_groups=selected_curate_groups))
}

row_tau = function(row) {
  is_nonzero = row > 0
  if (sum(is_nonzero) > 0) {
    exp_order = order(row[is_nonzero], decreasing = TRUE)
    curate_group_ordered = colnames(tc_curate_group)[is_nonzero][exp_order]
    highest = curate_group_ordered[1]
    order = paste(curate_group_ordered, collapse = "|")
    return(c(highest, order))
  } else {
    return(c(NA, NA))  # Return NA values if no nonzero elements
  }
}

curate_group2tau = function(tc_curate_group, rich.annotation = TRUE, transform_method) {
    if (rich.annotation) {
        cols = c("tau", "highest", "order")
    } else {
        cols = c("tau")
    }
    df_tau = data.frame(matrix(rep(NA, length(cols) * nrow(tc_curate_group)), nrow = nrow(tc_curate_group)))
    colnames(df_tau) = cols
    rownames(df_tau) = rownames(tc_curate_group)
    if (grepl('logn-', transform_method)) {
        tc_curate_group = exp(tc_curate_group)
    } else if (grepl('log2-', transform_method)) {
        tc_curate_group = 2**tc_curate_group
    } else if (grepl('lognp1-', transform_method)) {
        tc_curate_group = exp(tc_curate_group) - 1
    } else if (grepl('log2p1-', transform_method)) {
        tc_curate_group = 2**tc_curate_group - 1
    }
    tc_curate_group[tc_curate_group < 0] = 0
    xmax = apply(tc_curate_group, 1, function(x){max(x, na.rm=TRUE)})
    df_tau[,'tau'] = apply((1 - (tc_curate_group/xmax))/(ncol(tc_curate_group) - 1), 1, sum)
    if (rich.annotation) {
        tc_curate_group[is.na(tc_curate_group)] = 0
        results = t(apply(tc_curate_group, 1, row_tau))
        df_tau[ , c("highest", "order")] = results
    }
    return(df_tau)
}

check_mapping_rate = function(tc, sra, mapping_rate_cutoff) {
    if ('mapping_rate' %in% colnames(sra)) {
        cat(paste0('Mapping rate cutoff: ', mapping_rate_cutoff*100, '%\n'))
        is_mapping_good = (sra[['mapping_rate']] > mapping_rate_cutoff*100)
        is_mapping_good[is.na(is_mapping_good)] = TRUE
        if (any(!is_mapping_good)) {
            cat("Removed due to low mapping rate:\n")
            df_tmp = sra[!is_mapping_good,]
            for (i in rownames(df_tmp)) {
                sra_id = df_tmp[i,'run']
                mapping_rate = df_tmp[i,'mapping_rate']
                cat(paste0(sra_id, ': mapping rate = ', mapping_rate, '%\n'))
            }
            tc = tc[, colnames(tc) %in% sra[is_mapping_good, "run"], drop=FALSE]
        } else {
            cat("No entry removed due to low mapping rate.\n")
        }
        sra[!is_mapping_good, "exclusion"] = "low_mapping_rate"
    } else {
        cat('Mapping rate cutoff will not be applied.\n')
    }
    return(list(tc = tc, sra = sra))
}

check_within_curate_group_correlation = function(tc, sra, dist_method, min_dif, selected_curate_groups, one_out_per_iter = TRUE, correlation_threshold) {
    if (length(selected_curate_groups)==1) {
        cat('Only one curate_group category is available. Outlier removal will be skipped.\n')
        return(list(tc = tc, sra = sra))
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra2 = out[["sra"]]
    sra2[,'num_other_run_same_bp_curate_group'] = 0
    selected_curate_groups = selected_curate_groups[selected_curate_groups %in% unique(sra2[['curate_group']])]
    exclude_runs = c()
    for (sra_run in colnames(tc)) {
        is_sra = (sra2[['run']] == sra_run)
        my_curate_group = sra2[is_sra, "curate_group"]
        my_bioproject = sra2[is_sra, "bioproject"]
        is_not_my_bp = (sra2[['bioproject']] != my_bioproject)
        is_my_curate_group = (sra2[['curate_group']] == my_curate_group)
        run_other_bp = sra2[(is_not_my_bp | !is_my_curate_group), "run"]
        run_other_bp = run_other_bp[run_other_bp %in% colnames(tc)]
        tc_other_bp = tc[, run_other_bp]
        num_other_run_same_bp_curate_group = length(unique(sra2[(is_not_my_bp & is_my_curate_group), "bioproject"]))
        sra2[is_sra, "num_other_run_same_bp_curate_group"] = num_other_run_same_bp_curate_group
        sra2_other_bp <- sra2[sra2[['run']] %in% run_other_bp,]

        # If one curate_group is completely sourced from the same bioproject, we can't remove the whole bioproject for tc_ave_other_bp

        num_other_bp_same_curate_group = sum(sra2_other_bp[['curate_group']] == my_curate_group, na.rm=TRUE)
        if (num_other_bp_same_curate_group == 0) {
            tc_ave_other_bp = curate_group_mean(tc, sra2, selected_curate_groups)[['tc_ave']]
        } else {
            tc_ave_other_bp = curate_group_mean(tc_other_bp, sra2, selected_curate_groups)[['tc_ave']]
        }
        tc_ave = curate_group_mean(tc, sra2, selected_curate_groups)[['tc_ave']]
        coef = c()
        coef_other_bp = c()
        for (curate_group in selected_curate_groups) {
            tmp_coef = cor(tc[, sra_run], tc_ave[, curate_group], method=dist_method)
            if ((num_other_bp_same_curate_group == 0) & (curate_group == my_curate_group)) {
              tmp_coef_other_bp = cor(tc[, sra_run], tc_ave[, curate_group], method=dist_method)
            } else {
              tmp_coef_other_bp = cor(tc[, sra_run], tc_ave_other_bp[, curate_group], method=dist_method)
            }
            if (curate_group == my_curate_group) {
                tmp_coef = tmp_coef - min_dif
                tmp_coef_other_bp = tmp_coef_other_bp - min_dif
            }
            coef = c(coef, tmp_coef)
            coef_other_bp = c(coef_other_bp, tmp_coef_other_bp)
        }
        names(coef) = selected_curate_groups
        names(coef_other_bp) = selected_curate_groups
        if (max(coef, na.rm=TRUE) != coef[my_curate_group]) {
            cat('Registered as a candidate for exclusion. Better correlation to other categories:', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
        if (coef_other_bp[my_curate_group] < correlation_threshold) {
            cat('Registered as a candidate for exclusion. Low within-category correlation:', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
    }
    if (length(exclude_runs)) {
      if(one_out_per_iter == TRUE){
          cat("Excluding only one outlier per bioproject or same curate_group. \n")
          exclude_run_bps_and_curate_group = sra2[(sra2$run %in% exclude_runs), c("bioproject", "run", "curate_group")]
          first_bp_hit = exclude_run_bps_and_curate_group[match(unique(exclude_run_bps_and_curate_group$bioproject), exclude_run_bps_and_curate_group$bioproject),]
          first_same_curate_group_hit = exclude_run_bps_and_curate_group[match(unique(exclude_run_bps_and_curate_group$curate_group), exclude_run_bps_and_curate_group$curate_group),]
          # if a first_same_curate_group_hit is part of the same bioproject as the other removal candidates, ommit the same curate_group candidates
          if(any(first_same_curate_group_hit$bioproject %in% first_bp_hit$bioproject)){
              exclude_runs_tmp = c(first_bp_hit$run,first_same_curate_group_hit[!first_same_curate_group_hit$bioproject %in% first_bp_hit$bioproject]$run)
          }else{
              exclude_runs_tmp = c(first_bp_hit$run,first_same_curate_group_hit$run)
          }
          exclude_runs = unique(exclude_runs_tmp)
          # TODO This boolean vector should be all TRUE by definition. ???: exclude_run_bps_and_curate_group$run %in% exclude_runs
          exclude_bps = exclude_run_bps_and_curate_group[exclude_run_bps_and_curate_group$run %in% exclude_runs, "bioproject"]
      } else {
        exclude_run_bps = sra2[(sra2[['run']] %in% exclude_runs), c("bioproject", "run", "num_other_run_same_bp_curate_group")]
        exclude_bp_counts = data.frame(table(exclude_run_bps[['bioproject']]))
        exclude_run_bps = merge(exclude_run_bps, exclude_bp_counts, by.x = "bioproject", by.y = "Var1")
        exclude_run_bps = exclude_run_bps[order(exclude_run_bps[['num_other_run_same_bp_curate_group']], exclude_run_bps[['Freq']]),]
        rownames(exclude_run_bps) = 1:nrow(exclude_run_bps)
        min_other_run_same_bp_curate_group = exclude_run_bps[1, "num_other_run_same_bp_curate_group"]
        semimin_bp_count = exclude_run_bps[1, "Freq"]
        cat("minimum number of other BioProjects within curate_group:", min_other_run_same_bp_curate_group, "\n")
        cat("semi-minimum count of exclusion-candidate BioProjects:", semimin_bp_count, "\n")
        conditions = (exclude_run_bps[['Freq']] == semimin_bp_count)
        conditions = conditions & (exclude_run_bps[['num_other_run_same_bp_curate_group']] == min_other_run_same_bp_curate_group)
        exclude_bps = unique(exclude_run_bps[conditions, "bioproject"])
        exclude_runs = exclude_run_bps[(exclude_run_bps[['bioproject']] %in% exclude_bps), "run"]
      }
    }
    if (length(exclude_runs)) {
        cat('Partially removed BioProjects due to low within-curate_group correlation:', paste(exclude_bps, collapse=' '), '\n')
        cat('Removed Runs due to low within-curate_group correlation:', paste(exclude_runs, collapse=' '), '\n')
    }
    tc = tc[, !colnames(tc) %in% exclude_runs, drop=FALSE]
    sra[(sra[['run']] %in% exclude_runs), "exclusion"] = "low_within_curate_group_correlation"
    return(list(tc = tc, sra = sra))
}

batch_effect_subtraction = function(tc, sra, batch_effect_alg, transform_method, clip_negative) {
    if (ncol(tc)==1) {
        cat('Only 1 sample is available. Skipping batch effect subtraction.\n')
        return(list(tc = tc, sva = NULL))
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    out = sort_tc_and_metadata(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    out = remove_nonexpressed_gene(tc)
    tc = out[["tc_ex"]]
    tc_ne = out[["tc_ne"]]
    if (batch_effect_alg == "sva") {
        mod = try(model.matrix(~curate_group, data = sra))
        if ("try-error" %in% class(mod)) {
            return(list(tc = tc, sva = NULL))
        }
        mod0 = model.matrix(~1, data = sra)
        sva1 = try(sva(dat = as.matrix(tc), mod = mod, mod0 = mod0, B = 10))
        if ((class(sva1) != "try-error")) {
            cat("SVA correction was correctly performed.\n")
            tc = cleanY(y = tc, mod = mod, svs = sva1[['sv']])
            tc = rbind(tc, tc_ne)
        } else if ((class(sva1) == "try-error")) {
            cat("SVA correction failed.")
            # tc = rbind(tc, tc_ne) # Original tc shouldn't be returned as it is confusing. We may need a logic to retrieve the result of the previous round as final output.
        }
    } else if (batch_effect_alg == "combatseq") {
        bp_freq = as.data.frame(table(sra[,"bioproject"]))
        bp_freq_gt1 = bp_freq[bp_freq[,"Freq"]>1, "Var1"]
        bp_freq_eq1 = bp_freq[bp_freq[,"Freq"]==1, "Var1"]
        run_bp_freq_gt1 = sra[sra[,"bioproject"] %in% bp_freq_gt1, "run"]
        run_bp_freq_eq1 = sra[sra[,"bioproject"] %in% bp_freq_eq1, "run"]
        tc_combat = tc[, colnames(tc) %in% run_bp_freq_gt1]
        tcc_cn = colnames(tc_combat)
        batch = sra[(sra[,"bioproject"] %in% bp_freq_gt1),"bioproject"]
        group = sra[(sra[,"run"] %in% run_bp_freq_gt1),"curate_group"]
        tc_combat = try(ComBat_seq(as.matrix(tc_combat), batch=batch, group=group))
        if (class(tc_combat)[1] != "try-error") {
            cat("These runs are being removed, due to the bioproject only having 1 sample: \n")
            print(run_bp_freq_eq1)
            cat("Combatseq correction was correctly performed.\n")
            tc_combat = as.data.frame(tc_combat)
            colnames(tc_combat) = tcc_cn
            tc = rbind(tc_combat, tc_ne[, colnames(tc_combat)])
            sva1 = ''
        } else {
            stop("Combatseq correction failed.")
            # tc = rbind(tc, tc_ne) # Original tc shouldn't be returned as it is confusing. We may need a logic to retrieve the result of the previous round as final output.
            sva1 = ''
        }
    } else if (batch_effect_alg == "ruvseq") {
        x = as.factor(sra$curate_group)
        design = try(model.matrix(~curate_group, data = sra))
        if ("try-error" %in% class(design)) {
            return(list(tc = tc, sva = NULL))
        }
        y = DGEList(counts=as.matrix(tc+1), group=x)
        y = calcNormFactors(y, method="upperquartile")
        y = estimateGLMCommonDisp(y, design)
        y = estimateGLMTagwiseDisp(y, design)
        fit = glmFit(y, design)
        res = residuals(fit, type="deviance")
        seqUQ = betweenLaneNormalization(as.matrix(tc+1), which="upper", round = TRUE, offset = FALSE)
        controls = rep(TRUE,dim(as.matrix(tc+1))[1])
        batch_ruv_res = try(RUVr(seqUQ,controls,k=1,res)[[2]])
        if (class(batch_ruv_res)[1] != "try-error") {
            cat("RUVseq correction was correctly performed.\n")
            tc = rbind(batch_ruv_res, tc_ne)
            sva1 = ''
        } else {
            stop("RUVseq correction failed.")
            # tc = rbind(tc, tc_ne) # Original tc shouldn't be returned as it is confusing. We may need a logic to retrieve the result of the previous round as final output.
        }
    } else if (batch_effect_alg == "no") {
        cat("No batch effect correction was performed.\n")
        tc = rbind(tc, tc_ne)
        sva1 = ''
    } else {
        stop("Invalid batch effect correction algorithm.")
    }
    if (clip_negative) {
        if (endsWith(sub('-.*', '', transform_method), 'p1')) {
            is_negative = (tc < 0)
            txt = 'Number of negative values clipped to zero: %s/%s (%s x %s matrix)\n'
            val1 = format(sum(is_negative), big.mark=',', scientific=FALSE)
            val2 = format(prod(dim(tc)), big.mark=',', scientific=FALSE)
            val3 = format(dim(tc)[1], big.mark=',', scientific=FALSE)
            val4 = format(dim(tc)[2], big.mark=',', scientific=FALSE)
            cat(sprintf(txt, val1, val2, val3, val4))
            if (sum(is_negative) > 0) {
                tc[is_negative] = 0
            }
        } else {
            cat('`--clip_negative yes` is only applicable to `--norm log*p1-*`.\n')
        }
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
    curate_group_fac = factor(sra[, c("curate_group")])
    ann_label = data.frame(bioproject = bp_fac, curate_group = curate_group_fac)
    bp_col_uniq = unique(sra[order(sra[['bioproject']]), 'bp_color'])
    curate_group_col_uniq = unique(sra[order(sra[['curate_group']]), 'curate_group_color'])
    ann_color = list(bioproject = bp_col_uniq, curate_group = curate_group_col_uniq)
    breaks = c(0, seq(0.3, 1, 0.01))
    colnames(tc_dist_matrix) = sra[sra$run %in% colnames(tc_dist_matrix), 'run']
    head_half = substr(colnames(tc_dist_matrix), 1, 4)
    tail_half = substr(colnames(tc_dist_matrix), length(colnames(tc_dist_matrix))-3, length(colnames(tc_dist_matrix)))
    short_names = paste0(head_half, '..', tail_half)
    rownames(tc_dist_matrix) = ifelse(nchar(rownames(tc_dist_matrix))<=10, rownames(tc_dist_matrix), short_names)
    colnames(tc_dist_matrix) = ifelse(nchar(colnames(tc_dist_matrix))<=10, colnames(tc_dist_matrix), short_names)
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
    dend = as.dendrogram(hclust(tc_dist_dist))
    dend_colors = sra[order.dendrogram(dend),'curate_group_color']
    labels_colors(dend) = dend_colors
    dend_labels = sra[order.dendrogram(dend), 'run']
    dend = color_branches(dend, labels = dend_labels, col = dend_colors)
    dend = set(dend, "branches_lwd", 1)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    cex.xlab = min(fontsize, max(0.2, 0.5/log10(nrow(sra)), na.rm=TRUE), na.rm=TRUE)
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
      bg = sra[order.dendrogram(dend), 'curate_group_color'],
      fg = sra[order.dendrogram(dend), 'bp_color']
    )
}

draw_dendrogram_ggplot = function(sra, tc_dist_dist, fontsize = 7) {

  cg_col = unique(sra[,c('curate_group', 'curate_group_color')])
  bp_col = unique(sra[,c('bioproject', 'bp_color')])
  colnames(cg_col) = c('Group', 'Color')
  colnames(bp_col) = c('Group', 'Color')
  sra_colors = rbind(cg_col,bp_col)

  group_colors <-  data.table(Group=sra_colors$Group, Color=sra_colors$Color, key="Group")
  group_colors <- transpose(group_colors, make.names = "Group")
  colnames(tc_dist_dist)<-sra$run
  hc       <- hclust(tc_dist_dist)           # heirarchal clustering
  dendr    <- dendro_data(hc, type="rectangle") # convert for ggplot


  clust.df <- data.frame(label=sra$run, curate_group=factor(sra[sra$run %in% dendr$labels$label,'curate_group']), bioproject = factor(sra[sra$run %in% dendr$labels$label,'bioproject']))
  dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
  ggplot() +
    geom_segment(data=dendr$segments, aes(x=x, y=y, xend=xend, yend=yend), size = .8, show.legend = FALSE) +
    geom_segment(data=merge(dendr$segments[dendr$segments$yend==0,], dendr$labels[,c('label', 'curate_group', 'x')], by = 'x'), aes(x=x, y=y, xend=xend, yend=yend, color = curate_group), size = .8, show.legend = FALSE) +
    geom_text(data=dendr$labels, aes(x, y - .008, label=label, hjust=0, angle = 270, color = curate_group ), size=3, show.legend = FALSE) +
    scale_y_continuous(expand = c(.2, .1)) +
    geom_point(data = dendr$labels, aes(x,y, color = bioproject), size = 3, show.legend = FALSE) +
    geom_point(data = dendr$labels, aes(x,y, color = curate_group), size = 2, show.legend = FALSE) +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank()
    ) + scale_color_manual(values=group_colors)
}

draw_dendrogram_pvclust = function(sra, tc, nboot, pvclust_file, fontsize = 7) {
    dist_fun = function(x) {
        Dist(t(x), method = "pearson")
    }
    sp = sub(" ", "_", sra[['scientific_name']][1])
    if (file.exists(pvclust_file)) {
        if (file.info(pvclust_file)[['size']]) {
            cat("pvclust intermediate file found.\n")
            load(pvclust_file)
        }
    } else {
        cat("no pvclust intermediate file found. Start bootstrapping.\n")
        result = pvclust(tc, method.dist = dist_fun, method.hclust = "average", nboot = nboot, parallel = FALSE)  # UPGMA
        save(result, file = pvclust_file)
    }
    dend = as.dendrogram(result)
    dend_colors = sra[order.dendrogram(dend), 'curate_group_color']
    labels_colors(dend) = dend_colors
    dend_labels = sra[order.dendrogram(dend), 'run']
    dend = color_branches(dend, labels = dend_labels, col = dend_colors)
    dend = set(dend, "branches_lwd", 2)
    for (i in 1:ncol(tc)) {
        dend = dendrapply(dend, color_children2parent)
    }
    cex.xlab = min(0.2 + 1/log10(fontsize), na.rm=TRUE)
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
      bg = sra[order.dendrogram(dend), 'curate_group_color'],
      fg = sra[order.dendrogram(dend), 'bp_color']
    )
    text(result, print.num = FALSE, cex = 1, col.pv = "black")
}

draw_pca = function(sra, tc_dist_matrix, fontsize = 7) {
    pca = prcomp(tc_dist_matrix)
    xlabel = paste0("PC 1 (", round(summary(pca)[['importance']][2, 1] * 100, digits = 1), "%)")
    ylabel = paste0("PC 2 (", round(summary(pca)[['importance']][2, 2] * 100, digits = 1), "%)")
    plot(
      pca[['x']][, 1],
      pca[['x']][, 2],
      pch = 21,
      cex = 2,
      lwd = 1,
      bg = sra[['curate_group_color']],
      col = sra[['bp_color']],
      xlab = xlabel,
      ylab = ylabel,
      las = 1
    )
    # plot(pca$x[,1], pca$x[,2], pch=21, cex=2, lwd=2, bg=sra$curate_group_color, col=sra$bp_color, main=title,
    # xlab=xlabel, ylab=ylabel, las=1)
}

draw_mds = function(sra, tc_dist_dist, fontsize = 7) {
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
          bg = sra[['curate_group_color']],
          col = sra[['bp_color']],
          xlab = "MDS dimension 1",
          ylab = "MDS dimension 2",
          las = 1
        )
    }
}

draw_tsne = function(sra, tc, fontsize = 7) {
    perplexity = min(30, floor(nrow(sra)/4), na.rm=TRUE)
    try_out = tryCatch({
        Rtsne(
          as.matrix(t(tc)),
          theta = 0,
          check_duplicates = FALSE,
          verbose = FALSE,
          perplexity = perplexity,
          dims = 2
        )
    }, error = function(a) {
        return("t-SNE calculation failed.")
    })
    if (mode(try_out) == "character") {
        cat("t-SNE failed.\n")
        plot(c(0, 1), c(0, 1), ann=FALSE, bty="n", type="n", xaxt="n", yaxt="n")
        return(NULL)
    }
    out_tsne = try_out
    try_out = tryCatch({
        plot(
          out_tsne[['Y']][, 1],
          out_tsne[['Y']][, 2],
          pch = 21,
          cex = 2,
          lwd = 1,
          bg = sra[['curate_group_color']],
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
        out = tc_metadata_intersect(tc, sra)
        tc = out[["tc"]]
        sra = out[["sra"]]
        out = sort_tc_and_metadata(tc, sra)
        tc = out[["tc"]]
        sra = out[["sra"]]
        sra[['log10_total_spots']] = log10(sra[['total_spots']])
        sra[['log10_total_bases']] = log10(sra[['total_bases']])
        cols = c("curate_group","bioproject","lib_layout","lib_selection","instrument","mapping_rate",'log10_total_spots','log10_total_bases')
        label_cols = c("Group","BioProject","Library layout","Library selection","Instrument","Mapping rate",'Log10 total reads','Log10 total bases')
        if ('tmm_normalization_factor' %in% colnames(sra)) {
            cols = c(cols, 'tmm_normalization_factor')
            label_cols = c(label_cols, 'TMM normalization factor')
        }
        num_sv = sva_out[['n.sv']]
        if (num_sv==0) {
            cat('No surrogate variables found.\n')
            plot(c(0, 1), c(0, 1), ann=FALSE, bty="n", type="n", xaxt="n", yaxt="n")
            return(data.frame())
        }
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
    is_same_curate_group = outer(sra[['curate_group']], sra[['curate_group']], function(x, y) {
        x == y
    })
    plot(c(0.5, 4.5), c(0, 1), type = "n", xlab = "", ylab = "Pearson's correlation\ncoefficient", las = 1,
         xaxt = "n")
    boxplot(tc_dist_matrix[(!is_same_bp) & (!is_same_curate_group)], at = 1, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(is_same_bp) & (!is_same_curate_group)], at = 2, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(!is_same_bp) & (is_same_curate_group)], at = 3, add = TRUE, col = "gray", yaxt = "n")
    boxplot(tc_dist_matrix[(is_same_bp) & (is_same_curate_group)], at = 4, add = TRUE, col = "gray", yaxt = "n")
    labels = c("bw\nbw", "bw\nwi", "wi\nbw", "wi\nwi")
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = "Organ\nBioProject", padj = 0.5, hadj = 1, tick = FALSE)

}

draw_tau_histogram = function(tc, sra, selected_curate_groups, fontsize = 7, transform_method) {
    df_tau = curate_group2tau(curate_group_mean(tc, sra, selected_curate_groups)[['tc_ave']], rich.annotation = FALSE, transform_method)
    hist_out = hist(df_tau[['tau']], breaks = seq(0, 1, 0.05), las = 1, xlab = "Tau (expression specificity)",
                    ylab = "Gene count", main = "", col = "gray")
    num_noexp = sum(is.na(df_tau[['tau']]))
    num_all = nrow(df_tau)
    text_noexp = paste0("Excluded due to\nno expression:\n", num_noexp, "/", num_all, " genes")
    text(0, max(hist_out[['counts']], na.rm=TRUE) * 0.85, text_noexp, pos = 4)
}

draw_exp_level_histogram = function(tc, sra, selected_curate_groups, fontsize = 7, transform_method) {
    tc_curate_group = curate_group_mean(tc, sra, selected_curate_groups)[['tc_ave']]
    xmax = apply(tc_curate_group, 1, max)
    xmax[xmax < 0] = 0
    xmax[xmax > 15] = 15
    breaks = seq(0, 15, 1)
    hist_out = hist(xmax, breaks = breaks, las = 1, xlab = paste0("Max expression (",transform_method,")"), ylab = "Gene count",
                    main = "", col = "gray")
}

draw_legend = function(sra, new = TRUE, pos = "center", fontsize = 7, nlabel.in.col) {
    if (new) {
        plot.new()
    }
    curate_group_unique = unique(sra[['curate_group']])
    bp_unique = unique(sub(";.*", "", sra[['bioproject']]))
    curate_group_color_unique = unique(sra[['curate_group_color']])
    bp_color_unique = unique(sra[['bp_color']])
    ncol = ceiling((length(curate_group_unique) + length(bp_unique) + 2)/nlabel.in.col)
    legend_text = c("Organ", as.character(curate_group_unique), "", "BioProject", as.character(bp_unique))
    legend_color = c(rgb(1, 1, 1, 0), rep(rgb(1, 1, 1, 0), length(curate_group_color_unique)), rgb(1, 1, 1,
                                                                                             0), rgb(1, 1, 1, 0), bp_color_unique)
    legend_bg = c(rgb(1, 1, 1, 0), curate_group_color_unique, rgb(1, 1, 1, 0), rgb(1, 1, 1, 0), rep(rgb(1,
                                                                                                  1, 1, 0), length(bp_color_unique)))
    legend_font = c(2, rep(1, length(curate_group_color_unique)), 1, 2, rep(1, length(bp_color_unique)))
    legend(pos, legend = legend_text, pch = 21, lwd = 1, lty = 0, col = legend_color, pt.bg = legend_bg,
           text.font = legend_font, ncol = ncol, bty = "n")
}

save_plot = function(tc, sra, sva_out, dist_method, file, selected_curate_groups, fontsize = 7, transform_method, batch_effect_alg) {
    if (ncol(tc)==1) {
        cat('Only 1 sample is available. Skipping the plot.\n')
        return()
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    sra = add_color_to_metadata(sra, selected_curate_groups)
    out = sort_tc_and_metadata(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    pdf(file.path(dir_pdf, paste0(file, ".pdf")), height = 8, width = 7.2, fonts = "Helvetica", pointsize = fontsize)
    layout_matrix = matrix(c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1,
                             2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1,
                             1, 1, 1, 1, 1, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3,
                             3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9,
                             9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10,
                             10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), 14, 12,
                           byrow = TRUE)
    layout(layout_matrix)
    colnames(tc)<-sra$run
    ##
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
    draw_exp_level_histogram(tc, sra, selected_curate_groups, fontsize, transform_method)
    par(mar = c(4, 4, 1, 1))
    draw_tau_histogram(tc, sra, selected_curate_groups, fontsize, transform_method)
    par(mar = rep(0.1, 4))
    if (batch_effect_alg == 'sva') {
        df_r2 = draw_sva_summary(sva_out, tc, sra, fontsize)
        if (!all(is.na(df_r2))) {
            write.table(df_r2, file.path(dir_tsv, paste0(file, ".r2.tsv")), sep = "\t", row.names = FALSE)
        }
    }
    par(mar = rep(0.1, 4))
    draw_legend(sra, new = TRUE, pos = "center", fontsize = fontsize, nlabel.in.col = 8)
    graphics.off()
}


save_correlation = function( tc, sra, dist_method) {
    is_same_bp = outer(sra[['bioproject']], sra[['bioproject']], function(x, y) {
        x == y
    })
    is_same_curate_group = outer(sra[['curate_group']], sra[['curate_group']], function(x, y) {
        x == y
    })
    counter <- -1

  # Internal function to generate row names
    generate_row_names <- function() {
        counter <<- counter + 1
        return(paste0("round_", counter))
    }


  tc_dist_matrix = cor(tc, method = dist_method)
  tc_dist_matrix[is.na(tc_dist_matrix)] = 0

  # bw = between; wi = within; order: Bioproject -> Curate Group; tc_dist_bwwi: between BPs, within CGs
  tc_dist_bwbw = tc_dist_matrix[(!is_same_bp) & (!is_same_curate_group)]
  bwbw_mea = mean(tc_dist_bwbw)
  bwbw_med = median(tc_dist_bwbw)
  bwbw_var = var(tc_dist_bwbw)

  tc_dist_wibw = tc_dist_matrix[(is_same_bp) & (!is_same_curate_group)]
  wibw_mea = mean(tc_dist_wibw)
  wibw_med = median(tc_dist_wibw)
  wibw_var = var(tc_dist_wibw)

  tc_dist_bwwi = tc_dist_matrix[(!is_same_bp) & (is_same_curate_group)]
  bwwi_mea = mean(tc_dist_bwwi)
  bwwi_med = median(tc_dist_bwwi)
  bwwi_var = var(tc_dist_bwwi)

  tc_dist_wiwi = tc_dist_matrix[(is_same_bp) & (is_same_curate_group)]
  wiwi_mea = mean(tc_dist_wiwi)
  wiwi_med = median(tc_dist_wiwi)
  wiwi_var = var(tc_dist_wiwi)

  tc_dist_stats = c(bwbw_mea, bwbw_med, bwbw_var, wibw_mea, wibw_med, wibw_var, bwwi_mea, bwwi_med, bwwi_var, wiwi_mea, wiwi_med, wiwi_var)
  cat(tc_dist_stats)
  # Check if dataframe exists in environment, if not create it
  if (!exists("correlation_statistics", envir = .GlobalEnv)) {
    correlation_statistics <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean","bwbw_median","bwbw_variance", "wibw_mean","wibw_median","wibw_variance", "bwwi_mean","bwwi_median","bwwi_variance", "wiwi_mean","wiwi_median","wiwi_variance"))))
    rownames(correlation_statistics) <- generate_row_names()
    assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
  } else {
    correlation_statistics <- get("correlation_statistics", envir = .GlobalEnv)
    new_row <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean","bwbw_median","bwbw_variance", "wibw_mean","wibw_median","wibw_variance", "bwwi_mean","bwwi_median","bwwi_variance", "wiwi_mean","wiwi_median","wiwi_variance"))))
    rownames(new_row) <- generate_row_names()
    correlation_statistics <- rbind(correlation_statistics, new_row)
    assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
  }

  return(correlation_statistics)


}


transform_raw_to_fpkm = function(counts, effective_lengths, sra) {
    if ('tmm_library_size' %in% colnames(sra)) {
        cat('FPKM transformation with the original library sizes from amalgkit cstmm output.\n')
        cat('If --input_dir is specified with amalgkit cstmm output, resultant values will be TMM-FPKM.\n')
        cat('If --input_dir is specified with amalgkit merge output, resultant values will be non-TMM-corrected FPKM.\n')
        tmp = sra
        rownames(tmp) = tmp[['run']]
        library_sizes = tmp[colnames(counts),'tmm_library_size']
    } else {
        cat('FPKM transformation with the library sizes in the input files.\n')
        cat('Irrespective of --input_dir, resultant values will be non-TMM-corrected FPKM.\n')
        library_sizes = colSums(counts)
    }
    res = counts / effective_lengths / library_sizes * 1e+09 # 1e+09 = kb * million reads
    return(as.data.frame(res))
}

transform_raw_to_tpm = function(counts, effective_lengths) {
    x <- counts / effective_lengths
    res = t(t(x) * 1e+06 / colSums(x))
    return(res)
}

apply_transformation_logic = function(tc, tc_eff_length, transform_method, batch_effect_alg,
                                      step=c('before_batch','before_batch_plot','after_batch'), sra) {
    if (batch_effect_alg=='no') {
        if (step=='before_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step=='before_batch_plot') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step=='after_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        }
    } else if (batch_effect_alg=='sva') {
        if (step=='before_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step=='before_batch_plot') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step=='after_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        }
    } else if (batch_effect_alg=='ruvseq') {
        if (step=='before_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step=='before_batch_plot') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step=='after_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        }
    } else if (batch_effect_alg=='combatseq') {
        if (step=='before_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step=='before_batch_plot') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step=='after_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        }
    }
    if (bool_fpkm_tpm == TRUE) {
        if (grepl('fpkm', transform_method) ) {
            cat('Applying FPKM transformation.\n')
            tc = transform_raw_to_fpkm(tc, tc_eff_length[, colnames(tc)], sra)
        } else if (grepl('tpm', transform_method )) {
            cat('Applying TPM transformation.\n')
            tc = transform_raw_to_tpm(tc, tc_eff_length[, colnames(tc)])
        } else {
            cat('Applying neither FPKM nor TPM transformation.\n')
        }
    }
    if (bool_log == TRUE) {
        if (grepl('logn-', transform_method)) {
            cat('Applying log_n(x) normalization.\n')
            tc = log(tc)
        } else if (grepl('log2-', transform_method)) {
            cat('Applying log_2(x) normalization.\n')
            tc = log2(tc)
        } else if (grepl('lognp1-', transform_method)) {
            cat('Applying log_n(x+1) normalization.\n')
            tc = log(tc + 1)
        } else if (grepl('log2p1-', transform_method)) {
            cat('Applying log_2(x+1) normalization.\n')
            tc = log2(tc + 1)
        } else {
            cat('Applying no log normalization.\n')
        }
    }
    return(tc)
}

standardize_metadata_all = function(sra_all) {
    for (col in c('instrument','bioproject')) {
        is_missing = (sra_all[,col] == "")|(is.na(sra_all[,col]))
        sra_all[is_missing, col] = "not_provided"
    }
    return(sra_all)
}

get_species_metadata = function(sra_all, scientific_name, selected_curate_groups) {
    is_sp = (sra_all[,'scientific_name'] == scientific_name)
    is_curate_group = (sra_all[,'curate_group'] %in% selected_curate_groups)
    cat('Number of SRA runs for this species:', sum(is_sp), '\n')
    cat('Number of SRA runs for selected curate groups:', sum(is_curate_group), '\n')
    sra = sra_all[(is_sp & is_curate_group),]
    conditions = (sra[['exclusion']] == "no") & (!sra[['run']] %in% colnames(tc))
    if (any(conditions)) {
        cat("Failed quantification:", sra[conditions, "run"], "\n")
        sra[conditions, "exclusion"] = "failed_quantification"
    }
    return(sra)
}

exclude_inappropriate_sample_from_tc = function(tc, sra) {
    is_not_excluded = (sra[['exclusion']]=='no')
    cat('Number of non-excluded SRA runs (exclusion=="no"):', sum(is_not_excluded), '\n')
    tc = tc[,sra[is_not_excluded,'run'], drop=FALSE]
    return(tc)
}

exclude_inappropriate_sample_from_eff_length = function(tc_eff_length, tc) {
    row.names(tc_eff_length) = tc_eff_length[, 1]
    tc_eff_length = tc_eff_length[, colnames(tc)]
    return(tc_eff_length)
}

########cd############END OF FUNCTION DECLARATION####################################################

################################################# START OF SVA CORRECTION ########################################################


fontsize = 7

tc = read.table(est_counts_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, row.names = 1, check.names=FALSE)
tc_eff_length = read.table(eff_length_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, check.names=FALSE)

sra_all = read.table(metadata_path, sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "", stringsAsFactors = FALSE, check.names=FALSE)
sra_all = standardize_metadata_all(sra_all)

scientific_name = unique(sra_all[(sra_all[['run']] %in% colnames(tc)), "scientific_name"])
dir_curate = file.path(out_dir, 'curate')
dir_pdf = file.path(dir_curate, sub(" ", "_", scientific_name), 'plots')
dir.create(dir_pdf, showWarnings=FALSE, recursive=TRUE)
dir_rdata = file.path(dir_curate, sub(" ", "_", scientific_name), 'rdata')
dir.create(dir_rdata, showWarnings=FALSE, recursive=TRUE)
dir_tsv = file.path(dir_curate, sub(" ", "_", scientific_name), 'tables')
dir.create(dir_tsv, showWarnings=FALSE, recursive=TRUE)
setwd(dir_curate)
cat(log_prefix, "Working at:", getwd(), "\n")

sra = get_species_metadata(sra_all, scientific_name, selected_curate_groups)
tc = exclude_inappropriate_sample_from_tc(tc, sra)
out = sort_tc_and_metadata(tc, sra) ; tc = out[["tc"]] ; sra = out[["sra"]]
tc_eff_length = exclude_inappropriate_sample_from_eff_length(tc_eff_length, tc)
tc = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step='before_batch', sra=sra)
tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step='before_batch_plot', sra=sra)
is_input_zero = data.frame(tc_tmp==0, check.names=FALSE)
file_name = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.tc.tsv"))
write_table_with_index_name(df=tc_tmp, file_path=file_name, index_name='target_id')
out = curate_group_mean(tc_tmp, sra, selected_curate_groups)
tc_curate_group_uncorrected = out[['tc_ave']]
selected_curate_groups = out[['selected_curate_groups']]
file_name = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.curate_group.mean.tsv"))
write_table_with_index_name(df=tc_curate_group_uncorrected, file_path=file_name, index_name='target_id')


if (skip_curation_flag == TRUE) {
    cat("No curation requested, finishing early.\n")
    cat("Files created: \n")
    cat(file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.tc.tsv")) , "\n")
    cat(file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".uncorrected.curate_group.mean.tsv")), "\n")
    cat("Transformation applied: ", transform_method , "\n")
    cat(log_prefix, "Completed.\n")
    quit(save='no', status=0)
}

cat("Removing samples with mapping rate of 0.\n")
round = 0
sva_out = NULL
tc_batch_corrected = NULL
out = check_mapping_rate(tc, sra, 0)
tc = out[["tc"]]
sra = out[["sra"]]
tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step='before_batch_plot', sra=sra)
save_plot(tc_tmp, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".original"),
          selected_curate_groups, fontsize, transform_method, batch_effect_alg)
save_correlation(tc_tmp, sra, dist_method)
out = batch_effect_subtraction(tc, sra, batch_effect_alg, transform_method, clip_negative)
tc_batch_corrected = out[["tc"]]
sva_out = out[["sva"]]
if (!is.null(sva_out)) {
    file_name = paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg, ".", round, ".RData")
    save(sva_out, file = file.path(dir_rdata, file_name))
}
tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step='after_batch', sra=sra)
save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".original", ".", batch_effect_alg),
          selected_curate_groups, fontsize, transform_method, batch_effect_alg)

cat("Removing samples with the mapping rate smaller than", mapping_rate_cutoff, "\n")
round = 1
sva_out = NULL
tc_batch_corrected = NULL
out = check_mapping_rate(tc, sra, mapping_rate_cutoff)
tc = out[["tc"]]
sra = out[["sra"]]

tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step='before_batch_plot', sra=sra)
save_plot(tc_tmp, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".mapping_cutoff"),
          selected_curate_groups, fontsize, transform_method, batch_effect_alg)
out = batch_effect_subtraction(tc, sra, batch_effect_alg, transform_method, clip_negative)
tc_batch_corrected = out[["tc"]]
sva_out = out[["sva"]]
if (!is.null(sva_out)) {
    save(sva_out, file=file.path(dir_rdata, paste0(sub(" ", "_", scientific_name),".", batch_effect_alg,".", round, ".RData")))
}
tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step='after_batch', sra=sra)
save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, paste0(sub(" ", "_", scientific_name), ".", round, ".mapping_cutoff", ".", batch_effect_alg),
          selected_curate_groups, fontsize, transform_method, batch_effect_alg)
save_correlation(tc_batch_corrected_tmp, sra, dist_method)

round = 2
end_flag = 0
while (end_flag == 0) {
    cat("Iteratively checking within-curate_group correlation, round:", round, "\n")
    tc_cwtc = NULL
    num_run_before = sum(sra[['exclusion']] == "no")
    out = check_within_curate_group_correlation(tc, sra, dist_method, min_dif, selected_curate_groups, one_outlier_per_iteration, correlation_threshold)
    tc_cwtc = out[["tc"]]
    sra = out[["sra"]]
    num_run_after = sum(sra[['exclusion']] == "no")
    if ((num_run_before == num_run_after) | (plot_intermediate)) {
        sva_out = NULL
        tc_batch_corrected = NULL
        out = batch_effect_subtraction(tc[,colnames(tc_cwtc)], sra, batch_effect_alg, transform_method, clip_negative)
        tc_batch_corrected = out[["tc"]]
        sva_out = out[["sva"]]
        if (!is.null(sva_out)) {
            save(sva_out, file=file.path(dir_rdata, paste0(sub(" ", "_", scientific_name),".", batch_effect_alg, ".", round, ".RData")))
        }
        tc_cwtc_tmp = apply_transformation_logic(tc_cwtc, tc_eff_length, transform_method, batch_effect_alg, step='before_batch_plot', sra=sra)
        save_plot(tc_cwtc_tmp, sra, NULL, dist_method, paste0(sub(" ", "_", scientific_name), ".", round,
                                                          ".correlation_cutoff"), selected_curate_groups, fontsize, transform_method, batch_effect_alg)
        tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step='after_batch', sra=sra)
        file_base = paste0(sub(" ", "_", scientific_name), ".", round, ".correlation_cutoff", ".", batch_effect_alg)
        save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, file_base, selected_curate_groups, fontsize, transform_method, batch_effect_alg)
        save_correlation(tc_batch_corrected_tmp, sra, dist_method)

    }
    cat("Round:", round, ": # before =", num_run_before, ": # after =", num_run_after, "\n\n")
    if (num_run_before == num_run_after) {
        end_flag = 1
    }
    tc = tc_cwtc
    round = round + 1
}

cat("Finished checking within-curate_group correlation.\n")
if (batch_effect_alg != 'sva') {
    cat("Batch-effect removal algorithm is: ",batch_effect_alg, ". Applying transformation on final batch-removed counts.\n")
    tc_batch_corrected = tc_batch_corrected_tmp
}
if (maintain_zero) {
    cat('Any zero expression levels in the input will remain as zero-values in the output tables.\n')
    tc_batch_corrected = tc_batch_corrected[order(rownames(tc_batch_corrected)),]
    is_input_zero = is_input_zero[rownames(tc_batch_corrected),]
    stopifnot(all(rownames(is_input_zero)==rownames(tc_batch_corrected)))
    for (col in colnames(tc_batch_corrected)) {
        tc_batch_corrected[is_input_zero[[col]],col] = 0
    }
}
cat("Writing summary files for", scientific_name, "\n")
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".metadata.tsv"))
write.table(sra[,colnames(sra)!='index'], file=file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg , ".tc.tsv"))
write_table_with_index_name(df=tc_batch_corrected, file_path=file, index_name='target_id')
out = curate_group_mean(tc_batch_corrected, sra, selected_curate_groups)
tc_curate_group = out[['tc_ave']]
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg , ".curate_group.mean.tsv"))
write_table_with_index_name(df=tc_curate_group, file_path=file, index_name='target_id')
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg , ".correlation_statistics.tsv"))
write.table(correlation_statistics, file = file, row.names=TRUE, sep="\t", quotes=FALSE)
tc_tau = curate_group2tau(tc_curate_group, rich.annotation = TRUE, transform_method)
file = file.path(dir_tsv, paste0(sub(" ", "_", scientific_name), ".", batch_effect_alg , ".tau.tsv"))
write_table_with_index_name(df=tc_tau, file_path=file, index_name='target_id')
cat(log_prefix, "Completed.\n")
quit(save='no', status=0)