#!/usr/bin/env Rscript

detected_cores = tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
if (is.na(detected_cores) && (is.null(getOption("cores")) || is.na(getOption("cores")))) {
    options(cores = 1L)
}
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
    selected_sample_groups = c("root", "flower", "leaf")
    sample_group_colors = 'DEFAULT'
    transform_method = "log2p1-fpkm"
    one_outlier_per_iteration = as.logical(0)
    correlation_threshold = 0.3
    batch_effect_alg = 'sva'
    dist_method = "pearson"
    clip_negative = as.logical(1)
    maintain_zero = as.logical(1)
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
    skip_curation_flag = FALSE
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
    selected_sample_groups = strsplit(args[9], "\\|")[[1]]
    sample_group_colors = strsplit(args[10], ",")[[1]]
    transform_method = args[11]
    one_outlier_per_iteration = as.integer(args[12])
    correlation_threshold = as.numeric(args[13])
    batch_effect_alg = args[14]
    clip_negative = as.logical(as.integer(args[15]))
    maintain_zero = as.logical(as.integer(args[16]))
    r_util_path = file.path(args[17])
    skip_curation_flag = as.logical(as.integer(args[18]))
}
cat('est_counts_path:', est_counts_path, "\n")
cat('metadata_path:', metadata_path, "\n")
cat('out_dir:', out_dir, "\n")
cat('eff_length_path:', eff_length_path, "\n")
cat('dist_method:', dist_method, "\n")
cat('mapping_rate_cutoff:', mapping_rate_cutoff, "\n")
cat('min_dif:', min_dif, "\n")
cat('plot_intermediate:', plot_intermediate, "\n")
cat('selected_sample_groups:', selected_sample_groups, "\n")
cat('selected_sample_group_colors:', sample_group_colors, "\n")
cat('transform_method:', transform_method, "\n")
cat('one_outlier_per_iteration:', one_outlier_per_iteration, "\n")
cat('correlation_threshold:', correlation_threshold, "\n")
cat('batch_effect_alg:', batch_effect_alg, "\n")
cat('clip_negative:', clip_negative, "\n")
cat('maintain_zero:', maintain_zero, "\n")
cat('r_util_path:', r_util_path, "\n")
cat('skip_curation_flag:', skip_curation_flag, "\n")

source(r_util_path)

CURATE_FONT_SIZE_PT = 8
CURATE_FONT_FAMILY = 'Helvetica'

normalize_exclusion_values = function(exclusion_values) {
    normalized = as.character(exclusion_values)
    normalized[is.na(exclusion_values)] = NA_character_
    normalized = trimws(normalized)
    normalized = tolower(normalized)
    return(normalized)
}

is_non_excluded_flag = function(exclusion_values) {
    exclusion_norm = normalize_exclusion_values(exclusion_values)
    (!is.na(exclusion_norm)) & (exclusion_norm == 'no')
}

resolve_curate_fontsize = function(fontsize = CURATE_FONT_SIZE_PT) {
    return(CURATE_FONT_SIZE_PT)
}

if (batch_effect_alg == "ruvseq") {
    suppressWarnings(suppressPackageStartupMessages(library(RUVSeq, quietly = TRUE)))
} else {
    suppressWarnings(suppressPackageStartupMessages(library(sva, quietly = TRUE)))
}

tc_metadata_intersect = function(tc, sra) {
    sra_run = sra[['run']]
    tc = tc[, colnames(tc) %in% sra_run, drop = FALSE]
    sra = sra[sra[['run']] %in% colnames(tc),]
    return(list(tc = tc, sra = sra))
}

remove_nonexpressed_gene = function(tc) {
    gene_sum = apply(tc, 1, sum)
    tc_ex = tc[gene_sum > 0,]
    tc_ne = tc[gene_sum == 0,]
    return(list(tc_ex = tc_ex, tc_ne = tc_ne))
}

add_color_to_curate_metadata = function(sra, selected_sample_groups, sample_group_colors) {
    if (nrow(sra) == 0) {
        sra[['bp_color']] = character(0)
        sra[['sp_color']] = character(0)
        sra[['sample_group_color']] = character(0)
        return(sra)
    }
    sra = sra[, (!colnames(sra) %in% c("bp_color", "sp_color", "sample_group_color"))]
    if ('bioproject' %in% colnames(sra)) {
        bioproject = as.character(sra[['bioproject']])
        bioproject_u = sort(unique(bioproject))
    } else {
        bioproject_u = rep('PLACEHOLDER', nrow(sra))
    }
    scientific_name = as.character(sra[['scientific_name']])
    scientific_name_u = sort(unique(scientific_name))
    sample_group = as.character(sra[['sample_group']])
    sample_group_u = selected_sample_groups
    is_default_palette = (length(sample_group_colors) == 1) && (sample_group_colors == "DEFAULT")
    if ((!is_default_palette) && any(is.na(sample_group_colors))) {
        warning("Detected NA in sample_group_colors; falling back to DEFAULT colors.")
        sample_group_colors = "DEFAULT"
        is_default_palette = TRUE
    }
    if (is_default_palette) {
        if (length(selected_sample_groups) <= 8) {
            sample_group_color = brewer.pal(8, "Dark2")
            sample_group_color = sample_group_color[seq_len(length(selected_sample_groups))]
            bp_color = rainbow_hcl(length(bioproject_u), c = 50)
            sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
        } else if (length(selected_sample_groups) <= 12) {
            sample_group_color = brewer.pal(length(selected_sample_groups), "Paired")
            bp_color = rainbow_hcl(length(bioproject_u), c = 50)
            sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
        } else {
            sample_group_color = rainbow_hcl(length(selected_sample_groups), c = 100)
            bp_color = rainbow_hcl(length(bioproject_u), c = 50)
            sp_color = rainbow_hcl(length(scientific_name_u), c = 150)
        }
    } else {
        if (length(sample_group_colors) != length(selected_sample_groups)) {
            stop("Length of sample_group_colors must match length of selected_sample_groups")
        }
        sample_group_color = sample_group_colors
        bp_color = rainbow_hcl(length(bioproject_u), c = 50)
        sp_color = rainbow_hcl(length(scientific_name_u), c = 100)
    }
    df_sample_group = data.frame(sample_group = sample_group_u, sample_group_color = sample_group_color[seq_len(length(sample_group_u))], stringsAsFactors = FALSE)
    df_bp = data.frame(bioproject = bioproject_u, bp_color = bp_color[seq_len(length(bioproject_u))], stringsAsFactors = FALSE)
    df_sp = data.frame(scientific_name = scientific_name_u, sp_color = sp_color[seq_len(length(scientific_name_u))], stringsAsFactors = FALSE)
    sra = merge(sra, df_bp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_sp, sort = FALSE, all.y = FALSE)
    sra = merge(sra, df_sample_group, sort = FALSE, all.y = FALSE)
    return(sra)
}


sort_tc_and_metadata = function(tc, sra, sort_columns = c("sample_group", "scientific_name", "bioproject")) {
    for (column in rev(sort_columns)) {
        if (!column %in% colnames(sra)) {
            next
        }
        sra = sra[order(sra[[column]]),]
    }
    sra_intersection = sra[(sra[['run']] %in% colnames(tc)), 'run']
    tc = tc[, sra_intersection, drop = FALSE]
    return(list(tc = tc, sra = sra))
}

cleanY = function(y, mod, svs) {
    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    P = ncol(mod)
    return(y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P),]))
}

sample_group_mean = function(tc, sra, selected_sample_groups = NA, balance_bp = FALSE) {
    # if the data are SVA-corrected, balance_bp would not be necessary because project-specfic effects
    # were removed in SVA already.
    if ((ncol(tc) == 0) || (nrow(tc) == 0)) {
        tc_ave = data.frame(matrix(numeric(0), nrow = nrow(tc), ncol = 0))
        rownames(tc_ave) = rownames(tc)
        return(list(tc_ave = tc_ave, selected_sample_groups = selected_sample_groups))
    }
    if (all(is.na(selected_sample_groups))) {
        sp_sample_groups = unique(sra[['sample_group']])
    } else {
        sp_sample_groups = selected_sample_groups[selected_sample_groups %in% unique(sra[['sample_group']])]
    }
    tc_ave = data.frame(matrix(rep(NA, length(sp_sample_groups) * nrow(tc)), nrow = nrow(tc)))
    colnames(tc_ave) = sp_sample_groups
    rownames(tc_ave) = rownames(tc)
    sra_sample_group = sra[['sample_group']]
    sra_run = sra[['run']]
    sra_exclusion = sra[['exclusion']]
    sra_is_non_excluded = is_non_excluded_flag(sra_exclusion)
    tc_colnames = colnames(tc)
    is_run_in_tc = (sra_run %in% tc_colnames)
    if (balance_bp) {
        sra_bioproject = sra[['bioproject']]
        is_no_exclusion = sra_is_non_excluded
    }
    for (sample_group in sp_sample_groups) {
        is_sample_group = (sra_sample_group == sample_group)
        exclusion_sample_group = sra_is_non_excluded[is_sample_group]
        run_sample_group = sra_run[is_sample_group & is_run_in_tc]
        if (all(!exclusion_sample_group)) {
            warning_message = paste0('All samples of sample_group ', sample_group, ' are marked for exclusion. This sample_group will be omitted from further analysis.')
            selected_sample_groups = selected_sample_groups[!(selected_sample_groups %in% sample_group)]
            tc_ave = tc_ave[, !names(tc_ave) %in% c(sample_group)]
            warning(warning_message)
            next
        }
        if (length(run_sample_group) == 1) {
            exp_sample_group = tc[, run_sample_group]
        } else {
            if (balance_bp) {
                bps = unique(sra_bioproject[is_run_in_tc & is_sample_group & is_no_exclusion])
                df_tmp = data.frame(matrix(rep(NA, nrow(tc) * length(bps)), nrow = nrow(tc), ncol = length(bps)))
                colnames(df_tmp) = bps
                for (bp in bps) {
                    is_bp = (sra_bioproject == bp)
                    sra_ids = sra_run[is_bp & is_sample_group & is_no_exclusion]
                    tc_bp = tc[, sra_ids]
                    if (class(tc_bp) == "numeric") {
                        df_tmp[bp] = tc_bp
                    } else {
                        df_tmp[bp] = rowMeans(tc_bp)
                    }
                }
                exp_sample_group = rowMeans(df_tmp)
            } else {
                exp_sample_group = rowMeans(tc[, run_sample_group])
            }
        }
        tc_ave[, sample_group] = exp_sample_group
    }
    return(list(tc_ave = tc_ave, selected_sample_groups = selected_sample_groups))
}

row_tau = function(row) {
    is_nonzero = row > 0
    if (sum(is_nonzero) > 0) {
        exp_order = order(row[is_nonzero], decreasing = TRUE)
        sample_group_ordered = colnames(tc_sample_group)[is_nonzero][exp_order]
        highest = sample_group_ordered[1]
        order = paste(sample_group_ordered, collapse = "|")
        return(c(highest, order))
    } else {
        return(c(NA, NA))  # Return NA values if no nonzero elements
    }
}

sample_group2tau = function(tc_sample_group, rich.annotation = TRUE, transform_method) {
    if (rich.annotation) {
        cols = c("tau", "highest", "order")
    } else {
        cols = c("tau")
    }
    if ((nrow(tc_sample_group) == 0) || (ncol(tc_sample_group) == 0)) {
        df_tau = data.frame(matrix(nrow = 0, ncol = length(cols)))
        colnames(df_tau) = cols
        return(df_tau)
    }
    df_tau = data.frame(matrix(rep(NA, length(cols) * nrow(tc_sample_group)), nrow = nrow(tc_sample_group)))
    colnames(df_tau) = cols
    rownames(df_tau) = rownames(tc_sample_group)
    if (grepl('logn-', transform_method)) {
        tc_sample_group = exp(tc_sample_group)
    } else if (grepl('log2-', transform_method)) {
        tc_sample_group = 2**tc_sample_group
    } else if (grepl('lognp1-', transform_method)) {
        tc_sample_group = exp(tc_sample_group) - 1
    } else if (grepl('log2p1-', transform_method)) {
        tc_sample_group = 2**tc_sample_group - 1
    }
    tc_sample_group[tc_sample_group < 0] = 0
    xmax = apply(tc_sample_group, 1, function(x) { max(x, na.rm = TRUE) })
    df_tau[, 'tau'] = apply((1 - (tc_sample_group / xmax)) / (ncol(tc_sample_group) - 1), 1, sum)
    if (rich.annotation) {
        tc_sample_group[is.na(tc_sample_group)] = 0
        results = t(apply(tc_sample_group, 1, row_tau))
        df_tau[, c("highest", "order")] = results
    }
    return(df_tau)
}

check_mapping_rate = function(tc, sra, mapping_rate_cutoff) {
    if ('mapping_rate' %in% colnames(sra)) {
        cat(paste0('Mapping rate cutoff: ', mapping_rate_cutoff * 100, '%\n'))
        is_mapping_good = (sra[['mapping_rate']] > mapping_rate_cutoff * 100)
        is_mapping_good[is.na(is_mapping_good)] = TRUE
        if (any(!is_mapping_good)) {
            cat("Removed due to low mapping rate:\n")
            df_tmp = sra[!is_mapping_good,]
            for (i in rownames(df_tmp)) {
                sra_id = df_tmp[i, 'run']
                mapping_rate = df_tmp[i, 'mapping_rate']
                cat(paste0(sra_id, ': mapping rate = ', mapping_rate, '%\n'))
            }
            tc = tc[, colnames(tc) %in% sra[is_mapping_good, "run"], drop = FALSE]
        } else {
            cat("No entry removed due to low mapping rate.\n")
        }
        sra[!is_mapping_good, "exclusion"] = "low_mapping_rate"
    } else {
        cat('Mapping rate cutoff will not be applied.\n')
    }
    return(list(tc = tc, sra = sra))
}

check_within_sample_group_correlation = function(tc, sra, dist_method, min_dif, selected_sample_groups, one_out_per_iter = TRUE, correlation_threshold) {
    if (length(selected_sample_groups) == 1) {
        cat('Only one sample_group category is available. Outlier removal will be skipped.\n')
        return(list(tc = tc, sra = sra))
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra2 = out[["sra"]]
    if ((ncol(tc) == 0) || (nrow(sra2) == 0)) {
        cat('No sample is available. Outlier removal will be skipped.\n')
        return(list(tc = tc, sra = sra))
    }
    sra2[, 'num_other_run_same_bp_sample_group'] = 0
    selected_sample_groups = selected_sample_groups[selected_sample_groups %in% unique(sra2[['sample_group']])]
    if (length(selected_sample_groups) <= 1) {
        cat('Less than two sample_group categories are available after filtering. Outlier removal will be skipped.\n')
        return(list(tc = tc, sra = sra))
    }
    # This value is invariant inside the loop: tc, sra2 columns used by sample_group_mean(), and selected_sample_groups are unchanged.
    tc_ave_all = sample_group_mean(tc, sra2, selected_sample_groups)[['tc_ave']]
    coef_matrix_all = cor(tc, tc_ave_all, method = dist_method)
    if (is.null(dim(coef_matrix_all))) {
        coef_matrix_all = matrix(coef_matrix_all, nrow = 1, ncol = 1)
        rownames(coef_matrix_all) = colnames(tc)
        colnames(coef_matrix_all) = colnames(tc_ave_all)
    }
    sra2_run = sra2[['run']]
    sra2_sample_group = sra2[['sample_group']]
    sra2_bioproject = sra2[['bioproject']]
    tc_runs = colnames(tc)
    run_to_index = setNames(seq_len(nrow(sra2)), sra2_run)
    run_context_cache = list()
    tc_ave_other_bp_cache = list()
    coef_other_bp_cache = list()
    exclude_runs = c()
    for (sra_run in colnames(tc)) {
        is_sra = run_to_index[[sra_run]]
        my_sample_group = sra2_sample_group[[is_sra]]
        my_bioproject = sra2_bioproject[[is_sra]]
        run_context_key = paste0(my_sample_group, '\t', my_bioproject)
        if (!is.null(run_context_cache[[run_context_key]])) {
            run_other_bp = run_context_cache[[run_context_key]][['run_other_bp']]
            num_other_run_same_bp_sample_group = run_context_cache[[run_context_key]][['num_other_run_same_bp_sample_group']]
            num_other_bp_same_sample_group = run_context_cache[[run_context_key]][['num_other_bp_same_sample_group']]
            run_context_runs = run_context_cache[[run_context_key]][['run_context_runs']]
        } else {
            is_not_my_bp = (sra2_bioproject != my_bioproject)
            is_my_sample_group = (sra2_sample_group == my_sample_group)
            run_other_bp = sra2_run[(is_not_my_bp | !is_my_sample_group)]
            run_other_bp = run_other_bp[run_other_bp %in% tc_runs]
            num_other_run_same_bp_sample_group = length(unique(sra2_bioproject[is_not_my_bp & is_my_sample_group]))
            num_other_bp_same_sample_group = sum(is_not_my_bp & is_my_sample_group, na.rm = TRUE)
            run_context_runs = sra2_run[(!is_not_my_bp) & is_my_sample_group]
            run_context_runs = run_context_runs[run_context_runs %in% tc_runs]
            run_context_cache[[run_context_key]] = list(
                run_other_bp = run_other_bp,
                num_other_run_same_bp_sample_group = num_other_run_same_bp_sample_group,
                num_other_bp_same_sample_group = num_other_bp_same_sample_group,
                run_context_runs = run_context_runs
            )
        }
        sra2[is_sra, "num_other_run_same_bp_sample_group"] = num_other_run_same_bp_sample_group

        # If one sample_group is completely sourced from the same bioproject, we can't remove the whole bioproject for tc_ave_other_bp

        if (num_other_bp_same_sample_group == 0) {
            tc_ave_other_bp = tc_ave_all
        } else {
            cache_key = ifelse(length(run_other_bp) == 0, '__EMPTY__', paste(run_other_bp, collapse = '|'))
            if (!is.null(tc_ave_other_bp_cache[[cache_key]])) {
                tc_ave_other_bp = tc_ave_other_bp_cache[[cache_key]]
            } else {
                tc_other_bp = tc[, run_other_bp]
                tc_ave_other_bp = sample_group_mean(tc_other_bp, sra2, selected_sample_groups)[['tc_ave']]
                tc_ave_other_bp_cache[[cache_key]] = tc_ave_other_bp
            }
        }
        coef = coef_matrix_all[sra_run, selected_sample_groups]
        if (num_other_bp_same_sample_group == 0) {
            coef_other_bp = coef
        } else {
            if (is.null(coef_other_bp_cache[[run_context_key]])) {
                coef_other_bp_mat = cor(tc[, run_context_runs, drop = FALSE], tc_ave_other_bp, method = dist_method)
                if (is.null(dim(coef_other_bp_mat))) {
                    coef_other_bp_mat = matrix(coef_other_bp_mat, nrow = length(run_context_runs), ncol = ncol(tc_ave_other_bp))
                }
                if (is.null(rownames(coef_other_bp_mat))) {
                    rownames(coef_other_bp_mat) = run_context_runs
                }
                if (is.null(colnames(coef_other_bp_mat))) {
                    colnames(coef_other_bp_mat) = colnames(tc_ave_other_bp)
                }
                coef_other_bp_cache[[run_context_key]] = coef_other_bp_mat
            }
            coef_other_bp = coef_other_bp_cache[[run_context_key]][sra_run, selected_sample_groups, drop = TRUE]
        }
        names(coef) = selected_sample_groups
        names(coef_other_bp) = selected_sample_groups
        coef[my_sample_group] = coef[my_sample_group] - min_dif
        coef_other_bp[my_sample_group] = coef_other_bp[my_sample_group] - min_dif
        if (max(coef, na.rm = TRUE) != coef[my_sample_group]) {
            cat('Registered as a candidate for exclusion. Better correlation to other sample group(s):', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
        if (coef_other_bp[my_sample_group] < correlation_threshold) {
            cat('Registered as a candidate for exclusion. Low within-sample-group correlation:', sra_run, '\n')
            exclude_runs = c(exclude_runs, sra_run)
        }
    }
    if (length(exclude_runs)) {
        if (one_out_per_iter == TRUE) {
            cat("Excluding only one outlier per BioProject or same sample_group. \n")
            exclude_run_bps_and_sample_group = sra2[(sra2$run %in% exclude_runs), c("bioproject", "run", "sample_group")]
            first_bp_hit = exclude_run_bps_and_sample_group[match(unique(exclude_run_bps_and_sample_group$bioproject), exclude_run_bps_and_sample_group$bioproject),]
            first_same_sample_group_hit = exclude_run_bps_and_sample_group[match(unique(exclude_run_bps_and_sample_group$sample_group), exclude_run_bps_and_sample_group$sample_group),]
            # if a first_same_sample_group_hit is part of the same bioproject as the other removal candidates, ommit the same sample_group candidates
            if (any(first_same_sample_group_hit$bioproject %in% first_bp_hit$bioproject)) {
                exclude_runs_tmp = c(first_bp_hit$run, first_same_sample_group_hit[!first_same_sample_group_hit$bioproject %in% first_bp_hit$bioproject]$run)
            }else {
                exclude_runs_tmp = c(first_bp_hit$run, first_same_sample_group_hit$run)
            }
            exclude_runs = unique(as.character(exclude_runs_tmp))
            exclude_runs = exclude_runs[(!is.na(exclude_runs)) & (exclude_runs != '')]
            exclude_bps = unique(sra2[(sra2[['run']] %in% exclude_runs), "bioproject"])
        } else {
            exclude_run_bps = sra2[(sra2[['run']] %in% exclude_runs), c("bioproject", "run", "num_other_run_same_bp_sample_group")]
            exclude_bp_counts = data.frame(table(exclude_run_bps[['bioproject']]))
            exclude_run_bps = merge(exclude_run_bps, exclude_bp_counts, by.x = "bioproject", by.y = "Var1")
            exclude_run_bps = exclude_run_bps[order(exclude_run_bps[['num_other_run_same_bp_sample_group']], exclude_run_bps[['Freq']]),]
            rownames(exclude_run_bps) = 1:nrow(exclude_run_bps)
            min_other_run_same_bp_sample_group = exclude_run_bps[1, "num_other_run_same_bp_sample_group"]
            semimin_bp_count = exclude_run_bps[1, "Freq"]
            cat("Minimum number of other BioProjects within sample_group:", min_other_run_same_bp_sample_group, "\n")
            cat("Semi-minimum count of exclusion-candidate BioProjects:", semimin_bp_count, "\n")
            conditions = (exclude_run_bps[['Freq']] == semimin_bp_count)
            conditions = conditions & (exclude_run_bps[['num_other_run_same_bp_sample_group']] == min_other_run_same_bp_sample_group)
            exclude_bps = unique(exclude_run_bps[conditions, "bioproject"])
            exclude_runs = exclude_run_bps[(exclude_run_bps[['bioproject']] %in% exclude_bps), "run"]
        }
    }
    if (length(exclude_runs)) {
        cat('Partially removed BioProjects due to low within-sample_group correlation:', paste(exclude_bps, collapse = ' '), '\n')
        cat('Removed Runs due to low within-sample_group correlation:', paste(exclude_runs, collapse = ' '), '\n')
    }
    tc = tc[, !colnames(tc) %in% exclude_runs, drop = FALSE]
    sra[(sra[['run']] %in% exclude_runs), "exclusion"] = "low_within_sample_group_correlation"
    return(list(tc = tc, sra = sra))
}

batch_effect_subtraction = function(tc, sra, batch_effect_alg, transform_method, clip_negative) {
    if (ncol(tc) == 1) {
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
        mod = try(model.matrix(~sample_group, data = sra))
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
        bp_freq = as.data.frame(table(sra[, "bioproject"]))
        bp_freq_gt1 = bp_freq[bp_freq[, "Freq"] > 1, "Var1"]
        bp_freq_eq1 = bp_freq[bp_freq[, "Freq"] == 1, "Var1"]
        run_bp_freq_gt1 = sra[sra[, "bioproject"] %in% bp_freq_gt1, "run"]
        run_bp_freq_eq1 = sra[sra[, "bioproject"] %in% bp_freq_eq1, "run"]
        tc_combat = tc[, colnames(tc) %in% run_bp_freq_gt1]
        tcc_cn = colnames(tc_combat)
        batch = sra[(sra[, "bioproject"] %in% bp_freq_gt1), "bioproject"]
        group = sra[(sra[, "run"] %in% run_bp_freq_gt1), "sample_group"]
        tc_combat = try(ComBat_seq(as.matrix(tc_combat), batch = batch, group = group))
        if (class(tc_combat)[1] != "try-error") {
            cat("These runs are being removed, due to the bioproject only having 1 sample: \n")
            print(run_bp_freq_eq1)
            cat("Combatseq correction was correctly performed.\n")
            tc_combat = as.data.frame(tc_combat)
            colnames(tc_combat) = tcc_cn
            tc = rbind(tc_combat, tc_ne[, colnames(tc_combat)])
            sva1 = ''
        } else {
            cat("Combatseq correction failed. Trying again without group parameter. \n")
            tc_combat = tc[, colnames(tc) %in% run_bp_freq_gt1]
            tc_combat = try(ComBat_seq(as.matrix(tc_combat), batch = batch))
            if (class(tc_combat)[1] != "try-error") {
                cat("These runs are being removed, due to the bioproject only having 1 sample: \n")
                print(run_bp_freq_eq1)
                cat("Combatseq correction was correctly performed.\n")
                tc_combat = as.data.frame(tc_combat)
                colnames(tc_combat) = tcc_cn
                tc = rbind(tc_combat, tc_ne[, colnames(tc_combat)])
                sva1 = ''
            }
            else {
                stop("Combatseq correction failed.")
            }
        }
    } else if (batch_effect_alg == "ruvseq") {
        x = as.factor(sra$sample_group)
        design = try(model.matrix(~sample_group, data = sra))
        if ("try-error" %in% class(design)) {
            return(list(tc = tc, sva = NULL))
        }
        y = DGEList(counts = as.matrix(tc + 1), group = x)
        y = calcNormFactors(y, method = "upperquartile")
        y = estimateGLMCommonDisp(y, design)
        y = estimateGLMTagwiseDisp(y, design)
        fit = glmFit(y, design)
        res = residuals(fit, type = "deviance")
        seqUQ = betweenLaneNormalization(as.matrix(tc + 1), which = "upper", round = TRUE, offset = FALSE)
        controls = rep(TRUE, dim(as.matrix(tc + 1))[1])
        batch_ruv_res = try(RUVr(seqUQ, controls, k = 1, res)[[2]])
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
            val1 = format(sum(is_negative), big.mark = ',', scientific = FALSE)
            val2 = format(prod(dim(tc)), big.mark = ',', scientific = FALSE)
            val3 = format(dim(tc)[1], big.mark = ',', scientific = FALSE)
            val4 = format(dim(tc)[2], big.mark = ',', scientific = FALSE)
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

compute_dense_label_pt = function(num_labels, base_pt = 8, min_pt = 4, soft_limit = 20) {
    if (is.na(num_labels) || (num_labels <= 0)) {
        return(base_pt)
    }
    if (num_labels <= soft_limit) {
        return(base_pt)
    }
    scaled_pt = base_pt * sqrt(soft_limit / num_labels)
    max(min_pt, min(base_pt, scaled_pt))
}

shorten_run_labels = function(run_ids, max_len = 10, prefix_len = 4, suffix_len = 4) {
    run_ids = as.character(run_ids)
    short_ids = run_ids
    is_long = nchar(run_ids) > max_len
    if (any(is_long)) {
        head_half = substr(run_ids[is_long], 1, prefix_len)
        start_idx = pmax(1, nchar(run_ids[is_long]) - suffix_len + 1)
        tail_half = mapply(
            function(id, st) { substr(id, st, nchar(id)) },
            run_ids[is_long],
            start_idx,
            USE.NAMES = FALSE
        )
        short_ids[is_long] = paste0(head_half, '..', tail_half)
    }
    # Avoid duplicated labels, which can break tile indexing.
    make.unique(short_ids, sep = '.')
}

prepare_sample_correlation_heatmap_data = function(sra, tc_dist_matrix) {
    run_order = colnames(tc_dist_matrix)
    run_order = run_order[run_order %in% sra[['run']]]
    if (length(run_order) == 0) {
        return(NULL)
    }
    sra = sra[match(run_order, sra[['run']]), , drop = FALSE]
    bp_primary = sub(';.*', '', as.character(sra[['bioproject']]))
    sample_group_chr = as.character(sra[['sample_group']])
    run_chr = as.character(sra[['run']])
    order_idx = order(sample_group_chr, bp_primary, run_chr, na.last = TRUE)
    run_order = run_order[order_idx]
    sra = sra[order_idx, , drop = FALSE]
    tc_dist_matrix = as.matrix(tc_dist_matrix[run_order, run_order, drop = FALSE])
    colnames(tc_dist_matrix) = sra[sra$run %in% colnames(tc_dist_matrix), 'run']
    rownames(tc_dist_matrix) = colnames(tc_dist_matrix)
    short_names = shorten_run_labels(colnames(tc_dist_matrix), max_len = 10, prefix_len = 4, suffix_len = 4)
    rownames(tc_dist_matrix) = short_names
    colnames(tc_dist_matrix) = short_names
    sra[['bioproject_primary']] = bp_primary[order_idx]
    return(list(sra = sra, tc_dist_matrix = tc_dist_matrix))
}

draw_heatmap_base = function(sra, tc_dist_matrix, legend = TRUE, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    out = prepare_sample_correlation_heatmap_data(sra = sra, tc_dist_matrix = tc_dist_matrix)
    if (is.null(out)) {
        plot(c(0, 1), c(0, 1), type = 'n', ann = FALSE, bty = 'n', xaxt = 'n', yaxt = 'n')
        return(invisible(NULL))
    }
    sra = out[['sra']]
    tc_dist_matrix = out[['tc_dist_matrix']]
    n = ncol(tc_dist_matrix)
    if (n == 0) {
        plot(c(0, 1), c(0, 1), type = 'n', ann = FALSE, bty = 'n', xaxt = 'n', yaxt = 'n')
        return(invisible(NULL))
    }

    ann_gap = 0.25
    top_sg_y = n + 1 + ann_gap
    top_bp_y = n + 2 + ann_gap
    left_sg_x = 0 - ann_gap
    left_bp_x = -1 - ann_gap

    rdylbu11 = c(
        '#A50026', '#D73027', '#F46D43', '#FDAE61', '#FEE090', '#FFFFBF',
        '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#313695'
    )
    heat_breaks = c(0, seq(0.3, 1, 0.01))
    heat_colors = grDevices::colorRampPalette(rev(rdylbu11))(length(heat_breaks) - 1)
    value_to_color = function(x) {
        idx = findInterval(x, heat_breaks, all.inside = TRUE)
        heat_colors[idx]
    }

    bp_values = as.character(sra[['bioproject_primary']])
    sample_group_values = as.character(sra[['sample_group']])
    bp_color_values = as.character(sra[['bp_color']])
    sample_group_color_values = as.character(sra[['sample_group_color']])
    row_labels = rownames(tc_dist_matrix)
    col_labels = colnames(tc_dist_matrix)

    has_legend = isTRUE(legend)
    legend_pad = ifelse(has_legend, 2.1, 0)
    sample_label_pt = compute_dense_label_pt(num_labels = n, base_pt = fontsize, min_pt = 4, soft_limit = 20)
    cex_txt = sample_label_pt / fontsize
    cex_ann = 1
    longest_row_label = if (length(row_labels) > 0) max(nchar(row_labels), na.rm = TRUE) else 0
    row_label_pad = max(2.0, 0.22 * longest_row_label)
    left_axis_x = left_bp_x - row_label_pad
    row_label_x = left_axis_x + 0.12
    plot.new()
    plot.window(
        xlim = c(left_axis_x, n + 0.5 + legend_pad),
        ylim = c(0.5, top_bp_y + 1.6),
        xaxs = 'i',
        yaxs = 'i'
    )

    # Correlation tiles
    for (ri in seq_len(n)) {
        y_center = n - ri + 1
        y_bottom = y_center - 0.49
        y_top = y_center + 0.49
        for (ci in seq_len(n)) {
            x_center = ci
            x_left = x_center - 0.49
            x_right = x_center + 0.49
            rect(
                x_left, y_bottom, x_right, y_top,
                col = value_to_color(tc_dist_matrix[ri, ci]),
                border = NA
            )
        }
    }
    # Thin grid lines for readability
    for (k in seq(0.5, n + 0.5, by = 1)) {
        segments(0.5, k, n + 0.5, k, col = 'white', lwd = 0.4)
        segments(k, 0.5, k, n + 0.5, col = 'white', lwd = 0.4)
    }

    # Annotation strips (same size as heatmap cells)
    for (i in seq_len(n)) {
        yi = n - i + 1
        rect(i - 0.49, top_sg_y - 0.49, i + 0.49, top_sg_y + 0.49, col = sample_group_color_values[i], border = NA)
        rect(i - 0.49, top_bp_y - 0.49, i + 0.49, top_bp_y + 0.49, col = bp_color_values[i], border = NA)
        rect(left_sg_x - 0.49, yi - 0.49, left_sg_x + 0.49, yi + 0.49, col = sample_group_color_values[i], border = NA)
        rect(left_bp_x - 0.49, yi - 0.49, left_bp_x + 0.49, yi + 0.49, col = bp_color_values[i], border = NA)
    }

    # Labels
    text(x = seq_len(n), y = top_bp_y + 1.10, labels = col_labels, srt = 90, adj = c(1, 0.5), cex = cex_txt, col = 'black', xpd = TRUE)
    text(x = c(left_bp_x, left_sg_x), y = top_bp_y + 1.10, labels = c('BioProject', 'Sample group'),
         srt = 90, adj = c(1, 0.5), cex = cex_ann, col = 'black', xpd = TRUE)
    text(x = row_label_x, y = n:1, labels = row_labels, adj = c(0, 0.5), cex = cex_txt, col = 'black', xpd = TRUE)

    # Border around heatmap body
    rect(0.5, 0.5, n + 0.5, n + 0.5, border = 'black', lwd = 0.8)

    # Optional correlation colorbar
    if (has_legend) {
        bar_xleft = n + 1.0
        bar_xright = n + 1.25
        bar_ybottom = 0.5
        bar_ytop = n + 0.5
        nbar = 200
        ycuts = seq(bar_ybottom, bar_ytop, length.out = nbar + 1)
        vcuts = seq(0, 1, length.out = nbar)
        for (i in seq_len(nbar)) {
            rect(bar_xleft, ycuts[i], bar_xright, ycuts[i + 1], col = value_to_color(vcuts[i]), border = NA)
        }
        rect(bar_xleft, bar_ybottom, bar_xright, bar_ytop, border = 'black', lwd = 0.8)
        text(
            bar_xright + 0.14,
            bar_ytop + 1.25,
            labels = "Pearson's\ncorrelation\ncoefficient",
            adj = c(0, 1),
            cex = cex_ann,
            col = 'black',
            xpd = TRUE
        )
        tick_vals = c(0, 0.5, 1.0)
        tick_ys = bar_ybottom + (bar_ytop - bar_ybottom) * tick_vals
        for (i in seq_along(tick_vals)) {
            segments(bar_xright, tick_ys[i], bar_xright + 0.08, tick_ys[i], xpd = TRUE, col = 'black')
            text(bar_xright + 0.10, tick_ys[i], labels = sprintf('%.1f', tick_vals[i]), adj = c(0, 0.5), cex = cex_ann, col = 'black', xpd = TRUE)
        }
    }
    invisible(NULL)
}

build_numeric_heatmap_ggplot = function(mat, breaks, colors, fontsize = 8, legend_title = '') {
    fontsize = resolve_curate_fontsize(fontsize)
    mat = as.matrix(mat)
    if ((nrow(mat) == 0) || (ncol(mat) == 0)) {
        return(ggplot2::ggplot() + ggplot2::theme_void(base_size = fontsize, base_family = CURATE_FONT_FAMILY))
    }
    df_heat = as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df_heat) = c('row_name', 'col_name', 'value')
    row_levels = rev(rownames(mat))
    col_levels = colnames(mat)
    df_heat[['row_name']] = factor(df_heat[['row_name']], levels = row_levels)
    df_heat[['col_name']] = factor(df_heat[['col_name']], levels = col_levels)
    g = ggplot2::ggplot(df_heat, ggplot2::aes(x = col_name, y = row_name, fill = value)) +
        ggplot2::geom_tile(width = 0.98, height = 0.98) +
        ggplot2::scale_fill_gradientn(
            colors = colors,
            values = scales::rescale(breaks, from = range(breaks)),
            limits = range(breaks),
            oob = scales::squish,
            na.value = 'gray80',
            name = legend_title
        ) +
        ggplot2::theme_bw(base_size = fontsize, base_family = CURATE_FONT_FAMILY) +
        ggplot2::theme(
            text = ggplot2::element_text(size = fontsize, family = CURATE_FONT_FAMILY, color = 'black'),
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = ggplot2::element_blank()
        )
    return(g)
}

draw_heatmap = function(sra, tc_dist_matrix, legend = TRUE, show_colorbar = legend, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    g = build_heatmap_ggplot(
        sra = sra,
        tc_dist_matrix = tc_dist_matrix,
        legend = legend,
        show_colorbar = show_colorbar,
        fontsize = fontsize
    )
    draw_ggplot_in_current_plot_panel(g)
}

# ggplot2 heatmap for sample-correlation with row/column annotations.
# This keeps the same information content: correlation tiles + bioproject/sample_group strips on top and left.
build_heatmap_ggplot = function(sra, tc_dist_matrix, legend = TRUE, show_colorbar = legend, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    run_order = colnames(tc_dist_matrix)
    run_order = run_order[run_order %in% sra[['run']]]
    if (length(run_order) == 0) {
        g = ggplot2::ggplot() + ggplot2::theme_void(base_size = fontsize, base_family = CURATE_FONT_FAMILY)
        return(g)
    }
    sra = sra[match(run_order, sra[['run']]), , drop = FALSE]
    bp_primary = sub(';.*', '', as.character(sra[['bioproject']]))
    sample_group_chr = as.character(sra[['sample_group']])
    run_chr = as.character(sra[['run']])
    order_idx = order(sample_group_chr, bp_primary, run_chr, na.last = TRUE)
    run_order = run_order[order_idx]
    sra = sra[order_idx, , drop = FALSE]
    tc_dist_matrix = as.matrix(tc_dist_matrix[run_order, run_order, drop = FALSE])
    colnames(tc_dist_matrix) = sra[sra$run %in% colnames(tc_dist_matrix), 'run']
    rownames(tc_dist_matrix) = colnames(tc_dist_matrix)
    short_names = shorten_run_labels(colnames(tc_dist_matrix), max_len = 10, prefix_len = 4, suffix_len = 4)
    rownames(tc_dist_matrix) = short_names
    colnames(tc_dist_matrix) = short_names
    n = ncol(tc_dist_matrix)
    if (n == 0) {
        g = ggplot2::ggplot() + ggplot2::theme_void(base_size = fontsize, base_family = CURATE_FONT_FAMILY)
        return(g)
    }
    sample_label_pt = compute_dense_label_pt(num_labels = n, base_pt = fontsize, min_pt = 4, soft_limit = 20)

    # Heatmap matrix values
    row_names = rownames(tc_dist_matrix)
    col_names = colnames(tc_dist_matrix)
    df_heat = as.data.frame(as.table(tc_dist_matrix), stringsAsFactors = FALSE)
    colnames(df_heat) = c('row_label', 'col_label', 'value')
    df_heat[['row_idx']] = match(df_heat[['row_label']], row_names)
    df_heat[['col_idx']] = match(df_heat[['col_label']], col_names)
    df_heat[['x']] = df_heat[['col_idx']]
    df_heat[['y']] = n - df_heat[['row_idx']] + 1

    # Row/column annotations (same information as the previous heatmap implementation)
    bp_values = sub(';.*', '', as.character(sra[['bioproject']]))
    sample_group_values = as.character(sra[['sample_group']])
    annotation_df = data.frame(
        run = run_order,
        x = seq_len(n),
        y = n - seq_len(n) + 1,
        sample_group = sample_group_values,
        bioproject = bp_values,
        sample_group_color = as.character(sra[['sample_group_color']]),
        bioproject_color = as.character(sra[['bp_color']]),
        stringsAsFactors = FALSE
    )

    ann_gap = 0.25
    top_sg_y = n + 1 + ann_gap
    top_bp_y = n + 2 + ann_gap
    left_sg_x = 0 - ann_gap
    left_bp_x = -1 - ann_gap

    ann_points = rbind(
        data.frame(
            x = annotation_df[['x']],
            y = rep(top_sg_y, n),
            legend_group = 'Sample group',
            legend_label = annotation_df[['sample_group']],
            color = annotation_df[['sample_group_color']],
            stringsAsFactors = FALSE
        ),
        data.frame(
            x = annotation_df[['x']],
            y = rep(top_bp_y, n),
            legend_group = 'BioProject',
            legend_label = annotation_df[['bioproject']],
            color = annotation_df[['bioproject_color']],
            stringsAsFactors = FALSE
        ),
        data.frame(
            x = rep(left_sg_x, n),
            y = annotation_df[['y']],
            legend_group = 'Sample group',
            legend_label = annotation_df[['sample_group']],
            color = annotation_df[['sample_group_color']],
            stringsAsFactors = FALSE
        ),
        data.frame(
            x = rep(left_bp_x, n),
            y = annotation_df[['y']],
            legend_group = 'BioProject',
            legend_label = annotation_df[['bioproject']],
            color = annotation_df[['bioproject_color']],
            stringsAsFactors = FALSE
        )
    )
    bp_palette_df = unique(ann_points[ann_points[['legend_group']] == 'BioProject', c('legend_label', 'color')])
    bp_palette = bp_palette_df[['color']]
    names(bp_palette) = bp_palette_df[['legend_label']]
    sample_group_palette_df = unique(ann_points[ann_points[['legend_group']] == 'Sample group', c('legend_label', 'color')])
    sample_group_palette = sample_group_palette_df[['color']]
    names(sample_group_palette) = sample_group_palette_df[['legend_label']]

    # Match the previous heatmap settings:
    # color = "-RdYlBu2:71", breaks = c(0, seq(0.3, 1, 0.01)).
    # Use a smooth continuous colorbar while preserving the same value mapping.
    heat_breaks = c(0, seq(0.3, 1, 0.01))
    rdylbu11 = c(
        '#A50026', '#D73027', '#F46D43', '#FDAE61', '#FEE090', '#FFFFBF',
        '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#313695'
    )
    heat_colors = grDevices::colorRampPalette(rev(rdylbu11))(length(heat_breaks))

    g = ggplot2::ggplot() +
        ggplot2::geom_tile(
            data = df_heat,
            mapping = ggplot2::aes(x = x, y = y, fill = value),
            width = 0.98,
            height = 0.98
        ) +
        ggplot2::geom_tile(
            data = ann_points,
            mapping = ggplot2::aes(x = x, y = y),
            fill = ann_points[['color']],
            width = 0.98,
            height = 0.98,
            show.legend = FALSE
        ) +
        # Keep annotation legend entries while drawing annotation strips as tiles.
        ggplot2::geom_point(
            data = ann_points[ann_points[['legend_group']] == 'BioProject', , drop = FALSE],
            mapping = ggplot2::aes(x = x, y = y, color = legend_label),
            shape = 15,
            size = 0.01,
            alpha = 0,
            show.legend = legend
        ) +
        ggplot2::geom_point(
            data = ann_points[ann_points[['legend_group']] == 'Sample group', , drop = FALSE],
            mapping = ggplot2::aes(x = x, y = y, shape = legend_label),
            size = 0.01,
            alpha = 0,
            show.legend = legend
        ) +
        ggplot2::scale_fill_gradientn(
            colors = heat_colors,
            values = scales::rescale(heat_breaks, from = c(0, 1)),
            limits = c(0, 1),
            na.value = 'gray80',
            name = "Pearson's\ncorrelation\ncoefficient",
            breaks = c(0, 0.5, 1.0),
            guide = if (show_colorbar) ggplot2::guide_colorbar(order = 1, nbin = 256) else 'none'
        ) +
        ggplot2::scale_color_manual(
            values = bp_palette,
            name = "BioProject",
            guide = ggplot2::guide_legend(order = 2, override.aes = list(shape = 15, size = fontsize * 0.6, alpha = 1, colour = unname(bp_palette)))
        ) +
        ggplot2::scale_shape_manual(
            values = setNames(rep(15, length(sample_group_palette)), names(sample_group_palette)),
            name = "Sample group",
            guide = ggplot2::guide_legend(
                order = 3,
                override.aes = list(size = fontsize * 0.6, alpha = 1, colour = unname(sample_group_palette))
            )
        ) +
        ggplot2::scale_x_continuous(
            breaks = c(left_bp_x, left_sg_x, seq_len(n)),
            labels = c('BioProject', 'Sample group', col_names),
            expand = c(0, 0),
            position = 'top'
        ) +
        ggplot2::scale_y_continuous(
            breaks = c(seq_len(n), top_sg_y, top_bp_y),
            labels = c(rev(row_names), 'Sample group', 'BioProject'),
            expand = c(0, 0)
        ) +
        ggplot2::coord_fixed(
            xlim = c(left_bp_x - 0.5, n + 0.5),
            ylim = c(0.5, top_bp_y + 0.5),
            clip = 'on'
        ) +
        ggplot2::theme_bw(base_size = fontsize, base_family = CURATE_FONT_FAMILY) +
        ggplot2::theme(
            text = ggplot2::element_text(size = fontsize, family = CURATE_FONT_FAMILY, color = 'black'),
            panel.grid = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = sample_label_pt, color = 'black', angle = 90, hjust = 0, vjust = 0.5),
            axis.text.y = ggplot2::element_text(size = sample_label_pt, color = 'black'),
            legend.title = ggplot2::element_text(size = fontsize, color = 'black'),
            legend.text = ggplot2::element_text(size = fontsize, color = 'black'),
            legend.box = 'vertical',
            legend.key.height = grid::unit(0.14, 'in'),
            legend.spacing.y = grid::unit(0.04, 'in'),
            legend.position = if (legend || show_colorbar) 'right' else 'none',
            plot.margin = ggplot2::margin(6, 6, 6, 6, unit = 'pt')
        )

    if (!show_colorbar) {
        g = g + ggplot2::guides(fill = 'none')
    }
    if (!legend) {
        g = g + ggplot2::guides(color = 'none', shape = 'none')
    }
    return(g)
}

draw_heatmap_ggplot = function(sra, tc_dist_matrix, legend = TRUE, show_colorbar = legend, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    g = build_heatmap_ggplot(
        sra = sra,
        tc_dist_matrix = tc_dist_matrix,
        legend = legend,
        show_colorbar = show_colorbar,
        fontsize = fontsize
    )
    print(g)
    invisible(g)
}

save_heatmap_ggplot = function(sra, tc_dist_matrix, out_path, legend = TRUE, fontsize = 8, width = 3.6, height = 3.6,
                               legend_extra_width = 2.0, show_colorbar = legend) {
    fontsize = resolve_curate_fontsize(fontsize)
    g = build_heatmap_ggplot(
        sra = sra,
        tc_dist_matrix = tc_dist_matrix,
        legend = legend,
        show_colorbar = show_colorbar,
        fontsize = fontsize
    )
    out_width = if (legend || show_colorbar) width + legend_extra_width else width
    ggplot2::ggsave(filename = out_path, plot = g, width = out_width, height = height, units = 'in')
    invisible(g)
}

draw_dendrogram = function(sra, tc_dist_dist, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    dend = as.dendrogram(hclust(tc_dist_dist))
    dend_colors = sra[order.dendrogram(dend), 'sample_group_color']
    dend = apply_leaf_label_colors(dend, dend_colors)
    dend = apply_leaf_edge_colors(dend, dend_colors)
    for (i in seq_len(nrow(sra))) {
        dend = dendrapply(dend, color_children2parent)
    }
    dend = set_edge_lwd(dend, lwd = 1)
    sample_label_pt = compute_dense_label_pt(num_labels = nrow(sra), base_pt = fontsize, min_pt = 3.5, soft_limit = 18)
    sample_label_cex = sample_label_pt / fontsize
    plot(dend, las = 1, axes = FALSE, cex = sample_label_cex)
    axis(side = 2, line = 0, las = 1, cex.axis = 1)
    mtext("Distance", side = 2, line = 8.5, outer = FALSE, cex = 1)
    n = nrow(sra)
    symbols(
        1:n,
        rep(0, n),
        circles = rep(1, n),
        add = TRUE,
        inches = 0.02,
        xpd = TRUE,
        lwd = 1,
        bg = sra[order.dendrogram(dend), 'sample_group_color'],
        fg = sra[order.dendrogram(dend), 'bp_color']
    )
}

draw_dendrogram_ggplot = function(sra, tc_dist_dist, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)

    cg_col = unique(sra[, c('sample_group', 'sample_group_color')])
    bp_col = unique(sra[, c('bioproject', 'bp_color')])
    colnames(cg_col) = c('Group', 'Color')
    colnames(bp_col) = c('Group', 'Color')
    sra_colors = rbind(cg_col, bp_col)

    group_colors <- data.table(Group = sra_colors$Group, Color = sra_colors$Color, key = "Group")
    group_colors <- transpose(group_colors, make.names = "Group")
    colnames(tc_dist_dist) <- sra$run
    hc <- hclust(tc_dist_dist)           # heirarchal clustering
    dendr <- dendro_data(hc, type = "rectangle") # convert for ggplot


    clust.df <- data.frame(label = sra$run, sample_group = factor(sra[sra$run %in% dendr$labels$label, 'sample_group']), bioproject = factor(sra[sra$run %in% dendr$labels$label, 'bioproject']))
    dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
    ggplot() +
        geom_segment(data = dendr$segments, aes(x = x, y = y, xend = xend, yend = yend), size = .8, show.legend = FALSE) +
        geom_segment(data = merge(dendr$segments[dendr$segments$yend == 0,], dendr$labels[, c('label', 'sample_group', 'x')], by = 'x'), aes(x = x, y = y, xend = xend, yend = yend, color = sample_group), size = .8, show.legend = FALSE) +
        geom_text(data = dendr$labels, aes(x, y - .008, label = label, hjust = 0, angle = 270, color = sample_group), family = CURATE_FONT_FAMILY, size = fontsize / ggplot2::.pt, show.legend = FALSE) +
        scale_y_continuous(expand = c(.2, .1)) +
        geom_point(data = dendr$labels, aes(x, y, color = bioproject), size = 3, show.legend = FALSE) +
        geom_point(data = dendr$labels, aes(x, y, color = sample_group), size = 2, show.legend = FALSE) +
        theme(text = element_text(size = fontsize, family = CURATE_FONT_FAMILY),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_blank()
        ) +
        scale_color_manual(values = group_colors)
}

draw_dendrogram_pvclust = function(sra, tc, nboot, pvclust_file, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    cat("draw_dendrogram_pvclust(): bootstrap support is deprecated. Using hclust without bootstrap.\n")
    dend = as.dendrogram(hclust(calc_sample_distance(tc, method = 'pearson'), method = "average"))
    dend_colors = sra[order.dendrogram(dend), 'sample_group_color']
    dend = apply_leaf_label_colors(dend, dend_colors)
    dend = apply_leaf_edge_colors(dend, dend_colors)
    for (i in seq_len(ncol(tc))) {
        dend = dendrapply(dend, color_children2parent)
    }
    dend = set_edge_lwd(dend, lwd = 2)
    cex.xlab = min(0.2 + 1 / log10(fontsize), na.rm = TRUE)
    par(cex = cex.xlab)
    plot(dend, las = 1, ylab = "Distance", cex.axis = 1 / cex.xlab, cex.lab = 1 / cex.xlab)
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
        bg = sra[order.dendrogram(dend), 'sample_group_color'],
        fg = sra[order.dendrogram(dend), 'bp_color']
    )
}

draw_pca = function(sra, tc_dist_matrix, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    pca = prcomp(tc_dist_matrix)
    xlabel = paste0("PC 1 (", round(summary(pca)[['importance']][2, 1] * 100, digits = 1), "%)")
    ylabel = paste0("PC 2 (", round(summary(pca)[['importance']][2, 2] * 100, digits = 1), "%)")
    plot(
        pca[['x']][, 1],
        pca[['x']][, 2],
        pch = 21,
        cex = 2,
        lwd = 1,
        bg = sra[['sample_group_color']],
        col = sra[['bp_color']],
        xlab = xlabel,
        ylab = ylabel,
        las = 1
    )
    # plot(pca$x[,1], pca$x[,2], pch=21, cex=2, lwd=2, bg=sra$sample_group_color, col=sra$bp_color, main=title,
    # xlab=xlabel, ylab=ylabel, las=1)
}

draw_mds = function(sra, tc_dist_dist, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    mds_points = tryCatch({
        as.matrix(stats::cmdscale(tc_dist_dist, k = 2))
    }, error = function(a) {
        NULL
    })
    if (is.null(mds_points) || (ncol(mds_points) < 2)) {
        cat("MDS failed.\n")
        plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    } else {
        plot(
            mds_points[, 1],
            mds_points[, 2],
            pch = 21,
            cex = 2,
            lwd = 1,
            bg = sra[['sample_group_color']],
            col = sra[['bp_color']],
            xlab = "MDS dimension 1",
            ylab = "MDS dimension 2",
            las = 1
        )
    }
}

draw_tsne = function(sra, tc, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    num_samples = suppressWarnings(as.integer(min(nrow(sra), ncol(tc), na.rm = TRUE)))
    if (!is.finite(num_samples) || is.na(num_samples) || (num_samples < 4)) {
        cat("t-SNE skipped: at least 4 samples are required.\n")
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        return(NULL)
    }
    max_perplexity = floor((num_samples - 1) / 3)
    if (!is.finite(max_perplexity) || (max_perplexity < 1)) {
        cat("t-SNE skipped: unable to determine a valid perplexity.\n")
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        return(NULL)
    }
    perplexity = min(30, max_perplexity)
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
        plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
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
            bg = sra[['sample_group_color']],
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
    fontsize = resolve_curate_fontsize(fontsize)
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
        cols = c("sample_group", "bioproject", "lib_layout", "lib_selection", "instrument", "mapping_rate", 'log10_total_spots', 'log10_total_bases')
        label_cols = c("Sample group", "BioProject", "Library layout", "Library selection", "Instrument", "Mapping rate", 'Log10 total reads', 'Log10 total bases')
        if ('tmm_normalization_factor' %in% colnames(sra)) {
            cols = c(cols, 'tmm_normalization_factor')
            label_cols = c(label_cols, 'TMM normalization factor')
        }
        num_sv = sva_out[['n.sv']]
        if (num_sv == 0) {
            cat('No surrogate variables found.\n')
            plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
            text(
                x = 0.5,
                y = 0.58,
                labels = "No SV detected",
                cex = 1,
                font = 2
            )
            text(
                x = 0.5,
                y = 0.42,
                labels = "Adjusted R-squared not available",
                cex = 1
            )
            return(data.frame())
        }
        df = data.frame(matrix(NA, num_sv, length(cols)))
        colnames(df) = cols
        rownames(df) = paste0("SV", 1:nrow(df))
        for (i in seq_along(cols)) {
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
        g = build_numeric_heatmap_ggplot(mat = df2, breaks = breaks, colors = colors, fontsize = fontsize, legend_title = 'Adj. R-squared')
        draw_ggplot_in_current_plot_panel(g)
    }
    return(df)
}

draw_boxplot = function(sra, tc_dist_matrix, fontsize = 8) {
    fontsize = resolve_curate_fontsize(fontsize)
    is_same_bp = outer(sra[['bioproject']], sra[['bioproject']], function(x, y) {
        x == y
    })
    is_same_sample_group = outer(sra[['sample_group']], sra[['sample_group']], function(x, y) {
        x == y
    })
    group_bwbw = tc_dist_matrix[(!is_same_bp) & (!is_same_sample_group)]
    group_wibw = tc_dist_matrix[(is_same_bp) & (!is_same_sample_group)]
    group_bwwi = tc_dist_matrix[(!is_same_bp) & (is_same_sample_group)]
    group_wiwi = tc_dist_matrix[(is_same_bp) & (is_same_sample_group)]
    sanitize_group_values = function(x) {
        x = as.numeric(x)
        x[is.finite(x)]
    }
    group_bwbw = sanitize_group_values(group_bwbw)
    group_wibw = sanitize_group_values(group_wibw)
    group_bwwi = sanitize_group_values(group_bwwi)
    group_wiwi = sanitize_group_values(group_wiwi)
    draw_group_boxplot = function(values, at) {
        if (length(values) > 0) {
            boxplot(values, at = at, add = TRUE, col = "gray", yaxt = "n")
        }
    }
    plot(c(0.5, 4.5), c(0, 1), type = "n", xlab = "", ylab = "Pearson's correlation\ncoefficient", las = 1,
         xaxt = "n")
    draw_group_boxplot(group_bwbw, 1)
    draw_group_boxplot(group_wibw, 2)
    draw_group_boxplot(group_bwwi, 3)
    draw_group_boxplot(group_wiwi, 4)
    labels = c("bw\nbw", "bw\nwi", "wi\nbw", "wi\nwi")
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = "Sample group\nBioProject", padj = 0.5, hadj = 1, tick = FALSE)

    # Add mean PCC
    means <- c(
        if (length(group_bwbw) > 0) mean(group_bwbw, na.rm = TRUE) else NA_real_,
        if (length(group_wibw) > 0) mean(group_wibw, na.rm = TRUE) else NA_real_,
        if (length(group_bwwi) > 0) mean(group_bwwi, na.rm = TRUE) else NA_real_,
        if (length(group_wiwi) > 0) mean(group_wiwi, na.rm = TRUE) else NA_real_
    )
    finite_idx = which(is.finite(means))
    if (length(finite_idx) > 0) {
        points(finite_idx, means[finite_idx], col = "red", pch = 16)
    }
    if (all(is.finite(means[1:2]))) {
        lines(c(1, 2), means[1:2], col = "red")
        text(x = 1.5, y = max(means[1:2]) + 0.05, labels = round(abs(means[1] - means[2]), 2), col = "red", cex = 1)
    }
    if (all(is.finite(means[3:4]))) {
        lines(c(3, 4), means[3:4], col = "red")
        text(x = 3.5, y = max(means[3:4]) + 0.05, labels = round(abs(means[3] - means[4]), 2), col = "red", cex = 1)
    }

    labels = c("bw\nbw", "bw\nwi", "wi\nbw", "wi\nwi")
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = "Sample group\nBioProject", padj = 0.5, hadj = 1, tick = FALSE)

    # Add legend in the bottom left corner
    legend("bottomleft", legend = c("mean PCC", expression(Delta ~ "mean PCC")),
           pch = c(16, NA), col = c("red", "red"),
           lty = c(NA, 1), bty = "n", cex = 1,
           text.width = max(strwidth(c("mean PCC", expression(Delta ~ "mean PCC")))))  # Align text

}

save_correlation = function(tc, sra, dist_method, round, precomputed_tc_dist_matrix = NULL) {
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    if (!is.null(precomputed_tc_dist_matrix)) {
        run_order = colnames(precomputed_tc_dist_matrix)
        run_order = run_order[run_order %in% sra[['run']]]
        tc_dist_matrix = precomputed_tc_dist_matrix[run_order, run_order, drop = FALSE]
        sra = sra[match(run_order, sra[['run']]), , drop = FALSE]
    } else {
        tc_dist_matrix = NULL
    }
    if ((is.null(tc_dist_matrix) && (ncol(tc) <= 1)) || (!is.null(tc_dist_matrix) && (ncol(tc_dist_matrix) <= 1))) {
        cat('Not enough samples for correlation statistics. Recording NA values.\n')
        tc_dist_stats = rep(NA_real_, 12)
        if (!exists("correlation_statistics", envir = .GlobalEnv)) {
            correlation_statistics <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean", "bwbw_median", "bwbw_variance", "wibw_mean", "wibw_median", "wibw_variance", "bwwi_mean", "bwwi_median", "bwwi_variance", "wiwi_mean", "wiwi_median", "wiwi_variance"))))
            rownames(correlation_statistics) <- paste0("round_", round)
            assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
        } else {
            correlation_statistics <- get("correlation_statistics", envir = .GlobalEnv)
            new_row <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean", "bwbw_median", "bwbw_variance", "wibw_mean", "wibw_median", "wibw_variance", "bwwi_mean", "bwwi_median", "bwwi_variance", "wiwi_mean", "wiwi_median", "wiwi_variance"))))
            rownames(new_row) <- paste0("round_", round)
            correlation_statistics <- rbind(correlation_statistics, new_row)
            assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
        }
        return(correlation_statistics)
    }
    is_same_bp = outer(sra[['bioproject']], sra[['bioproject']], function(x, y) {
        x == y
    })
    is_same_sample_group = outer(sra[['sample_group']], sra[['sample_group']], function(x, y) {
        x == y
    })


    if (is.null(tc_dist_matrix)) {
        tc_dist_matrix = cor(tc, method = dist_method)
    }
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0

    # bw = between; wi = within; order: Bioproject -> Sample Group; tc_dist_bwwi: between BPs, within CGs
    tc_dist_bwbw = tc_dist_matrix[(!is_same_bp) & (!is_same_sample_group)]
    bwbw_mea = mean(tc_dist_bwbw)
    bwbw_med = median(tc_dist_bwbw)
    bwbw_var = var(tc_dist_bwbw)

    tc_dist_wibw = tc_dist_matrix[(is_same_bp) & (!is_same_sample_group)]
    wibw_mea = mean(tc_dist_wibw)
    wibw_med = median(tc_dist_wibw)
    wibw_var = var(tc_dist_wibw)

    tc_dist_bwwi = tc_dist_matrix[(!is_same_bp) & (is_same_sample_group)]
    bwwi_mea = mean(tc_dist_bwwi)
    bwwi_med = median(tc_dist_bwwi)
    bwwi_var = var(tc_dist_bwwi)

    tc_dist_wiwi = tc_dist_matrix[(is_same_bp) & (is_same_sample_group)]
    wiwi_mea = mean(tc_dist_wiwi)
    wiwi_med = median(tc_dist_wiwi)
    wiwi_var = var(tc_dist_wiwi)

    tc_dist_stats = c(bwbw_mea, bwbw_med, bwbw_var, wibw_mea, wibw_med, wibw_var, bwwi_mea, bwwi_med, bwwi_var, wiwi_mea, wiwi_med, wiwi_var)

    # Check if dataframe exists in environment, if not create it
    if (!exists("correlation_statistics", envir = .GlobalEnv)) {
        correlation_statistics <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean", "bwbw_median", "bwbw_variance", "wibw_mean", "wibw_median", "wibw_variance", "bwwi_mean", "bwwi_median", "bwwi_variance", "wiwi_mean", "wiwi_median", "wiwi_variance"))))
        rownames(correlation_statistics) <- paste0("round_", round)
        assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
    } else {
        correlation_statistics <- get("correlation_statistics", envir = .GlobalEnv)
        new_row <- data.frame(matrix(tc_dist_stats, ncol = 12, dimnames = list(NULL, c("bwbw_mean", "bwbw_median", "bwbw_variance", "wibw_mean", "wibw_median", "wibw_variance", "bwwi_mean", "bwwi_median", "bwwi_variance", "wiwi_mean", "wiwi_median", "wiwi_variance"))))
        rownames(new_row) <- paste0("round_", round)
        correlation_statistics <- rbind(correlation_statistics, new_row)
        assign("correlation_statistics", correlation_statistics, envir = .GlobalEnv)
    }

    return(correlation_statistics)


}


draw_tau_histogram = function(tc, sra, selected_sample_groups, fontsize = 8, transform_method, tc_sample_group = NULL) {
    fontsize = resolve_curate_fontsize(fontsize)
    if (is.null(tc_sample_group)) {
        tc_sample_group = sample_group_mean(tc, sra, selected_sample_groups)[['tc_ave']]
    }
    df_tau = sample_group2tau(tc_sample_group, rich.annotation = FALSE, transform_method)
    hist_out = hist(df_tau[['tau']], breaks = seq(0, 1, 0.05), las = 1, xlab = "Tau (expression specificity)",
                    ylab = "Gene count", main = "", col = "gray")
    num_noexp = sum(is.na(df_tau[['tau']]))
    num_all = nrow(df_tau)
    text_noexp = paste0("Excluded due to\nno expression:\n", num_noexp, "/", num_all, " genes")
    text(0, max(hist_out[['counts']], na.rm = TRUE) * 0.85, text_noexp, pos = 4, cex = 1)
}

draw_exp_level_histogram = function(tc, sra, selected_sample_groups, fontsize = 8, transform_method, tc_sample_group = NULL) {
    fontsize = resolve_curate_fontsize(fontsize)
    if (is.null(tc_sample_group)) {
        tc_sample_group = sample_group_mean(tc, sra, selected_sample_groups)[['tc_ave']]
    }
    xmax = apply(tc_sample_group, 1, max)
    xmax[xmax < 0] = 0
    xmax[xmax > 15] = 15
    breaks = seq(0, 15, 1)
    hist_out = hist(xmax, breaks = breaks, las = 1, xlab = paste0("Max expression (", transform_method, ")"), ylab = "Gene count",
                    main = "", col = "gray")
}

draw_legend = function(sra, new = TRUE, pos = "center", fontsize = 8, nlabel.in.col) {
    fontsize = resolve_curate_fontsize(fontsize)
    if (new) {
        plot.new()
    }
    sample_group_unique = unique(sra[['sample_group']])
    bp_unique = unique(sub(";.*", "", sra[['bioproject']]))
    sample_group_color_unique = unique(sra[['sample_group_color']])
    bp_color_unique = unique(sra[['bp_color']])
    ncol = ceiling((length(sample_group_unique) +
        length(bp_unique) +
        2) / nlabel.in.col)
    legend_text = c("Sample group", as.character(sample_group_unique), "", "BioProject", as.character(bp_unique))
    legend_color = c(rgb(1, 1, 1, 0), rep(rgb(1, 1, 1, 0), length(sample_group_color_unique)), rgb(1, 1, 1,
                                                                                                   0), rgb(1, 1, 1, 0), bp_color_unique)
    legend_bg = c(rgb(1, 1, 1, 0), sample_group_color_unique, rgb(1, 1, 1, 0), rgb(1, 1, 1, 0), rep(rgb(1,
                                                                                                        1, 1, 0), length(bp_color_unique)))
    legend_font = c(2, rep(1, length(sample_group_color_unique)), 1, 2, rep(1, length(bp_color_unique)))
    legend(pos, legend = legend_text, pch = 21, lwd = 1, lty = 0, col = legend_color, pt.bg = legend_bg,
           text.font = legend_font, cex = 1, pt.cex = 1, ncol = ncol, bty = "n")
}

save_plot = function(tc, sra, sva_out, dist_method, file, selected_sample_groups, sample_group_colors, fontsize = 8, transform_method, batch_effect_alg) {
    fontsize = resolve_curate_fontsize(fontsize)
    if (ncol(tc) <= 1) {
        cat('1 or fewer samples are available. Skipping the plot.\n')
        return(invisible(NULL))
    }
    out = tc_metadata_intersect(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    sra = add_color_to_curate_metadata(sra, selected_sample_groups, sample_group_colors)
    out = sort_tc_and_metadata(tc, sra)
    tc = out[["tc"]]
    sra = out[["sra"]]
    colnames(tc) <- sra$run
    tc_dist_matrix = cor(tc, method = dist_method)
    if (dist_method %in% c('pearson', 'spearman', 'kendall')) {
        tc_dist_dist = calc_sample_distance(
            tc,
            method = dist_method,
            na_fill = 1,
            epsilon = 1e-09,
            cor_mat = tc_dist_matrix
        )
    } else {
        tc_dist_dist = calc_sample_distance(tc, method = dist_method, na_fill = 1, epsilon = 1e-09)
    }
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    with_pdf_defaults(
        file = file.path(dir_pdf, paste0(file, ".pdf")),
        width = 7.2,
        height = 8,
        font_size = fontsize,
        font_family = CURATE_FONT_FAMILY,
        plot_fn = function(local_font_size, local_font_family) {
            layout_matrix = matrix(c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1,
                                     2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1,
                                     1, 1, 1, 1, 1, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 3, 3,
                                     3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9,
                                     9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10,
                                     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), 14, 12,
                                   byrow = TRUE)
            layout(layout_matrix)
            par(mar = c(6, 6, 1, 0))
            draw_dendrogram(sra, tc_dist_dist, local_font_size)
            par(mar = c(0, 0, 0, 0))
            draw_heatmap(sra, tc_dist_matrix, legend = FALSE, show_colorbar = TRUE, fontsize = local_font_size)
            par(mar = c(4, 4, 0.1, 1))
            draw_pca(sra, tc_dist_matrix, local_font_size)
            par(mar = c(4, 4, 0.1, 1))
            draw_mds(sra, tc_dist_dist, local_font_size)
            par(mar = c(4, 4, 0.1, 1))
            draw_tsne(sra, tc, local_font_size)
            par(mar = c(4, 5, 0.1, 1))
            draw_boxplot(sra, tc_dist_matrix, local_font_size)
            tc_sample_group = sample_group_mean(tc, sra, selected_sample_groups)[['tc_ave']]
            par(mar = c(4, 4, 1, 1))
            draw_exp_level_histogram(tc, sra, selected_sample_groups, local_font_size, transform_method, tc_sample_group = tc_sample_group)
            par(mar = c(4, 4, 1, 1))
            draw_tau_histogram(tc, sra, selected_sample_groups, local_font_size, transform_method, tc_sample_group = tc_sample_group)
            par(mar = rep(0.1, 4))
            if (batch_effect_alg == 'sva') {
                df_r2 = draw_sva_summary(sva_out, tc, sra, local_font_size)
                if (!all(is.na(df_r2))) {
                    write.table(df_r2, file.path(dir_tsv, paste0(file, ".r2.tsv")), sep = "\t", row.names = FALSE)
                }
            }
            par(mar = rep(0.1, 4))
            draw_legend(sra, new = TRUE, pos = "center", fontsize = local_font_size, nlabel.in.col = 8)
        }
    )
    return(invisible(tc_dist_matrix))
}

transform_raw_to_fpkm = function(counts, effective_lengths, sra) {
    if ('tmm_library_size' %in% colnames(sra)) {
        cat('FPKM transformation with the original library sizes from amalgkit cstmm output.\n')
        cat('If --input_dir is specified with amalgkit cstmm output, resultant values will be TMM-FPKM.\n')
        cat('If --input_dir is specified with amalgkit merge output, resultant values will be non-TMM-corrected FPKM.\n')
        tmp = sra
        rownames(tmp) = tmp[['run']]
        library_sizes = tmp[colnames(counts), 'tmm_library_size']
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
                                      step = c('before_batch', 'before_batch_plot', 'after_batch'), sra) {
    if (batch_effect_alg == 'no') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        }
    } else if (batch_effect_alg == 'sva') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        }
    } else if (batch_effect_alg == 'ruvseq') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        }
    } else if (batch_effect_alg == 'combatseq') {
        if (step == 'before_batch') {
            bool_fpkm_tpm = FALSE
            bool_log = FALSE
        } else if (step == 'before_batch_plot') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        } else if (step == 'after_batch') {
            bool_fpkm_tpm = TRUE
            bool_log = TRUE
        }
    }
    if (bool_fpkm_tpm == TRUE) {
        if (grepl('fpkm', transform_method)) {
            cat('Applying FPKM transformation.\n')
            tc = transform_raw_to_fpkm(tc, tc_eff_length[, colnames(tc)], sra)
        } else if (grepl('tpm', transform_method)) {
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
    for (col in c('instrument', 'bioproject')) {
        if (!col %in% colnames(sra_all)) {
            next
        }
        is_missing = (sra_all[, col] == "") | (is.na(sra_all[, col]))
        sra_all[is_missing, col] = "not_provided"
    }
    return(sra_all)
}

get_species_metadata = function(sra_all, scientific_name, selected_sample_groups) {
    is_sp = (sra_all[, 'scientific_name'] == scientific_name)
    is_sample_group = (sra_all[, 'sample_group'] %in% selected_sample_groups)
    cat('Number of SRA runs for this species:', sum(is_sp), '\n')
    cat('Number of SRA runs for selected sample groups:', sum(is_sample_group), '\n')
    sra = sra_all[(is_sp & is_sample_group),]
    conditions = is_non_excluded_flag(sra[['exclusion']]) & (!sra[['run']] %in% colnames(tc))
    if (any(conditions)) {
        cat("Failed quantification:", sra[conditions, "run"], "\n")
        sra[conditions, "exclusion"] = "failed_quantification"
    }
    return(sra)
}

exclude_inappropriate_sample_from_tc = function(tc, sra) {
    is_not_excluded = is_non_excluded_flag(sra[['exclusion']])
    cat('Number of non-excluded SRA runs (exclusion=="no"):', sum(is_not_excluded), '\n')
    tc = tc[, sra[is_not_excluded, 'run'], drop = FALSE]
    return(tc)
}

exclude_inappropriate_sample_from_eff_length = function(tc_eff_length, tc) {
    tc_eff_length = tc_eff_length[, colnames(tc), drop = FALSE]
    return(tc_eff_length)
}

initialize_round_summary = function() {
    data.frame(
        step = character(0),
        round = integer(0),
        reason = character(0),
        num_runs_before = integer(0),
        num_runs_after = integer(0),
        num_runs_removed = integer(0),
        removed_runs = character(0),
        stringsAsFactors = FALSE
    )
}

append_round_summary = function(round_summary, step, round, reason, runs_before, runs_after) {
    removed_runs = setdiff(runs_before, runs_after)
    removed_runs_txt = if (length(removed_runs) == 0) '' else paste(removed_runs, collapse = ' ')
    new_row = data.frame(
        step = step,
        round = as.integer(round),
        reason = reason,
        num_runs_before = as.integer(length(runs_before)),
        num_runs_after = as.integer(length(runs_after)),
        num_runs_removed = as.integer(length(removed_runs)),
        removed_runs = removed_runs_txt,
        stringsAsFactors = FALSE
    )
    rbind(round_summary, new_row)
}

write_curation_summaries = function(round_summary, sra, scientific_name, batch_effect_alg,
                                    dir_tsv, mapping_rate_cutoff, correlation_threshold,
                                    one_outlier_per_iteration, num_total_runs_species,
                                    num_runs_after_sample_group_filter, script_start_elapsed) {
    species_tag = gsub(" ", "_", scientific_name)
    file_round = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".curation_round_summary.tsv"))
    write.table(round_summary, file = file_round, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    is_kept = is_non_excluded_flag(sra[['exclusion']])
    kept_runs = sra[is_kept, 'run']
    excluded_runs = sra[!is_kept, 'run']
    kept_runs_txt = if (length(kept_runs) == 0) '' else paste(kept_runs, collapse = ' ')
    excluded_runs_txt = if (length(excluded_runs) == 0) '' else paste(excluded_runs, collapse = ' ')
    total_runtime_sec = as.numeric(proc.time()[["elapsed"]] - script_start_elapsed)
    final_summary = data.frame(
        scientific_name = scientific_name,
        batch_effect_alg = batch_effect_alg,
        mapping_rate_cutoff = mapping_rate_cutoff,
        correlation_threshold = correlation_threshold,
        one_outlier_per_iteration = as.logical(one_outlier_per_iteration),
        total_runtime_sec = round(total_runtime_sec, 6),
        num_total_runs_in_species = as.integer(num_total_runs_species),
        num_runs_after_sample_group_filter = as.integer(num_runs_after_sample_group_filter),
        num_runs_final_kept = as.integer(sum(is_kept)),
        num_runs_final_excluded = as.integer(sum(!is_kept)),
        final_kept_runs = kept_runs_txt,
        final_excluded_runs = excluded_runs_txt,
        stringsAsFactors = FALSE
    )
    file_final = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".curation_final_summary.tsv"))
    write.table(final_summary, file = file_final, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

########cd############END OF FUNCTION DECLARATION####################################################

################################################# START OF SVA CORRECTION ########################################################


fontsize = CURATE_FONT_SIZE_PT

tc = read.table(est_counts_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, row.names = 1, check.names = FALSE)
tc_eff_length = read.table(eff_length_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", fill = FALSE, row.names = 1, check.names = FALSE)

sra_all = read.table(metadata_path, sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
sra_all = standardize_metadata_all(sra_all)

scientific_name = sra_all[(sra_all[['run']] %in% colnames(tc)), "scientific_name"][1]
species_tag = gsub(" ", "_", scientific_name)
script_start_elapsed = as.numeric(proc.time()[["elapsed"]])
num_total_runs_species = sum(sra_all[, 'scientific_name'] == scientific_name, na.rm = TRUE)
dir_curate = file.path(out_dir, 'curate')
dir_pdf = file.path(dir_curate, species_tag, 'plots')
dir.create(dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir_rdata = file.path(dir_curate, species_tag, 'rdata')
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir_tsv = file.path(dir_curate, species_tag, 'tables')
dir.create(dir_tsv, showWarnings = FALSE, recursive = TRUE)
setwd(dir_curate)
cat(log_prefix, "Working at:", getwd(), "\n")

sra = get_species_metadata(sra_all, scientific_name, selected_sample_groups)
num_runs_after_sample_group_filter = nrow(sra)
round_summary = initialize_round_summary()
tc = exclude_inappropriate_sample_from_tc(tc, sra)
out = sort_tc_and_metadata(tc, sra); tc = out[["tc"]]; sra = out[["sra"]]
tc_eff_length = exclude_inappropriate_sample_from_eff_length(tc_eff_length, tc)
tc = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch', sra = sra)
tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
is_input_zero = data.frame(tc_tmp == 0, check.names = FALSE)
file_name = file.path(dir_tsv, paste0(species_tag, ".uncorrected.tc.tsv"))
write_table_with_index_name(df = tc_tmp, file_path = file_name, index_name = 'target_id')
original_sample_groups = selected_sample_groups
out = sample_group_mean(tc_tmp, sra, selected_sample_groups)
tc_sample_group_uncorrected = out[['tc_ave']]
selected_sample_groups = out[['selected_sample_groups']]
if (length(selected_sample_groups) != length(original_sample_groups)) {
    if (!(length(sample_group_colors) == 1 && sample_group_colors == "DEFAULT")) {
        sample_group_colors = sample_group_colors[match(selected_sample_groups, original_sample_groups)]
    }
}
file_name = file.path(dir_tsv, paste0(species_tag, ".uncorrected.sample_group.mean.tsv"))
write_table_with_index_name(df = tc_sample_group_uncorrected, file_path = file_name, index_name = 'target_id')

if (skip_curation_flag == TRUE) {
    cat("No curation requested, finishing early.\n")
    file_metadata = file.path(dir_tsv, paste0(species_tag, ".metadata.tsv"))
    write.table(sra[, colnames(sra) != 'index'], file = file_metadata, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    file_corrected_tc = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".tc.tsv"))
    write_table_with_index_name(df = tc_tmp, file_path = file_corrected_tc, index_name = 'target_id')
    file_corrected_mean = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".sample_group.mean.tsv"))
    write_table_with_index_name(df = tc_sample_group_uncorrected, file_path = file_corrected_mean, index_name = 'target_id')
    cat("Files created: \n")
    cat(file.path(dir_tsv, paste0(species_tag, ".uncorrected.tc.tsv")), "\n")
    cat(file.path(dir_tsv, paste0(species_tag, ".uncorrected.sample_group.mean.tsv")), "\n")
    cat(file_corrected_tc, "\n")
    cat(file_corrected_mean, "\n")
    cat(file_metadata, "\n")
    cat("Transformation applied: ", transform_method, "\n")
    round_summary = append_round_summary(
        round_summary = round_summary,
        step = 'skip_curation',
        round = -1,
        reason = 'skip_curation_requested',
        runs_before = colnames(tc),
        runs_after = colnames(tc)
    )
    write_curation_summaries(
        round_summary = round_summary,
        sra = sra,
        scientific_name = scientific_name,
        batch_effect_alg = batch_effect_alg,
        dir_tsv = dir_tsv,
        mapping_rate_cutoff = mapping_rate_cutoff,
        correlation_threshold = correlation_threshold,
        one_outlier_per_iteration = one_outlier_per_iteration,
        num_total_runs_species = num_total_runs_species,
        num_runs_after_sample_group_filter = num_runs_after_sample_group_filter,
        script_start_elapsed = script_start_elapsed
    )
    cat(log_prefix, "Completed.\n")
    quit(save = 'no', status = 0)
}

cat("Removing samples with mapping rate of 0.\n")
round = 0
sva_out = NULL
tc_batch_corrected = NULL
runs_before = colnames(tc)
out = check_mapping_rate(tc, sra, 0)
tc = out[["tc"]]
sra = out[["sra"]]
runs_after = colnames(tc)
round_summary = append_round_summary(
    round_summary = round_summary,
    step = 'mapping_rate_zero',
    round = round,
    reason = 'mapping_rate_below_or_equal_0',
    runs_before = runs_before,
    runs_after = runs_after
)
tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
tc_dist_matrix = save_plot(tc_tmp, sra, NULL, dist_method, paste0(species_tag, ".", round, ".original"),
                           selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
save_correlation(tc_tmp, sra, dist_method, round, precomputed_tc_dist_matrix = tc_dist_matrix)
out = batch_effect_subtraction(tc, sra, batch_effect_alg, transform_method, clip_negative)
tc_batch_corrected = out[["tc"]]
if (batch_effect_alg == "combatseq") {
    tc = tc[, colnames(tc_batch_corrected)]
}
sva_out = out[["sva"]]
if (!is.null(sva_out)) {
    file_name = paste0(species_tag, ".", batch_effect_alg, ".", round, ".RData")
    save(sva_out, file = file.path(dir_rdata, file_name))
}
tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step = 'after_batch', sra = sra)
save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, paste0(species_tag, ".", round, ".original", ".", batch_effect_alg),
          selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)

cat("Removing samples with the mapping rate smaller than", mapping_rate_cutoff, "\n")
round = 1
sva_out = NULL
tc_batch_corrected = NULL
runs_before = colnames(tc)
out = check_mapping_rate(tc, sra, mapping_rate_cutoff)
tc = out[["tc"]]
sra = out[["sra"]]
runs_after = colnames(tc)
round_summary = append_round_summary(
    round_summary = round_summary,
    step = 'mapping_rate_cutoff',
    round = round,
    reason = paste0('mapping_rate_below_or_equal_', mapping_rate_cutoff * 100, '_percent'),
    runs_before = runs_before,
    runs_after = runs_after
)

tc_tmp = apply_transformation_logic(tc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
save_plot(tc_tmp, sra, NULL, dist_method, paste0(species_tag, ".", round, ".mapping_cutoff"),
          selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
out = batch_effect_subtraction(tc, sra, batch_effect_alg, transform_method, clip_negative)
tc_batch_corrected = out[["tc"]]
if (batch_effect_alg == "combatseq") {
    tc = tc[, colnames(tc_batch_corrected)]
}
sva_out = out[["sva"]]
if (!is.null(sva_out)) {
    save(sva_out, file = file.path(dir_rdata, paste0(species_tag, ".", batch_effect_alg, ".", round, ".RData")))
}
tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step = 'after_batch', sra = sra)
tc_dist_matrix = save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, paste0(species_tag, ".", round, ".mapping_cutoff", ".", batch_effect_alg),
                           selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
save_correlation(tc_batch_corrected_tmp, sra, dist_method, round, precomputed_tc_dist_matrix = tc_dist_matrix)

round = 2
end_flag = 0
while (end_flag == 0) {
    cat("Iteratively checking within-sample_group correlation, round:", round, "\n")
    tc_cwtc = NULL
    runs_before = colnames(tc)
    num_run_before = sum(is_non_excluded_flag(sra[['exclusion']]))
    out = check_within_sample_group_correlation(tc, sra, dist_method, min_dif, selected_sample_groups, one_outlier_per_iteration, correlation_threshold)
    tc_cwtc = out[["tc"]]
    sra = out[["sra"]]
    runs_after = colnames(tc_cwtc)
    num_run_after = sum(is_non_excluded_flag(sra[['exclusion']]))
    round_summary = append_round_summary(
        round_summary = round_summary,
        step = 'correlation_iter',
        round = round,
        reason = 'low_within_sample_group_correlation',
        runs_before = runs_before,
        runs_after = runs_after
    )
    if ((num_run_before == num_run_after) | (plot_intermediate)) {
        sva_out = NULL
        tc_batch_corrected = NULL
        out = batch_effect_subtraction(tc[, colnames(tc_cwtc)], sra, batch_effect_alg, transform_method, clip_negative)
        tc_batch_corrected = out[["tc"]]
        if (batch_effect_alg == "combatseq") {
            tc = tc[, colnames(tc_batch_corrected)]
        }
        sva_out = out[["sva"]]
        if (!is.null(sva_out)) {
            save(sva_out, file = file.path(dir_rdata, paste0(species_tag, ".", batch_effect_alg, ".", round, ".RData")))
        }
        tc_cwtc_tmp = apply_transformation_logic(tc_cwtc, tc_eff_length, transform_method, batch_effect_alg, step = 'before_batch_plot', sra = sra)
        save_plot(tc_cwtc_tmp, sra, NULL, dist_method, paste0(species_tag, ".", round, ".correlation_cutoff"),
                  selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
        tc_batch_corrected_tmp = apply_transformation_logic(tc_batch_corrected, tc_eff_length, transform_method, batch_effect_alg, step = 'after_batch', sra = sra)
        file_base = paste0(species_tag, ".", round, ".correlation_cutoff", ".", batch_effect_alg)
        tc_dist_matrix = save_plot(tc_batch_corrected_tmp, sra, sva_out, dist_method, file_base, selected_sample_groups, sample_group_colors, fontsize, transform_method, batch_effect_alg)
        save_correlation(tc_batch_corrected_tmp, sra, dist_method, round, precomputed_tc_dist_matrix = tc_dist_matrix)
    }
    cat("Round:", round, ": # before =", num_run_before, ": # after =", num_run_after, "\n\n")
    if (num_run_before == num_run_after) {
        end_flag = 1
    }
    tc = tc_cwtc
    round = round + 1
}

cat("Finished checking within-sample_group correlation.\n")
if (batch_effect_alg != 'sva') {
    cat("Batch-effect removal algorithm is: ", batch_effect_alg, ". Applying transformation on final batch-removed counts.\n")
    tc_batch_corrected = tc_batch_corrected_tmp
}
if (maintain_zero) {
    cat('Any zero expression levels in the input will remain as zero-values in the output tables.\n')
    tc_batch_corrected = tc_batch_corrected[order(rownames(tc_batch_corrected)),]
    is_input_zero = is_input_zero[rownames(tc_batch_corrected),]
    stopifnot(all(rownames(is_input_zero) == rownames(tc_batch_corrected)))
    for (col in colnames(tc_batch_corrected)) {
        tc_batch_corrected[is_input_zero[[col]], col] = 0
    }
}
cat("Writing summary files for", scientific_name, "\n")
file = file.path(dir_tsv, paste0(species_tag, ".metadata.tsv"))
write.table(sra[, colnames(sra) != 'index'], file = file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
file = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".tc.tsv"))
write_table_with_index_name(df = tc_batch_corrected, file_path = file, index_name = 'target_id')
out = sample_group_mean(tc_batch_corrected, sra, selected_sample_groups)
tc_sample_group = out[['tc_ave']]
file = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".sample_group.mean.tsv"))
write_table_with_index_name(df = tc_sample_group, file_path = file, index_name = 'target_id')
file = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".correlation_statistics.tsv"))
write.table(correlation_statistics, file = file, row.names = TRUE, sep = "\t", quote = FALSE)
tc_tau = sample_group2tau(tc_sample_group, rich.annotation = TRUE, transform_method)
file = file.path(dir_tsv, paste0(species_tag, ".", batch_effect_alg, ".tau.tsv"))
write_table_with_index_name(df = tc_tau, file_path = file, index_name = 'target_id')
write_curation_summaries(
    round_summary = round_summary,
    sra = sra,
    scientific_name = scientific_name,
    batch_effect_alg = batch_effect_alg,
    dir_tsv = dir_tsv,
    mapping_rate_cutoff = mapping_rate_cutoff,
    correlation_threshold = correlation_threshold,
    one_outlier_per_iteration = one_outlier_per_iteration,
    num_total_runs_species = num_total_runs_species,
    num_runs_after_sample_group_filter = num_runs_after_sample_group_filter,
    script_start_elapsed = script_start_elapsed
)
cat(log_prefix, "Completed.\n")
quit(save = 'no', status = 0)
