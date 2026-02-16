get_singlecopy_bool_index = function(df_gc, spp_filled, percent_singlecopy_threshold = 50) {

    is_ge_singlecopy_threshold = function(x, num_sp, percent_singlecopy_threshold) {
        num_singlecopy_species = sum(x == 1)
        percent_singlecopy_species = num_singlecopy_species / num_sp * 100
        is_ge_singlecopy = percent_singlecopy_species >= percent_singlecopy_threshold
        return(is_ge_singlecopy)
    }

    num_sp = length(spp_filled)
    is_singlecopy = apply(df_gc[, spp_filled, drop = FALSE], 1, function(x) { is_ge_singlecopy_threshold(x, num_sp, percent_singlecopy_threshold) })
    num_sc = sum(is_singlecopy)
    txt = 'Number of single-copy orthogroups (>=%s percent species) detected for the %s species: %s\n'
    cat(sprintf(txt, percent_singlecopy_threshold, formatC(length(spp_filled), big.mark = ','), formatC(num_sc, big.mark = ',')))
    return(is_singlecopy)
}

impute_expression = function(dat, num_pc = 4, max_iter = 50, tol = 1e-06, strategy = 'em_pca') {
    impute_with_row_mean = function(mat) {
        mat = as.matrix(mat)
        row_means = rowMeans(mat, na.rm = TRUE)
        row_means[is.na(row_means)] = 0
        for (i in seq_len(nrow(mat))) {
            if (any(is.na(mat[i, ]))) {
                mat[i, is.na(mat[i, ])] = row_means[i]
            }
        }
        return(mat)
    }
    impute_with_iterative_pca = function(mat, n_pc, max_iter, tol) {
        mat = as.matrix(mat)
        missing_mask = is.na(mat)
        if (!any(missing_mask)) {
            return(mat)
        }
        mat_imputed = impute_with_row_mean(mat)
        for (iter in seq_len(max_iter)) {
            col_means = colMeans(mat_imputed)
            centered = sweep(mat_imputed, 2, col_means, '-')
            pca_fit = tryCatch(
                stats::prcomp(centered, center = FALSE, scale. = FALSE, rank. = n_pc),
                error = function(e) { NULL }
            )
            if (is.null(pca_fit)) {
                return(NULL)
            }
            k = min(n_pc, ncol(pca_fit[['x']]), ncol(pca_fit[['rotation']]))
            if (k < 1) {
                return(NULL)
            }
            scores = pca_fit[['x']][, seq_len(k), drop = FALSE]
            loadings = pca_fit[['rotation']][, seq_len(k), drop = FALSE]
            reconstructed = scores %*% t(loadings)
            reconstructed = sweep(reconstructed, 2, col_means, '+')
            old_vals = mat_imputed[missing_mask]
            new_vals = reconstructed[missing_mask]
            if (any(!is.finite(new_vals))) {
                return(NULL)
            }
            mat_imputed[missing_mask] = new_vals
            delta = max(abs(old_vals - new_vals), na.rm = TRUE)
            if (!is.finite(delta) || (delta < tol)) {
                break
            }
        }
        return(mat_imputed)
    }
    fit_nipals_pca = function(mat, n_pc, max_iter, tol) {
        x = as.matrix(mat)
        n_obs = nrow(x)
        n_var = ncol(x)
        max_rank = min(n_obs, n_var)
        n_comp = min(max(1, n_pc), max_rank)
        if (n_comp < 1) {
            return(NULL)
        }
        scores = matrix(0, nrow = n_obs, ncol = n_comp)
        loadings = matrix(0, nrow = n_var, ncol = n_comp)
        colnames(scores) = paste0('PC', seq_len(n_comp))
        colnames(loadings) = paste0('PC', seq_len(n_comp))
        x_resid = x
        for (comp_idx in seq_len(n_comp)) {
            col_var = apply(x_resid, 2, var)
            col_var[!is.finite(col_var)] = -Inf
            init_col = which.max(col_var)
            if ((length(init_col) == 0) || !is.finite(col_var[init_col])) {
                return(NULL)
            }
            t_vec = x_resid[, init_col]
            if ((!all(is.finite(t_vec))) || (sum(t_vec^2) <= .Machine$double.eps)) {
                candidate_cols = which(colSums(x_resid^2) > .Machine$double.eps)
                if (length(candidate_cols) == 0) {
                    break
                }
                t_vec = x_resid[, candidate_cols[1]]
            }
            if ((!all(is.finite(t_vec))) || (sum(t_vec^2) <= .Machine$double.eps)) {
                break
            }
            p_vec = rep(0, n_var)
            for (iter_idx in seq_len(max_iter)) {
                denom_t = sum(t_vec^2)
                if ((!is.finite(denom_t)) || (denom_t <= .Machine$double.eps)) {
                    break
                }
                p_vec = as.vector(crossprod(x_resid, t_vec) / denom_t)
                p_norm = sqrt(sum(p_vec^2))
                if ((!is.finite(p_norm)) || (p_norm <= .Machine$double.eps)) {
                    break
                }
                p_vec = p_vec / p_norm
                t_new = as.vector(x_resid %*% p_vec)
                if (!all(is.finite(t_new))) {
                    return(NULL)
                }
                delta = max(abs(t_new - t_vec), na.rm = TRUE)
                t_vec = t_new
                if ((!is.finite(delta)) || (delta < tol)) {
                    break
                }
            }
            if ((!all(is.finite(t_vec))) || (!all(is.finite(p_vec)))) {
                return(NULL)
            }
            scores[, comp_idx] = t_vec
            loadings[, comp_idx] = p_vec
            x_resid = x_resid - tcrossprod(t_vec, p_vec)
        }
        valid_comp = which(
            (colSums(abs(scores)) > .Machine$double.eps) &
            (colSums(abs(loadings)) > .Machine$double.eps)
        )
        if (length(valid_comp) == 0) {
            return(NULL)
        }
        return(list(
            scores = scores[, valid_comp, drop = FALSE],
            loadings = loadings[, valid_comp, drop = FALSE]
        ))
    }
    impute_with_nipals = function(mat, n_pc, max_iter, tol) {
        mat = as.matrix(mat)
        missing_mask = is.na(mat)
        if (!any(missing_mask)) {
            return(mat)
        }
        mat_imputed = impute_with_row_mean(mat)
        for (iter in seq_len(max_iter)) {
            col_means = colMeans(mat_imputed)
            centered = sweep(mat_imputed, 2, col_means, '-')
            nipals_fit = fit_nipals_pca(
                mat = centered,
                n_pc = n_pc,
                max_iter = max(20, max_iter),
                tol = tol
            )
            if (is.null(nipals_fit)) {
                return(NULL)
            }
            reconstructed = nipals_fit[['scores']] %*% t(nipals_fit[['loadings']])
            reconstructed = sweep(reconstructed, 2, col_means, '+')
            old_vals = mat_imputed[missing_mask]
            new_vals = reconstructed[missing_mask]
            if (any(!is.finite(new_vals))) {
                return(NULL)
            }
            mat_imputed[missing_mask] = new_vals
            delta = max(abs(old_vals - new_vals), na.rm = TRUE)
            if (!is.finite(delta) || (delta < tol)) {
                break
            }
        }
        return(mat_imputed)
    }
    normalize_strategy = function(strategy) {
        if (is.null(strategy) || (length(strategy) == 0) || is.na(strategy[[1]]) || (strategy[[1]] == '')) {
            return('em_pca')
        }
        strategy_norm = tolower(as.character(strategy[[1]]))
        if (strategy_norm %in% c('em_pca', 'emppca', 'em_ppca', 'em-ppca', 'iterative_pca', 'iterative-pca', 'ppca')) {
            return('em_pca')
        }
        if (strategy_norm %in% c('nipals', 'nipals_pca', 'nipals-pca')) {
            return('nipals')
        }
        if (strategy_norm %in% c('row_mean', 'row-mean', 'mean')) {
            return('row_mean')
        }
        stop(sprintf("Unknown missing-value strategy: %s", strategy[[1]]))
    }

    dat = as.matrix(dat)
    strategy = normalize_strategy(strategy)
    is_all_na_row = apply(dat, 1, function(x) { all(is.na(x)) })
    tmp = dat[!is_all_na_row, , drop = FALSE]
    txt = 'Number of removed rows with all NA values in the expression matrix: %s\n'
    cat(sprintf(txt, formatC(sum(is_all_na_row), big.mark = ',')))
    num_na = sum(is.na(tmp))
    num_sp = ncol(tmp)
    num_gene = nrow(tmp)
    num_all = num_sp * num_gene
    txt = 'Imputing %s missing values in a total of %s observations (%s genes x %s samples).\n'
    cat(sprintf(txt, formatC(num_na, big.mark = ','), formatC(num_all, big.mark = ','), formatC(num_gene, big.mark = ','), formatC(num_sp, big.mark = ',')))
    if ((nrow(tmp) == 0) || (ncol(tmp) == 0)) {
        cat('No data available for imputation. Returning empty matrix.\n')
        return(tmp)
    }
    if (strategy == 'row_mean') {
        cat('Missing-value strategy: row_mean.\n')
        imputed_dat = impute_with_row_mean(tmp)
    } else {
        cat(sprintf('Missing-value strategy: %s.\n', strategy))
        max_pc = min(nrow(tmp) - 1, ncol(tmp) - 1)
        if (is.na(max_pc) || max_pc < 1) {
            cat('Too few observations for PCA. Falling back to row-mean imputation.\n')
            imputed_dat = impute_with_row_mean(tmp)
        } else {
            if (num_pc > max_pc) {
                num_pc = max_pc
            }
            if (strategy == 'nipals') {
                imputed_dat = impute_with_nipals(tmp, n_pc = num_pc, max_iter = max_iter, tol = tol)
            } else {
                imputed_dat = impute_with_iterative_pca(tmp, n_pc = num_pc, max_iter = max_iter, tol = tol)
            }
            if (is.null(imputed_dat)) {
                cat('Base PCA imputation failed. Falling back to row-mean imputation.\n')
                imputed_dat = impute_with_row_mean(tmp)
            }
        }
    }
    num_negative = sum(imputed_dat < 0)
    txt = 'Number of negative values clipped to zero in the imputed expression matrix: %s\n'
    cat(sprintf(txt, formatC(num_negative, big.mark = ',')))
    imputed_dat[imputed_dat < 0] = 0
    return(imputed_dat)
}

write_table_with_index_name = function(df, file_path, index_name = 'target_id', sort = TRUE) {
    df_index = data.frame(placeholder_name = rownames(df), stringsAsFactors = FALSE)
    colnames(df_index) = index_name
    df = cbind(df_index, df)
    if (sort) {
        df = df[order(df[[index_name]]),]
    }
    write.table(df, file = file_path, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
}

format_genus_species_label = function(x) {
    x = as.character(x)
    x = gsub('_', ' ', x)
    sapply(strsplit(x, '\\s+'), function(tokens) {
        tokens = tokens[tokens != '']
        if (length(tokens) >= 2) {
            paste(tokens[1], paste(tokens[-1], collapse = ' '), sep = '\n')
        } else if (length(tokens) == 1) {
            tokens[1]
        } else {
            ''
        }
    })
}

save_exclusion_plot = function(df, out_path, font_size, y_label = "Count") {
    data_summary = aggregate(cbind(count = exclusion) ~ scientific_name + exclusion, df, length)
    data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN = sum)
    data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
    g = ggplot(data_summary, aes(x = scientific_name, y = count, fill = exclusion))
    g = g + geom_bar(stat = "identity")
    g = g + scale_x_discrete(labels = format_genus_species_label)
    g = g + labs(x = "", y = y_label, fill = "exclusion")
    g = g + theme_bw(base_size = font_size, base_family = 'Helvetica')
    g = g + theme(
        axis.text = element_text(size = font_size, color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = font_size, color = 'black'),
        #panel.grid.major.y=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = font_size, color = 'black'),
        legend.position = 'bottom',
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(rep(0.1, 4), "cm")
    )
    num_spp = length(unique(df[['scientific_name']]))
    plot_width = max(3.6, 0.11 * num_spp)
    ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')
}
