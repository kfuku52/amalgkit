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

impute_expression = function(dat, num_pc = 4, max_iter = 50, tol = 1e-06) {
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

    dat = as.matrix(dat)
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
    max_pc = min(nrow(tmp) - 1, ncol(tmp) - 1)
    if (is.na(max_pc) || max_pc < 1) {
        cat('Too few observations for PCA. Falling back to row-mean imputation.\n')
        imputed_dat = impute_with_row_mean(tmp)
    } else {
        if (num_pc > max_pc) {
            num_pc = max_pc
        }
        imputed_dat = impute_with_iterative_pca(tmp, n_pc = num_pc, max_iter = max_iter, tol = tol)
        if (is.null(imputed_dat)) {
            cat('Base PCA imputation failed. Falling back to row-mean imputation.\n')
            imputed_dat = impute_with_row_mean(tmp)
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

save_exclusion_plot = function(df, out_path, font_size) {
    data_summary = aggregate(cbind(count = exclusion) ~ scientific_name + exclusion, df, length)
    data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN = sum)
    data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
    g = ggplot(data_summary, aes(x = scientific_name, y = count, fill = exclusion))
    g = g + geom_bar(stat = "identity")
    g = g + scale_x_discrete(labels = format_genus_species_label)
    g = g + labs(x = "", y = "Count", fill = "exclusion")
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
