#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(edgeR, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly = TRUE)))

mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, 'debug', 'batch')

if (mode == "debug") {
    dir_work = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_gfe_pipeline/amalgkit_out'
    dir_count = file.path(dir_work, "merge")
    file_orthogroup_table = file.path(dir_work, 'cstmm', 'cstmm_multispecies_busco_table.tsv')
    file_genecount = file.path(dir_work, 'cstmm', 'cstmm_orthogroup_genecount.tsv')
    dir_cstmm = file.path(dir_work, "cstmm")
    mode_tmm = 'multi_species'
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
    setwd(dir_work)
} else if (mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    dir_count = args[1]
    file_orthogroup_table = args[2]
    file_genecount = args[3]
    dir_cstmm = args[4]
    mode_tmm = args[5]
    r_util_path = args[6]
}
source(r_util_path)

get_spp_filled = function(dir_count, df_gc = NA) {
    sciname_dirs = list.dirs(dir_count, full.names = FALSE, recursive = FALSE)
    spp_filled = c()
    for (sciname_dir in sciname_dirs) {
        count_files = list.files(path = file.path(dir_count, sciname_dir), pattern = ".*est_counts\\.tsv")
        if (length(count_files) == 1) {
            spp_filled = c(spp_filled, count_files)
        } else {
            warning(paste0('Multiple or no est_counts files were detected for ', sciname_dir, ': ', paste(count_files, collapse = ', ')))
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
    infile = list.files(path = sciname_path, pattern = ".*est_counts\\.tsv")
    if (length(infile) > 1) {
        stop(paste0("Multiple *count.tsv files found: ", sp, "\n"))
    } else if (length(infile) == 0) {
        warning(paste0("Skipping. No *est_counts.tsv files found: ", sp, "\n"))
        return(NULL)
    }
    infile_path = file.path(sciname_path, infile[1])
    cat('Input file found, reading:', infile[1], '\n')
    dat = read.delim(infile_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
    dat = dat[, (colnames(dat) != 'length'), drop = FALSE]
    colnames(dat) = paste(sp, colnames(dat), sep = '_')
    return(dat)
}

get_uncorrected = function(dir_count, file_genecount = NA) {
    if (is.na(file_genecount)) {
        df_gc = NA
    } else {
        df_gc = read.table(file_genecount, header = TRUE, sep = '\t', check.names = FALSE, quote = '', comment.char = '')
        rownames(df_gc) = df_gc[['orthogroup_id']]
        df_gc[, 'orthogroup_id'] = NULL
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
    df_gc = read.table(file_genecount, header = TRUE, sep = '\t', check.names = FALSE, quote = '', comment.char = '')
    rownames(df_gc) = df_gc[['orthogroup_id']]
    df_gc[, 'orthogroup_id'] = NULL
    df_og = read.table(file_orthogroup_table, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, quote = '', comment.char = '')
    spp_filled = get_spp_filled(dir_count, df_gc)
    is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
    df_singleog = df_og[is_singlecopy, spp_filled, drop = FALSE]
    df_sog = df_singleog
    for (sp in spp_filled) {
        if (!sp %in% names(uncorrected)) {
            next
        }
        df_sog = merge(df_sog, uncorrected[[sp]], by.x = sp, by.y = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE)
    }
    df_sog = df_sog[, -(1:length(spp_filled))]
    rownames(df_sog) = rownames(df_singleog)
    return(df_sog)
}

get_df_nonzero = function(df_sog, imputation = TRUE) {
    is_no_count_col = apply(df_sog, 2, function(x) { sum(x, na.rm = TRUE) == 0 })
    txt = 'Removing %s out of %s samples whose read mapping values are all zero.\n'
    cat(sprintf(txt, formatC(sum(is_no_count_col), big.mark = ','), formatC(ncol(df_sog), big.mark = ',')))
    df_nonzero = df_sog[, !is_no_count_col]
    if (imputation) {
        df_nonzero = impute_expression(df_nonzero)
    } else {
        is_na_containing_row = apply(df_sog, 1, function(x) { any(is.na(x)) })
        txt = 'Removing %s out of %s orthogroups because missing values are observed in at least one species.\n'
        cat(sprintf(txt, formatC(sum(is_na_containing_row), big.mark = ','), formatC(nrow(df_sog), big.mark = ',')))
        df_nonzero = df_sog[!is_na_containing_row,]
    }
    return(df_nonzero)
}

create_eff_length_symlink = function(dir_count, dir_cstmm, sp) {
    path_sp = file.path(dir_count, sp)
    eff_length_files = list.files(path = path_sp, pattern = ".*eff_length\\.tsv")
    if (length(eff_length_files) == 1) {
        path_target = file.path(path_sp, eff_length_files[1])
        path_link = file.path(dir_cstmm, sp, eff_length_files[1])
        cat('Copying file from', path_target, 'to', path_link, '\n')
        file.copy(from = path_target, to = path_link, overwrite = TRUE)
    } else {
        warning(paste0('No eff_length.tsv file found: ', path_sp))
    }
}

append_tmm_stats_to_metadata = function(df_metadata, cnf_out2) {

    my_fun = function(x) {
        split_sample_name = strsplit(x, '_')[[1]]
        run_name = paste(split_sample_name[3:length(split_sample_name)], collapse = '_')
        return(run_name)
    }

    df_nf = cnf_out2[[2]]
    df_nf[['sample']] = rownames(df_nf)
    df_nf[['scientific_name']] = df_nf[['sample']]
    df_nf[['scientific_name']] = sub('_', 'PLACEHOLDER', df_nf[['scientific_name']])
    df_nf[['scientific_name']] = sub('_.*', '', df_nf[['scientific_name']])
    df_nf[['scientific_name']] = sub('PLACEHOLDER', ' ', df_nf[['scientific_name']])
    df_nf[['run']] = sapply(df_nf[['sample']], function(x) { my_fun(x) })
    df_nf = df_nf[, c('scientific_name', 'run', 'lib.size', 'norm.factors')]
    colnames(df_nf) = c('scientific_name', 'run', 'tmm_library_size', 'tmm_normalization_factor')
    out_cols = c(colnames(df_metadata), colnames(df_nf)[3:ncol(df_nf)])
    df_metadata = merge(df_metadata, df_nf, by = c('scientific_name', 'run'), sort = FALSE, all.x = TRUE, all.y = FALSE)
    df_metadata = df_metadata[, out_cols]
    filled_mapping_rate = df_metadata[['mapping_rate']]
    filled_mapping_rate[is.na(filled_mapping_rate)] = -999
    df_metadata[((filled_mapping_rate == 0) & (df_metadata[['exclusion']] == 'no')), 'exclusion'] = 'no_mapping'
    df_metadata[((!df_metadata[['run']] %in% df_nf[['run']]) & (df_metadata[['exclusion']] == 'no')), 'exclusion'] = 'no_cstmm_output'
    df_metadata[(is.na(df_metadata[['tmm_normalization_factor']]) & (df_metadata[['exclusion']] == 'no')), 'exclusion'] = 'cstmm_failed'
    return(df_metadata)
}

plot_norm_factor_histogram = function(df_metadata, font_size = 8) {
    tmp = df_metadata[(!is.na(df_metadata[['tmm_normalization_factor']])),]
    x_limit = max(abs(log2(tmp[['tmm_normalization_factor']])), na.rm = TRUE)
    for (fill_by in c('scientific_name', 'sample_group')) {
        g = ggplot2::ggplot(tmp) +
            geom_histogram(aes(x = log2(tmm_normalization_factor), fill = !!rlang::sym(fill_by)), position = "stack", alpha = 0.7, bins = 40) +
            theme_bw(base_size = font_size, base_family = 'Helvetica') +
            xlim(c(-x_limit, x_limit)) +
            labs(x = 'log2(TMM normalization factor)', y = 'Count') +
            guides(fill = guide_legend(ncol = 1)) +
            theme(
                axis.text = element_text(size = font_size, color = 'black'),
                axis.title = element_text(size = font_size, color = 'black'),
                #panel.grid.major.y=element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                legend.title = element_blank(),
                legend.text = element_text(size = font_size, color = 'black'),
                rect = element_rect(fill = "transparent"),
                plot.margin = unit(rep(0.1, 4), "cm")
            )
        out_path = file.path(dir_cstmm, paste0('cstmm_normalization_factor_histogram.', fill_by, '.pdf'))
        ggsave(out_path, plot = g, width = 4.8, height = 2.4, units = 'in')
    }
}

plot_norm_factor_scatter = function(df_metadata, font_size = 8) {
    tmp = df_metadata[(!is.na(df_metadata[['tmm_normalization_factor']])),]
    x_limit = max(abs(log2(tmp[['tmm_normalization_factor']])), na.rm = TRUE)
    g = ggplot2::ggplot(tmp, aes(x = log10(tmm_library_size), y = log2(tmm_normalization_factor), fill = scientific_name, color = sample_group)) +
        geom_point(shape = 21, alpha = 0.7) +
        scale_fill_hue(l = 65) +
        scale_color_hue(l = 45) +
        theme_bw(base_size = font_size, base_family = 'Helvetica') +
        ylim(c(-x_limit, x_limit)) +
        labs(x = 'log10(Library size)', y = 'log2(TMM normalization factor)') +
        guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) +
        theme(
            axis.text = element_text(size = font_size, color = 'black'),
            axis.title = element_text(size = font_size, color = 'black'),
            #panel.grid.major.y=element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.title = element_text(size = font_size, color = 'black'),
            legend.text = element_text(size = font_size, color = 'black'),
            rect = element_rect(fill = "transparent"),
            plot.margin = unit(rep(0.1, 4), "cm")
        )
    ggsave(file.path(dir_cstmm, 'cstmm_normalization_factor_scatter.pdf'), plot = g, width = 4.8, height = 2.0, units = 'in')
}

get_library_sizes = function(df_nonzero, uncorrected) {
    df_libsize = data.frame(library_size = rep(NA, ncol(df_nonzero)))
    rownames(df_libsize) = colnames(df_nonzero)
    for (sp in names(uncorrected)) {
        for (sample in colnames(uncorrected[[sp]])) {
            if (sample %in% rownames(df_libsize)) {
                df_libsize[sample, 'library_size'] = sum(uncorrected[[sp]][[sample]], na.rm = TRUE)
            }
        }
    }
    library_sizes = df_libsize[['library_size']]
    return(library_sizes)
}

save_mean_expression_boxplot = function(df_nonzero, cnf_out2, uncorrected, corrected, font_size = 8) {
    df_nonzero_tmm = df_nonzero
    for (col in colnames(df_nonzero_tmm)) {
        tmm_normalization_factor = cnf_out2[[2]][col, 'norm.factors']
        df_nonzero_tmm[, col] = df_nonzero_tmm[, col] / tmm_normalization_factor # manually apply TMM normalization factors
    }
    mean_before = apply(df_nonzero, 2, function(x) { mean(x, na.rm = TRUE) })
    mean_after = apply(df_nonzero_tmm, 2, function(x) { mean(x, na.rm = TRUE) })
    var_before = round(var(mean_before), digits = 1)
    var_after = round(var(mean_after), digits = 1)
    txt = 'Across-species among-sample variance of mean single-copy gene raw counts before and after TMM normalization:'
    cat(txt, var_before, 'and', var_after, '\n')

    mean_sra_before = c()
    mean_sra_after = c()
    for (sp in names(corrected)) {
        mean_sra_before = c(mean_sra_before, apply(uncorrected[[sp]], 2, function(x) { mean(x, na.rm = TRUE) }))
        mean_sra_after = c(mean_sra_before, apply(corrected[[sp]], 2, function(x) { mean(x, na.rm = TRUE) }))
    }
    var_sra_before = round(var(mean_sra_before), digits = 1)
    var_sra_after = round(var(mean_sra_after), digits = 1)
    txt = 'Across-species among-sample variance of all-gene mean raw counts before and after TMM normalization:'
    cat(txt, var_sra_before, 'and', var_sra_after, '\n')

    values = c(mean_before, mean_after)
    labels = c(rep('Raw\ncounts', length(mean_before)), rep('TMM-\ncorrected\ncounts', length(mean_after)))
    df = data.frame(labels = labels, values = values)

    sra_values = c(mean_sra_before, mean_sra_after)
    sra_labels = c(rep('Raw\ncounts', length(mean_sra_before)), rep('TMM-\ncorrected\ncounts', length(mean_sra_after)))
    sra_df = data.frame(labels = sra_labels, values = sra_values)

    ps = list()
    ps[[1]] = ggplot(df, aes(x = labels, y = values)) +
        geom_boxplot(outlier.shape = NA) +
        ylab('Mean count of single-copy genes')
    ps[[2]] = ggplot(sra_df, aes(x = labels, y = values)) +
        geom_boxplot(outlier.shape = NA) +
        ylab('Mean count of all genes')

    for (i in 1:length(ps)) {
        ps[[i]] = ps[[i]] +
            ylim(0, 1000) +
            theme_bw(base_size = font_size, base_family = 'Helvetica') +
            theme(
                axis.text = element_text(size = font_size, color = 'black'),
                axis.title = element_text(size = font_size, color = 'black'),
                legend.text = element_text(size = font_size, color = 'black'),
                legend.title = element_text(size = font_size, color = 'black'),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                rect = element_rect(fill = "transparent"),
                plot.margin = unit(rep(0.1, 4), "cm"),
            )
    }

    filename = file.path(dir_cstmm, 'cstmm_mean_expression_boxplot.pdf')
    grDevices::pdf(file = filename, height = 3.6, width = 3.6, family = 'Helvetica', pointsize = font_size)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
    print(ps[[1]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(ps[[2]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grDevices::dev.off()
}

save_corrected_output_files = function(uncorrected) {
    corrected = list()
    for (sp in names(uncorrected)) {
        dat = uncorrected[[sp]]
        df_nf_sp = cnf_out2[[2]][startsWith(rownames(cnf_out2[[2]]), sp),]
        for (i in 1:length(df_nf_sp[, 1])) {
            SRR = as.character(row.names(df_nf_sp[i,]))
            tmm_normalization_factor = as.double(df_nf_sp[i, "norm.factors"])
            dat[, SRR] = dat[, SRR] / tmm_normalization_factor # manually apply TMM normalization factors
        }
        dat_out = cbind(target_id = rownames(dat), dat)
        rownames(dat_out) = NULL
        colnames(dat_out) = sub(paste0(sp, '_'), '', colnames(dat_out))
        dir_cstmm_sp = file.path(dir_cstmm, sp)
        if (!file.exists(dir_cstmm_sp)) {
            dir.create(dir_cstmm_sp)
        }
        create_eff_length_symlink(dir_count, dir_cstmm, sp)
        file_path = file.path(dir_cstmm_sp, paste0(sp, "_cstmm_counts.tsv"))
        write.table(dat_out, file_path, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
        corrected[[sp]] = dat_out
        rownames(corrected[[sp]]) = corrected[[sp]][['target_id']]
        corrected[[sp]][, 'target_id'] = NULL
    }
    return(corrected)
}

print_excluded_run_summary = function(df_metadata) {
    exclusion_reasons = sort(unique(df_metadata[['exclusion']]))
    for (reason in exclusion_reasons) {
        num_run = sum(df_metadata[['exclusion']] == reason)
        if (reason == 'no') {
            txt_final = paste0('Number of retained runs (exclusion = "no"): ', num_run, '\n')
        } else {
            txt = paste0('Number of excluded runs (exclusion = "', reason, '"): ', num_run, '\n')
            cat(txt)
        }
    }
    cat(txt_final)
}


if (mode_tmm == 'single_species') {
    uncorrected = get_uncorrected(dir_count = dir_count, file_genecount = NA)
    stopifnot(length(names(uncorrected)) == 1)
    sp = names(uncorrected)[[1]]
    df_sog = uncorrected[[sp]]
    df_nonzero = get_df_nonzero(df_sog)
} else if (mode_tmm == 'multi_species') {
    uncorrected = get_uncorrected(dir_count = dir_count, file_genecount = file_genecount)
    df_sog = get_df_exp_single_copy_ortholog(file_genecount, file_orthogroup_table, dir_count, uncorrected)
    df_nonzero = get_df_nonzero(df_sog)
}

libraray_sizes = get_library_sizes(df_nonzero, uncorrected)
cnf_in = edgeR::DGEList(counts = df_nonzero, lib.size = libraray_sizes)
cat('Round 1: Performing TMM normalization to determine the appropriate baseline.\n')
cnf_out1 = edgeR::calcNormFactors(cnf_in, method = 'TMM', refColumn = NULL)
x = cnf_out1[[2]][['norm.factors']]
cat('Round 1: Median TMM normalization factor =', median(x), '\n')
median_value = sort(x)[ceiling(length(x) / 2)]
median_index = (1:length(x))[x == median_value]

cat('Round 2: Performing TMM normalization for output.\n')
cnf_out2 = edgeR::calcNormFactors(cnf_in, method = 'TMM', refColumn = median_index)
cat('Round 2: Median TMM normalization factor =', median(cnf_out2[[2]][['norm.factors']]), '\n')

path_metadata = file.path(dir_count, 'metadata.tsv')
df_metadata = read.table(path_metadata, header = TRUE, sep = '\t', check.names = FALSE, quote = '', comment.char = '')
df_metadata = append_tmm_stats_to_metadata(df_metadata, cnf_out2)
df_metadata = df_metadata[, !startsWith(colnames(df_metadata), 'Unnamed')]
out_path = file.path(dir_cstmm, 'metadata.tsv')
write.table(df_metadata, out_path, row.names = FALSE, sep = '\t', quote = FALSE)
print_excluded_run_summary(df_metadata)
plot_norm_factor_histogram(df_metadata = df_metadata)
plot_norm_factor_scatter(df_metadata = df_metadata)
corrected = save_corrected_output_files(uncorrected)
save_mean_expression_boxplot(df_nonzero, cnf_out2, uncorrected, corrected, font_size = 8)

cat(sprintf('Number of SRA samples for exclusion potting: %s\n', formatC(nrow(df_metadata), format = 'd', big.mark = ',')))
out_path = file.path(dir_cstmm, 'cstmm_exclusion.pdf')
save_exclusion_plot(df = df_metadata, out_path = out_path, font_size = 8)

cat('cstmm.r completed!\n')
