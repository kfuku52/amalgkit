#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly = TRUE)))
mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, 'debug', 'batch')
if (mode == "debug") {
    dir_work = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
    setwd(dir_work)
    dir_merge = file.path(dir_work, 'merge')
    file_metadata = file.path(dir_merge, 'metadata.tsv')
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
} else if (mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    dir_merge = args[1]
    file_metadata = args[2]
    r_util_path = args[3]
}
source(r_util_path)
font_size = resolve_plot_font_size(8)
font_family = resolve_plot_font_family('Helvetica')

normalize_exclusion_values = function(exclusion_values) {
    normalized = as.character(exclusion_values)
    normalized[is.na(exclusion_values)] = NA_character_
    normalized = trimws(normalized)
    normalized = tolower(normalized)
    return(normalized)
}

is_non_excluded_flag = function(exclusion_values) {
    normalized = normalize_exclusion_values(exclusion_values)
    return((!is.na(normalized)) & (normalized == 'no'))
}

df = read.table(file_metadata, sep = '\t', header = TRUE, quote = '', comment.char = '', check.names = FALSE)
is_non_excluded = is_non_excluded_flag(df[['exclusion']])
is_excluded = !is_non_excluded
is_mapping_rate_available = (!is.na(df[['mapping_rate']]))
cat(sprintf('Number of non-excluded SRA samples: %s\n', formatC(sum(is_non_excluded), format = 'd', big.mark = ',')))
cat(sprintf('Number of excluded SRA samples: %s\n', formatC(sum(is_excluded), format = 'd', big.mark = ',')))
cat(sprintf('Number of SRA samples with available mapping rates: %s\n', formatC(sum(is_mapping_rate_available), format = 'd', big.mark = ',')))
if ('fastp_duplication_rate' %in% colnames(df)) {
    is_fastp_duplication_rate_available = (!is.na(df[['fastp_duplication_rate']]))
    cat(sprintf('Number of SRA samples with available fastp duplication rates: %s\n',
                formatC(sum(is_fastp_duplication_rate_available), format = 'd', big.mark = ',')))
}
if ('fastp_insert_size_peak' %in% colnames(df)) {
    is_fastp_insert_size_peak_available = (!is.na(df[['fastp_insert_size_peak']]))
    cat(sprintf('Number of SRA samples with available fastp insert size peaks: %s\n',
                formatC(sum(is_fastp_insert_size_peak_available), format = 'd', big.mark = ',')))
}

normalize_hist_breaks = function(values, bin_breaks = NULL, bins = 40) {
    if (!is.null(bin_breaks)) {
        out_breaks = sort(unique(as.numeric(bin_breaks)))
        if (length(out_breaks) >= 2) {
            return(out_breaks)
        }
    }
    values = values[is.finite(values)]
    if (length(values) == 0) {
        return(NULL)
    }
    value_min = min(values)
    value_max = max(values)
    if (value_min == value_max) {
        # Widen the range when all values are identical so multiple bins are shown.
        pad = max(1, abs(value_min) * 0.1)
        value_min = value_min - pad
        value_max = value_max + pad
    }
    out_breaks = pretty(c(value_min, value_max), n = max(3, bins))
    out_breaks = sort(unique(out_breaks))
    if (length(out_breaks) < 2) {
        step = max(1, ceiling(abs(value_min) * 0.05))
        center = mean(c(value_min, value_max))
        out_breaks = seq(center - 3 * step, center + 3 * step, by = step)
    }
    return(out_breaks)
}

save_frequency_plot = function(df, value_col, out_path, x_label, font_size = 8, bins = 40,
                                bin_breaks = NULL, x_breaks = NULL, x_limits = NULL, x_angle = 90) {
    font_size = resolve_plot_font_size(font_size)
    font_family = resolve_plot_font_family('Helvetica')
    if (!(value_col %in% colnames(df))) {
        cat(sprintf('%s column not found. Skipping plot generation.\n', value_col))
        return(NULL)
    }
    is_non_excluded = is_non_excluded_flag(df[['exclusion']])
    is_available = (!is.na(df[[value_col]]))
    df2 = df[(is_non_excluded & is_available), , drop = FALSE]
    cat(sprintf('Number of SRA samples for %s plotting: %s\n', value_col, formatC(nrow(df2), format = 'd', big.mark = ',')))
    if (nrow(df2) == 0) {
        cat(sprintf('No data available for %s. Skipping plot generation.\n', value_col))
        return(NULL)
    }
    df2[['scientific_name']] = as.character(df2[['scientific_name']])
    df2[['scientific_name']][is.na(df2[['scientific_name']]) | (df2[['scientific_name']] == '')] = 'Unknown species'
    df2[['scientific_name']] = factor(df2[['scientific_name']], levels = sort(unique(df2[['scientific_name']])))
    bin_breaks = normalize_hist_breaks(df2[[value_col]], bin_breaks = bin_breaks, bins = bins)
    g = ggplot(data = df2)
    if (is.null(bin_breaks)) {
        g = g + geom_histogram(
            aes(x = .data[[value_col]], fill = scientific_name),
            bins = bins,
            color = 'black',
            linewidth = 0.2,
            position = 'stack'
        )
    } else {
        g = g + geom_histogram(
            aes(x = .data[[value_col]], fill = scientific_name),
            breaks = bin_breaks,
            color = 'black',
            linewidth = 0.2,
            position = 'stack'
        )
    }
    if (!is.null(x_breaks)) {
        g = g + scale_x_continuous(breaks = x_breaks)
    }
    if (!is.null(x_limits)) {
        g = g + coord_cartesian(xlim = x_limits)
    }
    g = g + theme_bw(base_size = font_size, base_family = font_family)
    x_hjust = ifelse(x_angle == 0, 0.5, 1.0)
    g = g + labs(x = x_label, y = 'Sample count', fill = 'Species')
    g = g + build_standard_ggplot_theme(
        font_size = font_size,
        font_family = font_family,
        x_angle = x_angle,
        x_hjust = x_hjust,
        legend_position = 'bottom'
    ) + theme(
        legend.key.height = unit(0.24, 'cm'),
        legend.key.width = unit(0.35, 'cm')
    )
    ggsave(out_path, plot = g, width = 3.6, height = 3.6, units = 'in')
}

read_est_count_summary = function(file_path) {
    dat = tryCatch(
        read.table(file_path, sep = '\t', header = TRUE, quote = '', comment.char = '', check.names = FALSE),
        error = function(e) {
            warning(sprintf('Failed to read est_counts file. Skipping %s: %s', file_path, e$message))
            return(NULL)
        }
    )
    if (is.null(dat)) {
        return(NULL)
    }
    sample_cols = colnames(dat)
    if ('target_id' %in% sample_cols) {
        sample_cols = sample_cols[sample_cols != 'target_id']
    }
    if (length(sample_cols) == 0) {
        warning(sprintf('No sample columns found in est_counts file. Skipping %s', file_path))
        return(NULL)
    }
    dat_numeric = dat[, sample_cols, drop = FALSE]
    for (col in sample_cols) {
        dat_numeric[[col]] = suppressWarnings(as.numeric(dat_numeric[[col]]))
    }
    mean_before = vapply(dat_numeric, function(x) mean(x, na.rm = TRUE), numeric(1))
    lib_size = vapply(dat_numeric, function(x) sum(x, na.rm = TRUE), numeric(1))
    species_tag = basename(dirname(file_path))
    species_name = gsub('_', ' ', species_tag)
    summary_df = data.frame(
        scientific_name = rep(species_name, length(sample_cols)),
        run = sample_cols,
        mean_before = as.numeric(mean_before),
        lib_size = as.numeric(lib_size),
        stringsAsFactors = FALSE
    )
    return(summary_df)
}

collect_est_count_summaries = function(dir_merge) {
    est_count_files = sort(list.files(
        path = dir_merge,
        pattern = '_est_counts\\.tsv$',
        recursive = TRUE,
        full.names = TRUE
    ))
    if (length(est_count_files) == 0) {
        return(data.frame(
            scientific_name = character(0),
            run = character(0),
            mean_before = numeric(0),
            lib_size = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    summaries = vector('list', length(est_count_files))
    keep_index = 0
    for (file_path in est_count_files) {
        tmp = read_est_count_summary(file_path)
        if (is.null(tmp) || (nrow(tmp) == 0)) {
            next
        }
        keep_index = keep_index + 1
        summaries[[keep_index]] = tmp
    }
    if (keep_index == 0) {
        return(data.frame(
            scientific_name = character(0),
            run = character(0),
            mean_before = numeric(0),
            lib_size = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    out = do.call(rbind, summaries[seq_len(keep_index)])
    row_key = paste(out[['scientific_name']], out[['run']], sep = '|||')
    is_duplicated = duplicated(row_key)
    if (any(is_duplicated)) {
        duplicated_keys = sort(unique(row_key[is_duplicated]))
        warning(sprintf('Detected duplicated species/run rows across est_counts tables; keeping first entries: %s',
                        paste(duplicated_keys, collapse = ', ')))
        out = out[!is_duplicated, , drop = FALSE]
    }
    return(out)
}

save_mean_expression_boxplot = function(df, dir_merge, font_size = 8) {
    font_size = resolve_plot_font_size(font_size)
    font_family = resolve_plot_font_family('Helvetica')
    summary_df = collect_est_count_summaries(dir_merge = dir_merge)
    if (nrow(summary_df) == 0) {
        cat('No est_counts tables were found for mean expression plotting. Skipping merge_mean_expression_boxplot.pdf\n')
        return(NULL)
    }
    if (all(c('scientific_name', 'run') %in% colnames(df))) {
        metadata_cols = intersect(c('scientific_name', 'run', 'exclusion'), colnames(df))
        metadata_df = unique(df[, metadata_cols, drop = FALSE])
        summary_df = merge(
            summary_df,
            metadata_df,
            by = c('scientific_name', 'run'),
            sort = FALSE,
            all.x = TRUE,
            all.y = FALSE
        )
    } else if ('run' %in% colnames(df)) {
        metadata_cols = intersect(c('run', 'exclusion'), colnames(df))
        metadata_df = unique(df[, metadata_cols, drop = FALSE])
        summary_df = merge(summary_df, metadata_df, by = 'run', sort = FALSE, all.x = TRUE, all.y = FALSE)
    } else {
        summary_df[['exclusion']] = NA_character_
    }
    if ('exclusion' %in% colnames(summary_df)) {
        has_metadata_match = !is.na(summary_df[['exclusion']])
        if (sum(has_metadata_match) == 0) {
            cat('No est_counts runs matched metadata rows. Skipping merge_mean_expression_boxplot.pdf\n')
            return(NULL)
        }
        exclusion_norm = normalize_exclusion_values(summary_df[['exclusion']])
        is_non_excluded = has_metadata_match & (!is.na(exclusion_norm)) & (exclusion_norm == 'no')
        summary_df = summary_df[is_non_excluded, , drop = FALSE]
    }
    cat(sprintf('Number of SRA samples for mean_expression plotting: %s\n',
                formatC(nrow(summary_df), format = 'd', big.mark = ',')))
    if (nrow(summary_df) == 0) {
        cat('No non-excluded samples available for mean expression plotting. Skipping merge_mean_expression_boxplot.pdf\n')
        return(NULL)
    }
    valid_libsize = summary_df[['lib_size']]
    valid_libsize = valid_libsize[is.finite(valid_libsize) & (valid_libsize > 0)]
    if (length(valid_libsize) == 0) {
        scale_factor = rep(1.0, nrow(summary_df))
    } else {
        ref_libsize = stats::median(valid_libsize)
        if (!is.finite(ref_libsize) || (ref_libsize <= 0)) {
            ref_libsize = 1.0
        }
        scale_factor = summary_df[['lib_size']] / ref_libsize
        bad_scale = (!is.finite(scale_factor)) | (scale_factor <= 0)
        scale_factor[bad_scale] = 1.0
    }
    mean_before = as.numeric(summary_df[['mean_before']])
    mean_after = mean_before / scale_factor
    var_before = round(var(mean_before, na.rm = TRUE), digits = 1)
    var_after = round(var(mean_after, na.rm = TRUE), digits = 1)
    cat('Among-sample variance of mean all-gene raw counts before and after library-size normalization:',
        var_before, 'and', var_after, '\n')
    values = c(mean_before, mean_after)
    labels = c(
        rep('Raw\ncounts', length(mean_before)),
        rep('Library-size\ncorrected\ncounts', length(mean_after))
    )
    plot_df = data.frame(labels = labels, values = values, stringsAsFactors = FALSE)
    finite_values = plot_df[['values']]
    finite_values = finite_values[is.finite(finite_values)]
    finite_values = finite_values[finite_values >= 0]
    if (length(finite_values) == 0) {
        ymax = 1
    } else {
        ymax = max(finite_values)
        if ((!is.finite(ymax)) || (ymax <= 0)) {
            ymax = 1
        } else {
            ymax = ymax * 1.05
        }
    }
    g = ggplot(plot_df, aes(x = labels, y = values)) +
        geom_boxplot(outlier.shape = NA) +
        coord_cartesian(ylim = c(0, ymax)) +
        theme_bw(base_size = font_size, base_family = font_family) +
        labs(x = '', y = 'Mean count of all genes') +
        build_standard_ggplot_theme(font_size = font_size, font_family = font_family)
    out_path = file.path(dir_merge, 'merge_mean_expression_boxplot.pdf')
    ggsave(out_path, plot = g, width = 2.6, height = 3.6, units = 'in')
}

insert_axis_breaks = function(values) {
    values = values[is.finite(values)]
    if (length(values) == 0) {
        return(NULL)
    }
    value_min = floor(min(values, na.rm = TRUE))
    value_max = ceiling(max(values, na.rm = TRUE))
    value_range = value_max - value_min
    if (value_range <= 0) {
        by = ifelse(value_max <= 250, 10, 25)
    } else if (value_range <= 200) {
        by = 25
    } else if (value_range <= 500) {
        by = 50
    } else if (value_range <= 1000) {
        by = 100
    } else {
        by = 200
    }
    min_tick = floor(value_min / by) * by
    max_tick = ceiling(value_max / by) * by
    if (min_tick == max_tick) {
        min_tick = min_tick - 2 * by
        max_tick = max_tick + 2 * by
    }
    ticks = seq(min_tick, max_tick, by = by)
    if (length(ticks) < 3) {
        ticks = seq(min_tick - by, max_tick + by, by = by)
    }
    return(ticks)
}

df2 = df[(is_non_excluded & is_mapping_rate_available), , drop = FALSE]
cat(sprintf('Number of SRA samples for mapping_rate potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = mapping_rate), outlier.size = 0.3)
g = g + scale_x_discrete(labels = format_genus_species_label)
g = g + ylim(0, 100)
g = g + theme_bw(base_size = font_size, base_family = font_family)
g = g + labs(x = '', y = 'Mapping rate')
g = g + build_standard_ggplot_theme(
    font_size = font_size,
    font_family = font_family,
    x_angle = 90,
    legend_title_blank = TRUE
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_mapping_rate.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[is_non_excluded, , drop = FALSE]
cat(sprintf('Number of SRA samples for total_spots potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = total_spots), outlier.size = 0.3)
g = g + scale_x_discrete(labels = format_genus_species_label)
g = g + coord_cartesian(ylim = c(0, NA))
g = g + theme_bw(base_size = font_size, base_family = font_family)
g = g + labs(x = '', y = 'Total spots')
g = g + build_standard_ggplot_theme(
    font_size = font_size,
    font_family = font_family,
    x_angle = 90,
    legend_title_blank = TRUE
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_total_spots.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[is_non_excluded, , drop = FALSE]
cat(sprintf('Number of SRA samples for total_bases potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = total_bases), outlier.size = 0.3)
g = g + scale_x_discrete(labels = format_genus_species_label)
g = g + coord_cartesian(ylim = c(0, NA))
g = g + theme_bw(base_size = font_size, base_family = font_family)
g = g + labs(x = '', y = 'Total bases')
g = g + build_standard_ggplot_theme(
    font_size = font_size,
    font_family = font_family,
    x_angle = 90,
    legend_title_blank = TRUE
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_total_bases.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[is_non_excluded, , drop = FALSE]
cat(sprintf('Number of SRA samples for lib_layout potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
has_scientific_name = !is.na(df2[['scientific_name']]) & (trimws(as.character(df2[['scientific_name']])) != '')
has_lib_layout = !is.na(df2[['lib_layout']]) & (trimws(as.character(df2[['lib_layout']])) != '')
df2_layout = df2[(has_scientific_name & has_lib_layout), , drop = FALSE]
if (nrow(df2_layout) == 0) {
    cat('No non-excluded samples with scientific_name/lib_layout values. Skipping merge_library_layout.pdf\n')
} else {
    data_summary = aggregate(cbind(count = lib_layout) ~ scientific_name + lib_layout, df2_layout, length)
    data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN = sum)
    data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
    g = ggplot(data_summary, aes(x = scientific_name, y = count, fill = lib_layout))
    g = g + geom_bar(stat = "identity")
    g = g + scale_x_discrete(labels = format_genus_species_label)
    g = g + labs(x = "", y = "Count", fill = "Library layout")
    g = g + theme_bw(base_size = font_size, base_family = font_family)
    g = g + build_standard_ggplot_theme(
        font_size = font_size,
        font_family = font_family,
        x_angle = 90,
        legend_position = 'bottom',
        legend_title_blank = TRUE
    )
    num_spp = length(unique(df[['scientific_name']]))
    plot_width = max(3.6, 0.11 * num_spp)
    out_path = file.path(dir_merge, 'merge_library_layout.pdf')
    ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')
}

save_mean_expression_boxplot(
    df = df,
    dir_merge = dir_merge,
    font_size = font_size
)

save_frequency_plot(
    df = df,
    value_col = 'fastp_duplication_rate',
    out_path = file.path(dir_merge, 'merge_fastp_duplication_rate_histogram.pdf'),
    x_label = 'Fastp duplication rate (%)',
    font_size = font_size,
    bin_breaks = seq(0, 100, by = 10),
    x_breaks = seq(0, 100, by = 10),
    x_limits = c(0, 100),
    x_angle = 0
)

if ('fastp_insert_size_peak' %in% colnames(df)) {
    insert_df = df[(is_non_excluded & (!is.na(df[['fastp_insert_size_peak']]))), , drop = FALSE]
    insert_breaks = if (nrow(insert_df) > 0) insert_axis_breaks(insert_df[['fastp_insert_size_peak']]) else NULL
    insert_limits = if (is.null(insert_breaks)) NULL else c(min(insert_breaks), max(insert_breaks))
} else {
    insert_breaks = NULL
    insert_limits = NULL
}
save_frequency_plot(
    df = df,
    value_col = 'fastp_insert_size_peak',
    out_path = file.path(dir_merge, 'merge_fastp_insert_size_peak_histogram.pdf'),
    x_label = 'Fastp insert size peak',
    font_size = font_size,
    bin_breaks = insert_breaks,
    x_breaks = insert_breaks,
    x_limits = insert_limits,
    x_angle = 0
)

cat(sprintf('Number of SRA samples for exclusion potting: %s\n', formatC(nrow(df), format = 'd', big.mark = ',')))
out_path = file.path(dir_merge, 'merge_exclusion.pdf')
save_exclusion_plot(df = df, out_path = out_path, font_size = font_size, y_label = "Sample count")

cat('merge.r completed!\n')
