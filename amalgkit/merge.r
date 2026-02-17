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

df = read.table(file_metadata, sep = '\t', header = TRUE, quote = '', comment.char = '', check.names = FALSE)
is_excluded = (!df[['exclusion']] == 'no')
is_mapping_rate_available = (!is.na(df[['mapping_rate']]))
cat(sprintf('Number of non-excluded SRA samples: %s\n', formatC(sum(!is_excluded), format = 'd', big.mark = ',')))
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
    is_excluded = (!df[['exclusion']] == 'no')
    is_available = (!is.na(df[[value_col]]))
    df2 = df[((!is_excluded) & (is_available)),]
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

df2 = df[((!is_excluded) & (is_mapping_rate_available)),]
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

df2 = df[((!is_excluded)),]
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

df2 = df[((!is_excluded)),]
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

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for lib_layout potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
data_summary = aggregate(cbind(count = lib_layout) ~ scientific_name + lib_layout, df2, length)
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
    insert_df = df[((!is_excluded) & (!is.na(df[['fastp_insert_size_peak']]))),]
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
