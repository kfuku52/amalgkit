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
font_size = 8

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

save_frequency_plot = function(df, value_col, out_path, x_label, font_size = 8, bins = 40,
                                bin_breaks = NULL, x_breaks = NULL, x_limits = NULL, x_angle = 90) {
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
    g = ggplot(data = df2)
    if (is.null(bin_breaks)) {
        g = g + geom_histogram(aes(x = .data[[value_col]]), bins = bins, color = 'black', fill = 'gray70')
    } else {
        g = g + geom_histogram(aes(x = .data[[value_col]]), breaks = bin_breaks, color = 'black', fill = 'gray70')
    }
    if (!is.null(x_breaks)) {
        g = g + scale_x_continuous(breaks = x_breaks)
    }
    if (!is.null(x_limits)) {
        g = g + coord_cartesian(xlim = x_limits)
    }
    g = g + theme_bw(base_size = font_size)
    x_hjust = ifelse(x_angle == 0, 0.5, 1.0)
    g = g + labs(x = x_label, y = 'Frequency')
    g = g + theme(
        axis.text = element_text(size = font_size, color = 'black'),
        axis.text.x = element_text(angle = x_angle, hjust = x_hjust, vjust = 0.5),
        axis.title = element_text(size = font_size, color = 'black'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = font_size, color = 'black'),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(rep(0.1, 4), "cm")
    )
    ggsave(out_path, plot = g, width = 3.6, height = 3.6, units = 'in')
}

insert_axis_breaks = function(values) {
    value_min = floor(min(values, na.rm = TRUE))
    value_max = ceiling(max(values, na.rm = TRUE))
    value_range = value_max - value_min
    if (value_range <= 200) {
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
    return(seq(min_tick, max_tick, by = by))
}

df2 = df[((!is_excluded) & (is_mapping_rate_available)),]
cat(sprintf('Number of SRA samples for mapping_rate potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = mapping_rate), outlier.size = 0.3)
g = g + ylim(0, 100)
g = g + theme_bw(base_size = font_size)
g = g + labs(x = '', y = 'Mapping rate')
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
    rect = element_rect(fill = "transparent"),
    plot.margin = unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_mapping_rate.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for total_spots potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = total_spots), outlier.size = 0.3)
g = g + scale_y_log10()
g = g + theme_bw(base_size = font_size)
g = g + labs(x = '', y = 'Total spots')
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
    rect = element_rect(fill = "transparent"),
    plot.margin = unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_total_spots.pdf')
ggsave(out_path, plot = g, width = plot_width, height = 3.6, units = 'in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for total_bases potting: %s\n', formatC(nrow(df2), format = 'd', big.mark = ',')))
g = ggplot(data = df2)
g = g + geom_boxplot(aes(x = scientific_name, y = total_bases), outlier.size = 0.3)
g = g + scale_y_log10()
g = g + theme_bw(base_size = font_size)
g = g + labs(x = '', y = 'Total bases')
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
    rect = element_rect(fill = "transparent"),
    plot.margin = unit(rep(0.1, 4), "cm")
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
g = g + labs(x = "", y = "Count", fill = "Library layout")
g = g + theme_bw(base_size = font_size)
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
save_exclusion_plot(df = df, out_path = out_path, font_size = font_size)

cat('merge.r completed!\n')
