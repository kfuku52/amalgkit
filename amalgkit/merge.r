#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly=TRUE)))
mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')
if (mode=="debug") {
    dir_work = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
    setwd(dir_work)
    dir_merge = file.path(dir_work, 'merge')
    file_metadata = file.path(dir_merge, 'metadata.tsv')
    r_util_path = '/Users/kf/Dropbox/repos/amalgkit/amalgkit/util.r'
} else if (mode=="batch") {
    args = commandArgs(trailingOnly=TRUE)
    dir_merge = args[1]
    file_metadata = args[2]
    r_util_path = args[3]
}
source(r_util_path)
font_size = 8

df = read.table(file_metadata, sep='\t', header=TRUE, quote='', comment.char='', check.names=FALSE)
is_excluded = (!df[['exclusion']]=='no')
is_mapping_rate_available = (!is.na(df[['mapping_rate']]))
cat(sprintf('Number of non-excluded SRA samples: %s\n', formatC(sum(!is_excluded), format='d', big.mark=',')))
cat(sprintf('Number of excluded SRA samples: %s\n', formatC(sum(is_excluded), format='d', big.mark=',')))
cat(sprintf('Number of SRA samples with available mapping rates: %s\n', formatC(sum(is_mapping_rate_available), format='d', big.mark=',')))

df2 = df[((!is_excluded)&(is_mapping_rate_available)),]
cat(sprintf('Number of SRA samples for mapping_rate potting: %s\n', formatC(nrow(df2), format='d', big.mark=',')))
g = ggplot(data=df2)
g = g + geom_boxplot(aes(x=scientific_name, y=mapping_rate), outlier.size=0.3)
g = g + ylim(0, 100)
g = g + theme_bw(base_size=font_size)
g = g + labs(x='', y='Mapping rate')
g = g + theme(
    axis.text=element_text(size=font_size, color='black'),
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title=element_text(size=font_size, color='black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    legend.title=element_blank(),
    legend.text=element_text(size=font_size, color='black'),
    rect=element_rect(fill="transparent"),
    plot.margin=unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_mapping_rate.pdf')
ggsave(out_path, plot=g, width=plot_width, height=3.6, units='in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for total_spots potting: %s\n', formatC(nrow(df2), format='d', big.mark=',')))
g = ggplot(data=df2)
g = g + geom_boxplot(aes(x=scientific_name, y=total_spots), outlier.size=0.3)
g = g + scale_y_log10()
g = g + theme_bw(base_size=font_size)
g = g + labs(x='', y='Total spots')
g = g + theme(
    axis.text=element_text(size=font_size, color='black'),
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title=element_text(size=font_size, color='black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    legend.title=element_blank(),
    legend.text=element_text(size=font_size, color='black'),
    rect=element_rect(fill="transparent"),
    plot.margin=unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_total_spots.pdf')
ggsave(out_path, plot=g, width=plot_width, height=3.6, units='in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for total_bases potting: %s\n', formatC(nrow(df2), format='d', big.mark=',')))
g = ggplot(data=df2)
g = g + geom_boxplot(aes(x=scientific_name, y=total_bases), outlier.size=0.3)
g = g + scale_y_log10()
g = g + theme_bw(base_size=font_size)
g = g + labs(x='', y='Total bases')
g = g + theme(
    axis.text=element_text(size=font_size, color='black'),
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title=element_text(size=font_size, color='black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    legend.title=element_blank(),
    legend.text=element_text(size=font_size, color='black'),
    rect=element_rect(fill="transparent"),
    plot.margin=unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_total_bases.pdf')
ggsave(out_path, plot=g, width=plot_width, height=3.6, units='in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for lib_layout potting: %s\n', formatC(nrow(df2), format='d', big.mark=',')))
data_summary = aggregate(cbind(count=lib_layout) ~ scientific_name + lib_layout, df2, length)
data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN=sum)
data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
g = ggplot(data_summary, aes(x=scientific_name, y=count, fill=lib_layout))
g = g + geom_bar(stat = "identity")
g = g + labs(x = "", y = "Count", fill = "Library layout")
g = g + theme_bw(base_size=font_size)
g = g + theme(
    axis.text=element_text(size=font_size, color='black'),
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title=element_text(size=font_size, color='black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    legend.title=element_blank(),
    legend.text=element_text(size=font_size, color='black'),
    legend.position='bottom',
    rect=element_rect(fill="transparent"),
    plot.margin=unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_library_layout.pdf')
ggsave(out_path, plot=g, width=plot_width, height=3.6, units='in')

df2 = df[((!is_excluded)),]
cat(sprintf('Number of SRA samples for curate_group potting: %s\n', formatC(nrow(df2), format='d', big.mark=',')))
data_summary = aggregate(cbind(count=curate_group) ~ scientific_name + curate_group, df2, length)
data_summary[['total']] = ave(data_summary[['count']], data_summary[['scientific_name']], FUN=sum)
data_summary[['proportion']] = data_summary[['count']] / data_summary[['total']]
g = ggplot(data_summary, aes(x=scientific_name, y=count, fill=curate_group))
g = g + geom_bar(stat = "identity")
g = g + labs(x = "", y = "Count", fill = "curate_group")
g = g + theme_bw(base_size=font_size)
g = g + theme(
    axis.text=element_text(size=font_size, color='black'),
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    axis.title=element_text(size=font_size, color='black'),
    #panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    legend.title=element_blank(),
    legend.text=element_text(size=font_size, color='black'),
    legend.position='bottom',
    rect=element_rect(fill="transparent"),
    plot.margin=unit(rep(0.1, 4), "cm")
)
num_spp = length(unique(df[['scientific_name']]))
plot_width = max(3.6, 0.11 * num_spp)
out_path = file.path(dir_merge, 'merge_curate_group.pdf')
ggsave(out_path, plot=g, width=plot_width, height=3.6, units='in')

cat('merge.r completed!\n')
