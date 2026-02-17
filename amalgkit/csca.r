#!/usr/bin/env Rscript

detected_cores = tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
if (is.na(detected_cores) && (is.null(getOption("cores")) || is.na(getOption("cores")))) {
    options(cores = 1L)
}
suppressWarnings(suppressPackageStartupMessages(library(Rtsne, quietly = TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quietly = TRUE)))
options(stringsAsFactors = FALSE)
correction_labels = c('uncorrected', 'corrected')

save_ggplot_grid = function(plots, filename, width, height, nrow, ncol, font_size = 8, font_family = 'Helvetica') {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family(font_family)
    with_pdf_defaults(
        file = filename,
        width = width,
        height = height,
        font_size = font_size,
        font_family = font_family,
        plot_fn = function(local_font_size, local_font_family) {
            grid::grid.newpage()
            grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow, ncol = ncol)))
            for (idx in seq_along(plots)) {
                row = ((idx - 1) %/% ncol) + 1
                col = ((idx - 1) %% ncol) + 1
                print(plots[[idx]], vp = grid::viewport(layout.pos.row = row, layout.pos.col = col))
            }
        }
    )
}

debug_mode = ifelse(length(commandArgs(trailingOnly = TRUE)) == 1, "debug", "batch")
font_size = 8
species_shape_threshold = 12

if (debug_mode == "debug") {
    developer = 'kf'
    if (developer == 'mf') {
        selected_sample_groups = c('root', 'flower', 'leaf')
        dir_work = '/home/s229181/projects/amalgkit_paper/Plant_set'
        #selected_sample_groups = c('Amicoumacin_low','Amicoumacin_high','anaerobic','Azithromycin_low','Azithromycin_high','control','ciprofloxacin_low','ciprofloxacin_high','colistin_low','colistin_high','H2O2','heat','meropenem_low','meropenem_high','NaCl','Anaerobic','IPTG','biofilm_medium','acid','H202','butyrate','aerobic','Ampicillin','Vancomycin','ciprofloxacin','colistin','glucose','Glucose','rifampicin','probiotic','cleaner','ampicillin','tetracycline','pediocin','glycerol','Pyruvate','Ca2+','Glycerol','H2o2','anhydrotetracycline','TB47','treated','iron','Lactate','rion','phage','Ag','biofilm','MIC','AZM','citrate','NaNO2','Acetate','sucrose','coumermycin','copper','mitomycin','arabinose','Cefotaxime','Cellulose','vancomycin','mupirocin','galactose','macrophages','tobramycin')
        sample_group_colors = 'DEFAULT'
        #dir_work = '/home/s229181/projects/amalgkit_paper/Prokaryote_set'
        dir_csca_input_table = file.path(dir_work, 'csca/csca_input_symlinks')
        file_orthogroup = file.path(dir_work, 'csca/multispecies_busco_table.tsv')
        file_genecount = file.path(dir_work, 'csca/multispecies_genecount.tsv')
        r_util_path = '/home/s229181/projects/amalgkit_paper/amalgkit/amalgkit/util.r'
        dir_csca = file.path(dir_work, 'csca')
        batch_effect_alg = 'sva'
        missing_strategy = 'em_pca'
    } else if (developer == 'kf') {
        selected_sample_groups = c('root', 'flower', 'leaf')
        sample_group_colors = c('#d95f02ff', '#1b9e77ff', '#7570b3ff')
        dir_work = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/data/evolutionary_transcriptomics/20230527_amalgkit/amalgkit_out'
        dir_csca_input_table = file.path(dir_work, 'csca/csca_input_symlinks')
        file_orthogroup = file.path(dir_work, 'csca/multispecies_busco_table.tsv')
        file_genecount = file.path(dir_work, 'csca/multispecies_genecount.tsv')
        r_util_path = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/repos/amalgkit/amalgkit/util.r'
        dir_csca = file.path(dir_work, 'csca')
        batch_effect_alg = 'sva'
        missing_strategy = 'em_pca'
    }
} else if (debug_mode == "batch") {
    args = commandArgs(trailingOnly = TRUE)
    selected_sample_groups = strsplit(args[1], "\\|")[[1]]
    sample_group_colors = strsplit(args[2], ",")[[1]]
    dir_work = args[3]
    dir_csca_input_table = args[4]
    file_orthogroup = args[5]
    file_genecount = args[6]
    r_util_path = args[7]
    dir_csca = args[8]
    batch_effect_alg = args[9]
    if (length(args) >= 10) {
        missing_strategy = args[10]
    } else {
        missing_strategy = 'em_pca'
    }
}
source(r_util_path)
resolve_csca_font_size = function(font_size = PLOT_FONT_SIZE_PT) {
    resolve_plot_font_size(font_size)
}
resolve_csca_font_family = function(font_family = PLOT_FONT_FAMILY) {
    resolve_plot_font_family(font_family)
}
compute_dense_label_cex = function(num_labels, base_pt = PLOT_FONT_SIZE_PT, min_pt = 4, soft_limit = 20) {
    dense_pt = compute_dense_label_pt(
        num_labels = num_labels,
        base_pt = base_pt,
        min_pt = min_pt,
        soft_limit = soft_limit
    )
    dense_pt / base_pt
}
font_size = resolve_csca_font_size(font_size)
font_family = resolve_csca_font_family('Helvetica')
setwd(dir_csca)
cat('selected_sample_groups:', selected_sample_groups, "\n")
cat('selected_sample_group_colors:', sample_group_colors, "\n")
cat('dir_work:', dir_work, "\n")
cat('dir_csca_input_table:', dir_csca_input_table, "\n")
cat('file_orthogroup:', file_orthogroup, "\n")
cat('file_genecount:', file_genecount, "\n")
cat('r_util_path:', r_util_path, "\n")
cat('dir_csca:', dir_csca, "\n")
cat('batch_effect_alg:', batch_effect_alg, "\n")
cat('missing_strategy:', missing_strategy, "\n")

sort_labels = function(df_label, label_orders) {
    if ((nrow(df_label) == 0) || (length(label_orders) == 0)) {
        return(df_label[0, , drop = FALSE])
    }
    label_key = sub(' ', '_', paste(df_label[['scientific_name']], df_label[['sample_group']], sep = '_'))
    row_idx_by_key = split(seq_len(nrow(df_label)), label_key)
    ordered_rows = unlist(row_idx_by_key[label_orders], use.names = FALSE)
    if (length(ordered_rows) == 0) {
        return(df_label[0, , drop = FALSE])
    }
    df_label[ordered_rows, , drop = FALSE]
}

propagate_uniform_edge_colors = function(node) {
    if (is.leaf(node) || (length(node) == 0)) {
        return(node)
    }
    for (i in seq_along(node)) {
        node[[i]] = propagate_uniform_edge_colors(node[[i]])
    }
    if (length(node) != 2) {
        return(node)
    }
    child1_color = attr(node[[1]], "edgePar")[['col']]
    child2_color = attr(node[[2]], "edgePar")[['col']]
    if (is.null(child1_color) || is.null(child2_color)) {
        return(node)
    }
    if (is.na(child1_color) || is.na(child2_color)) {
        return(node)
    }
    if (identical(child1_color, child2_color)) {
        edge_par = attr(node, "edgePar")
        if (is.null(edge_par)) {
            edge_par = list()
        }
        edge_par[['col']] = child1_color
        attr(node, "edgePar") = edge_par
    }
    return(node)
}

distance_from_correlation_matrix = function(cor_mat, na_fill = 1, epsilon = 0) {
    cor_mat = as.matrix(cor_mat)
    dist_mat = 1 - cor_mat
    dist_mat[!is.finite(dist_mat)] = na_fill
    dist_mat = (dist_mat + t(dist_mat)) / 2
    dist_mat = dist_mat + epsilon
    diag(dist_mat) = 0
    return(as.dist(dist_mat))
}

build_tc_plot_cache = function(tc) {
    tc_cor_plot = as.matrix(suppressWarnings(cor(tc, method = 'pearson')))
    tc_cor_plot[is.na(tc_cor_plot)] = 0
    tc_nonnegative = tc
    tc_nonnegative[tc_nonnegative < 0] = 0
    tc_cor_boxplot = as.matrix(suppressWarnings(cor(tc_nonnegative, method = 'pearson')))
    tc_cor_boxplot[is.na(tc_cor_boxplot)] = 0
    tc_cor_pairwise = as.matrix(suppressWarnings(cor(tc, method = 'pearson', use = 'pairwise.complete.obs')))
    list(
        tc_cor_plot = tc_cor_plot,
        tc_cor_boxplot = tc_cor_boxplot,
        tc_dist_dendrogram = distance_from_correlation_matrix(tc_cor_pairwise, na_fill = 1, epsilon = 0),
        tc_dist_mds = distance_from_correlation_matrix(tc_cor_pairwise, na_fill = 1, epsilon = 0.000000001)
    )
}

build_averaged_plot_cache = function(averaged_orthologs) {
    cache = list()
    for (d in correction_labels) {
        cache[[d]] = build_tc_plot_cache(averaged_orthologs[[d]])
    }
    return(cache)
}

get_cached_plot_data = function(averaged_plot_cache, correction, cache_key, fallback = NULL) {
    if (is.null(averaged_plot_cache)) {
        return(fallback)
    }
    correction_cache = averaged_plot_cache[[correction]]
    if (is.null(correction_cache)) {
        return(fallback)
    }
    cached_value = correction_cache[[cache_key]]
    if (is.null(cached_value)) {
        return(fallback)
    }
    cached_value
}

for_each_averaged_correction = function(averaged_orthologs, df_color_averaged, fn) {
    for (d in correction_labels) {
        fn(
            correction = d,
            tc = averaged_orthologs[[d]],
            df_label = df_color_averaged
        )
    }
}

format_csca_heatmap_labels = function(labels) {
    labels = as.character(labels)
    sapply(labels, function(label) {
        tokens = strsplit(label, "_", fixed = TRUE)[[1]]
        if (length(tokens) >= 3) {
            genus_species = paste(tokens[1], tokens[2], sep = " ")
            sample_group = paste(tokens[3:length(tokens)], collapse = "_")
            return(paste(genus_species, sample_group, sep = "\n"))
        }
        gsub("_", " ", label)
    })
}

get_species_shape_map = function(species_values) {
    species_values = as.character(species_values)
    species_unique = unique(species_values[!is.na(species_values)])
    # Use a broad pool so species remain distinguishable in larger datasets.
    shape_pool = c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14)
    if (length(species_unique) == 0) {
        return(setNames(integer(0), character(0)))
    }
    shape_values = rep(shape_pool, length.out = length(species_unique))
    names(shape_values) = species_unique
    return(shape_values)
}

get_point_shapes_for_species = function(species_values, species_shape_map = NULL, default_shape = 16) {
    species_values = as.character(species_values)
    if (is.null(species_shape_map)) {
        species_shape_map = get_species_shape_map(species_values)
    }
    pch_values = unname(species_shape_map[species_values])
    pch_values[is.na(pch_values)] = default_shape
    return(pch_values)
}

get_species_encoding_info = function(species_values, shape_threshold = species_shape_threshold) {
    species_values = as.character(species_values)
    species_unique = unique(species_values[!is.na(species_values)])
    n_species = length(species_unique)
    list(
        species_unique = species_unique,
        n_species = n_species,
        use_shape = (n_species <= shape_threshold)
    )
}

get_species_centroids = function(x, y, species_values) {
    df = data.frame(
        x = as.numeric(x),
        y = as.numeric(y),
        scientific_name = as.character(species_values),
        stringsAsFactors = FALSE
    )
    df = df[is.finite(df[['x']]) & is.finite(df[['y']]) & !is.na(df[['scientific_name']]), , drop = FALSE]
    if (nrow(df) == 0) {
        return(data.frame(x = numeric(0), y = numeric(0), scientific_name = character(0), stringsAsFactors = FALSE))
    }
    aggregate(cbind(x, y) ~ scientific_name, data = df, FUN = median)
}

add_species_centroid_labels_base = function(x, y, species_values, cex = 1, col = 'black') {
    centroids = get_species_centroids(x = x, y = y, species_values = species_values)
    if (nrow(centroids) == 0) {
        return(invisible(NULL))
    }
    text(
        x = centroids[['x']],
        y = centroids[['y']],
        labels = centroids[['scientific_name']],
        cex = cex,
        col = col,
        pos = 3,
        offset = 0.2,
        xpd = TRUE
    )
    invisible(NULL)
}

add_species_centroid_labels_ggplot = function(g, df, x_col, y_col, font_size = 8) {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family('Helvetica')
    centroids = get_species_centroids(
        x = df[[x_col]],
        y = df[[y_col]],
        species_values = df[['scientific_name']]
    )
    if (nrow(centroids) == 0) {
        return(g)
    }
    g + ggplot2::geom_text(
        data = centroids,
        mapping = ggplot2::aes(x = x, y = y, label = scientific_name),
        inherit.aes = FALSE,
        family = font_family,
        size = font_size / ggplot2::.pt,
        color = 'black',
        check_overlap = TRUE,
        vjust = -0.4,
        show.legend = FALSE
    )
}

get_species_base_scatter_style = function(df_label, default_shape = 16, fallback_color = 'gray40') {
    species_info = get_species_encoding_info(df_label[['scientific_name']])
    point_color = as.character(df_label[['sample_group_color']])
    point_color[is.na(point_color)] = fallback_color
    point_shape = if (species_info[['use_shape']]) {
        get_point_shapes_for_species(df_label[['scientific_name']], default_shape = default_shape)
    } else {
        rep(default_shape, nrow(df_label))
    }
    return(list(
        species_info = species_info,
        point_color = point_color,
        point_shape = point_shape
    ))
}

draw_species_scatter_base = function(x, y, df_label, xlab, ylab, cex = 2, lwd = 1, las = 1) {
    style = get_species_base_scatter_style(df_label = df_label)
    plot(
        x,
        y,
        pch = style[['point_shape']],
        cex = cex,
        lwd = lwd,
        col = style[['point_color']],
        xlab = xlab,
        ylab = ylab,
        las = las
    )
    if (!style[['species_info']][['use_shape']]) {
        add_species_centroid_labels_base(
            x = x,
            y = y,
            species_values = df_label[['scientific_name']]
        )
    }
    invisible(style)
}

add_species_style_layers_ggplot = function(g, df, x_col, y_col, sample_group_colors,
                                           font_size = 8, point_size = 0.7, point_alpha = 1, na_rm = TRUE) {
    species_info = get_species_encoding_info(df[['scientific_name']])
    species_shapes = get_species_shape_map(df[['scientific_name']])
    g = g + ggplot2::scale_color_manual(values = sample_group_colors)
    if (species_info[['use_shape']]) {
        g = g +
            ggplot2::geom_point(ggplot2::aes(shape = scientific_name), size = point_size, alpha = point_alpha, na.rm = na_rm) +
            ggplot2::scale_shape_manual(values = species_shapes) +
            ggplot2::labs(color = 'Sample group', shape = 'Species') +
            ggplot2::guides(
                color = ggplot2::guide_legend(order = 1),
                shape = ggplot2::guide_legend(order = 2)
            )
    } else {
        g = g +
            ggplot2::geom_point(size = point_size, alpha = point_alpha, na.rm = na_rm) +
            ggplot2::labs(color = 'Sample group') +
            ggplot2::guides(color = ggplot2::guide_legend(order = 1))
        g = add_species_centroid_labels_ggplot(g, df, x_col = x_col, y_col = y_col, font_size = font_size)
    }
    return(g)
}

get_sample_group_palette = function(df, group_col = 'sample_group', color_col = 'sample_group_color',
                                     fallback_color = 'gray40') {
    tmp = unique(df[, c(group_col, color_col), drop = FALSE])
    tmp[[group_col]] = as.character(tmp[[group_col]])
    tmp[[color_col]] = as.character(tmp[[color_col]])
    tmp = tmp[!is.na(tmp[[group_col]]), , drop = FALSE]
    tmp[[color_col]][is.na(tmp[[color_col]])] = fallback_color
    palette = tmp[[color_col]]
    names(palette) = tmp[[group_col]]
    palette
}

compute_scatter_axis_limits = function(xvals, yvals, require_variation = FALSE, pad_ratio = 0.01, flat_half_span = 0.5) {
    finite_x = as.numeric(xvals[is.finite(xvals)])
    finite_y = as.numeric(yvals[is.finite(yvals)])
    if ((length(finite_x) == 0) || (length(finite_y) == 0)) {
        return(list(ok = FALSE, reason = 'no_finite_data'))
    }
    xmin = min(finite_x)
    xmax = max(finite_x)
    ymin = min(finite_y)
    ymax = max(finite_y)
    if (require_variation && ((xmin == xmax) || (ymin == ymax))) {
        return(list(ok = FALSE, reason = 'no_variation'))
    }
    if (xmin == xmax) {
        xmin = xmin - flat_half_span
        xmax = xmax + flat_half_span
    }
    if (ymin == ymax) {
        ymin = ymin - flat_half_span
        ymax = ymax + flat_half_span
    }
    xpad = (xmax - xmin) * pad_ratio
    ypad = (ymax - ymin) * pad_ratio
    return(list(
        ok = TRUE,
        xlim = c(xmin - xpad, xmax + xpad),
        ylim = c(ymin - ypad, ymax + ypad)
    ))
}

add_scatter_text_theme = function(g, font_size = 8) {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family('Helvetica')
    g + build_standard_ggplot_theme(font_size = font_size, font_family = font_family)
}

build_unaveraged_scatter_plot = function(df, x_col, y_col, x_label, y_label, axis_limits, sample_group_colors,
                                         font_size = 8, point_size = 0.7, point_alpha = 1, na_rm = TRUE) {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family('Helvetica')
    g = ggplot(df, ggplot2::aes(x = !!rlang::sym(x_col), y = !!rlang::sym(y_col), color = sample_group))
    g = g + ggplot2::theme_bw(base_size = font_size, base_family = font_family)
    g = add_species_style_layers_ggplot(
        g = g,
        df = df,
        x_col = x_col,
        y_col = y_col,
        sample_group_colors = sample_group_colors,
        font_size = font_size,
        point_size = point_size,
        point_alpha = point_alpha,
        na_rm = na_rm
    )
    g = g +
        ggplot2::xlab(x_label) +
        ggplot2::ylab(y_label) +
        ggplot2::coord_cartesian(xlim = axis_limits[['xlim']], ylim = axis_limits[['ylim']])
    add_scatter_text_theme(g, font_size = font_size)
}

save_scatter_plot_pdf = function(g, filename, fail_prefix = 'Plot could not be computed for file',
                                 width = 4.25, height = 2.15) {
    tryCatch(
        {
            ggplot2::ggsave(file = filename, g, height = height, width = width)
        },
        error = function(cond) {
            message(paste(fail_prefix, filename))
        }
    )
}

append_columns_by_key = function(left_df, right_df, key_col) {
    left_key = as.character(left_df[[key_col]])
    right_key = as.character(right_df[[key_col]])
    right_idx = match(left_key, right_key)
    value_cols = setdiff(colnames(right_df), key_col)
    for (col in value_cols) {
        left_df[[col]] = right_df[[col]][right_idx]
    }
    left_df
}

draw_empty_panel = function() {
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
}

compute_tsne_embedding = function(tc, seed = 1, max_perplexity = 30, min_samples = 4,
                                  dims = 2, theta = 0, check_duplicates = FALSE, verbose = FALSE) {
    num_samples = ncol(tc)
    perplexity = min(max_perplexity, floor(num_samples / 4), na.rm = TRUE)
    if ((num_samples < min_samples) || (!is.finite(perplexity)) || (perplexity <= 0)) {
        return(list(
            ok = FALSE,
            reason = 'insufficient_samples',
            num_samples = num_samples,
            perplexity = perplexity
        ))
    }
    set.seed(seed)
    out_tsne = try(
        Rtsne(
            as.matrix(t(tc)),
            theta = theta,
            check_duplicates = check_duplicates,
            verbose = verbose,
            perplexity = perplexity,
            dims = dims
        ),
        silent = TRUE
    )
    if ('try-error' %in% class(out_tsne)) {
        return(list(
            ok = FALSE,
            reason = 'tsne_failed',
            num_samples = num_samples,
            perplexity = perplexity
        ))
    }
    Y = as.matrix(out_tsne$Y)
    if (ncol(Y) < 2) {
        return(list(
            ok = FALSE,
            reason = 'insufficient_dimensions',
            num_samples = num_samples,
            perplexity = perplexity
        ))
    }
    return(list(
        ok = TRUE,
        Y = Y,
        num_samples = num_samples,
        perplexity = perplexity
    ))
}

draw_multisp_heatmap = function(tc, df_label, tc_dist_matrix = NULL, font_size = 8, show_colorbar = FALSE) {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family('Helvetica')
    if (is.null(tc_dist_matrix)) {
        tc_dist_matrix = as.matrix(suppressWarnings(cor(tc, method = 'pearson')))
    } else {
        tc_dist_matrix = as.matrix(tc_dist_matrix)
    }
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    n = ncol(tc_dist_matrix)
    sample_label_pt = compute_dense_label_pt(num_labels = n, base_pt = font_size, min_pt = 4, soft_limit = 20)
    if (n == 0) {
        g = ggplot2::ggplot() + ggplot2::theme_void(base_size = font_size, base_family = font_family)
        draw_ggplot_in_current_plot_panel(g)
        return(invisible(NULL))
    }
    row_names = rownames(tc_dist_matrix)
    col_names = colnames(tc_dist_matrix)
    if (is.null(row_names) || length(row_names) != n) {
        row_names = as.character(seq_len(n))
    }
    if (is.null(col_names) || length(col_names) != n) {
        col_names = as.character(seq_len(n))
    }
    row_label_text = format_csca_heatmap_labels(row_names)
    col_label_text = format_csca_heatmap_labels(col_names)
    label_key = paste0(gsub(' ', '_', as.character(df_label[['scientific_name']])), '_', as.character(df_label[['sample_group']]))
    label_index = match(col_names, label_key)
    sample_group_color = as.character(df_label[['sample_group_color']])[label_index]
    sp_color = as.character(df_label[['sp_color']])[label_index]
    sample_group_color[is.na(sample_group_color)] = 'gray70'
    sp_color[is.na(sp_color)] = 'gray40'
    rownames(tc_dist_matrix) = row_names
    colnames(tc_dist_matrix) = col_names
    df_heat = as.data.frame(as.table(tc_dist_matrix), stringsAsFactors = FALSE)
    colnames(df_heat) = c('row_label', 'col_label', 'value')
    df_heat[['row_idx']] = match(df_heat[['row_label']], row_names)
    df_heat[['col_idx']] = match(df_heat[['col_label']], col_names)
    df_heat[['x']] = df_heat[['col_idx']]
    df_heat[['y']] = n - df_heat[['row_idx']] + 1

    ann_gap = 0.25
    top_sample_group_y = n + 1 + ann_gap
    top_species_y = n + 2 + ann_gap
    left_sample_group_x = 0 - ann_gap
    left_species_x = -1 - ann_gap
    rdylbu11 = c(
        '#A50026', '#D73027', '#F46D43', '#FDAE61', '#FEE090', '#FFFFBF',
        '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#313695'
    )
    heat_breaks = c(0, seq(0.3, 1, 0.01))
    heat_colors = grDevices::colorRampPalette(rev(rdylbu11))(length(heat_breaks))
    ann_points = rbind(
        data.frame(
            x = seq_len(n),
            y = rep(top_sample_group_y, n),
            color = sample_group_color,
            stringsAsFactors = FALSE
        ),
        data.frame(
            x = seq_len(n),
            y = rep(top_species_y, n),
            color = sp_color,
            stringsAsFactors = FALSE
        ),
        data.frame(
            x = rep(left_sample_group_x, n),
            y = n - seq_len(n) + 1,
            color = sample_group_color,
            stringsAsFactors = FALSE
        ),
        data.frame(
            x = rep(left_species_x, n),
            y = n - seq_len(n) + 1,
            color = sp_color,
            stringsAsFactors = FALSE
        )
    )
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
            inherit.aes = FALSE,
            show.legend = FALSE
        ) +
        ggplot2::scale_fill_gradientn(
            colors = heat_colors,
            values = scales::rescale(heat_breaks, from = c(0, 1)),
            limits = c(0, 1),
            oob = scales::squish,
            na.value = 'gray80',
            name = "Pearson's\ncorrelation\ncoefficient",
            breaks = c(0, 0.5, 1.0),
            guide = if (show_colorbar) ggplot2::guide_colorbar(order = 1, nbin = 256) else 'none'
        ) +
        ggplot2::scale_x_continuous(
            breaks = c(left_species_x, left_sample_group_x, seq_len(n)),
            labels = c('Species', 'Sample group', col_label_text),
            expand = c(0, 0),
            position = 'top'
        ) +
        ggplot2::scale_y_continuous(
            breaks = c(seq_len(n), top_sample_group_y, top_species_y),
            labels = c(rev(row_label_text), 'Sample group', 'Species'),
            expand = c(0, 0)
        ) +
        ggplot2::coord_fixed(
            xlim = c(left_species_x - 0.5, n + 0.5),
            ylim = c(0.5, top_species_y + 0.5),
            clip = 'on'
        ) +
        ggplot2::theme_bw(base_size = font_size, base_family = font_family) +
        ggplot2::theme(
            text = ggplot2::element_text(size = font_size, family = font_family, color = 'black'),
            panel.grid = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = sample_label_pt, family = font_family, color = 'black', angle = 90, hjust = 0, vjust = 0.5),
            axis.text.y = ggplot2::element_text(size = sample_label_pt, family = font_family, color = 'black'),
            axis.ticks = ggplot2::element_blank(),
            legend.position = if (show_colorbar) 'right' else 'none',
            legend.title = ggplot2::element_text(size = font_size, family = font_family, color = 'black'),
            legend.text = ggplot2::element_text(size = font_size, family = font_family, color = 'black'),
            legend.key.height = grid::unit(0.14, 'in'),
            legend.spacing.y = grid::unit(0.04, 'in'),
            panel.border = ggplot2::element_blank(),
            plot.background = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(6, 6, 6, 6, unit = 'pt')
        )
    draw_ggplot_in_current_plot_panel(g)
    invisible(NULL)
}

draw_multisp_dendrogram = function(tc, df_label, df_metadata,
                                   label_cex = 1, axis_cex = 1,
                                   nboot = NULL, pvclust_file = NULL, tc_dist_dist = NULL) {
    if (!is.null(nboot) || !is.null(pvclust_file)) {
        cat('Bootstrap-based dendrogram support is deprecated. Using hclust without bootstrap.\n')
    }
    colnames(tc) = sub("_.*", "", sub('_', ' ', colnames(tc)))
    if (is.null(tc_dist_dist)) {
        tc_dist_dist = calc_sample_distance(tc, method = 'pearson')
    }
    attr(tc_dist_dist, 'Labels') = colnames(tc)
    hc = hclust(tc_dist_dist, method = 'average')
    dend = as.dendrogram(hc)
    dend_colors = df_label[order.dendrogram(dend), 'sample_group_color']
    label_colors = df_label[order.dendrogram(dend), 'sp_color']
    dend = apply_leaf_label_colors(dend, label_colors)
    dend = apply_leaf_edge_colors(dend, dend_colors)
    dend = propagate_uniform_edge_colors(dend)
    dend = set_edge_lwd(dend, lwd = 2)
    plot(dend, las = 1, yaxt = 'n', cex = label_cex)
    axis(2, las = 1, cex.axis = axis_cex)
    mtext(text = 'Distance', side = 2, line = 4, cex = axis_cex)

    n = ncol(tc)
    f = 100
    sample_group_unique = unique(df_metadata['sample_group'])
    sp_unique = unique(df_metadata[['scientific_name']])
    bp_unique = unique(df_metadata[['bioproject']])
    sample_group_color_unique = unique(df_metadata[['sample_group_color']])
    sp_color_unique = unique(df_metadata[['sp_color']])
    bp_color_unique = unique(df_metadata[['bp_color']])
    legend_text = c(as.character(sample_group_unique), "", as.character(sp_unique), "", as.character(bp_unique))
    legend_bg = c(sample_group_color_unique, "white", sp_color_unique, "white", bp_color_unique)
    legend_fg = c(rep("black", length(sample_group_color_unique)), "white", rep("black", length(sp_color_unique)), "white", rep("black", length(bp_color_unique)))
    #plot.new() ; par(mar=c(0,0,0,0))
    #legend("center", legend=legend_text, cex=1, pch=22, lty=0, lwd=1, pt.bg=legend_bg, col=legend_fg)
}

draw_multisp_pca = function(tc, df_label, tc_dist_matrix = NULL) {
    if (is.null(tc_dist_matrix)) {
        tc_dist_matrix = as.matrix(suppressWarnings(cor(tc, method = 'pearson')))
    } else {
        tc_dist_matrix = as.matrix(tc_dist_matrix)
    }
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    set.seed(1)
    pca = prcomp(tc_dist_matrix)
    xlabel = paste0("PC 1 (", round(summary(pca)$importance[2, 1] * 100, digits = 1), "%)")
    ylabel = paste0("PC 2 (", round(summary(pca)$importance[2, 2] * 100, digits = 1), "%)")
    draw_species_scatter_base(
        x = pca[['x']][, 1],
        y = pca[['x']][, 2],
        df_label = df_label,
        xlab = xlabel,
        ylab = ylabel,
        cex = 2,
        lwd = 1,
        las = 1
    )
}

draw_multisp_mds = function(tc, df_label, tc_dist_dist = NULL) {
    if (is.null(tc_dist_dist)) {
        tc_dist_dist = calc_sample_distance(tc, method = 'pearson', na_fill = 1, epsilon = 0.000000001)
    }
    set.seed(1)
    mds_points = tryCatch(
    { as.matrix(stats::cmdscale(tc_dist_dist, k = 2)) },
    error = function(a) { NULL }
    )
    if (is.null(mds_points) || (ncol(mds_points) < 2)) {
        cat('MDS failed.\n')
        draw_empty_panel()
    } else {
        draw_species_scatter_base(
            x = mds_points[, 1],
            y = mds_points[, 2],
            df_label = df_label,
            xlab = "MDS dimension 1",
            ylab = "MDS dimension 2",
            cex = 2,
            lwd = 1,
            las = 1
        )
    }
}

draw_multisp_tsne = function(tc, df_label) {
    tsne_result = compute_tsne_embedding(tc = tc)
    if (!tsne_result[['ok']]) {
        if (tsne_result[['reason']] == 'insufficient_samples') {
            cat('t-SNE skipped: insufficient number of samples.\n')
        } else {
            cat('t-SNE failed.\n')
        }
        draw_empty_panel()
        return(invisible(NULL))
    }
    try_out = tryCatch(
    {
        draw_species_scatter_base(
            x = tsne_result[['Y']][, 1],
            y = tsne_result[['Y']][, 2],
            df_label = df_label,
            xlab = "t-SNE dimension 1",
            ylab = "t-SNE dimension 2",
            cex = 2,
            lwd = 1,
            las = 1
        )
    },
        error = function(a) { return("t-SNE plot failed.") }
    )
    if (mode(try_out) == "character") {
        cat('t-SNE failed.\n')
        draw_empty_panel()
    }
}

draw_multisp_legend = function(df_label) {
    cex_axis = 1
    sample_group_unique = as.character(df_label$sample_group[!duplicated(df_label$sample_group)])
    sp_unique = as.character(df_label$scientific_name[!duplicated(df_label$scientific_name)])
    sample_group_color_unique = as.character(df_label$sample_group_color[!duplicated(df_label$sample_group_color)])
    species_info = get_species_encoding_info(sp_unique)
    species_shape_map = get_species_shape_map(sp_unique)
    sp_shape_unique = get_point_shapes_for_species(sp_unique, species_shape_map = species_shape_map)
    toumei = rgb(1, 1, 1, 0)
    if (species_info[['use_shape']]) {
        legend_text = c('Sample group', as.character(sample_group_unique), "", 'Species', as.character(sp_unique))
        legend_pch = c(16, rep(16, length(sample_group_color_unique)), 16, 16, sp_shape_unique)
        legend_col = c(toumei, sample_group_color_unique, toumei, toumei, rep('black', length(sp_unique)))
        legend_font = c(2, rep(1, length(sample_group_color_unique)), 1, 2, rep(3, length(sp_unique)))
    } else {
        legend_text = c('Sample group', as.character(sample_group_unique), "", 'Species', 'Centroid labels on plots')
        legend_pch = c(16, rep(16, length(sample_group_color_unique)), 16, 16, 16)
        legend_col = c(toumei, sample_group_color_unique, toumei, toumei, 'black')
        legend_font = c(2, rep(1, length(sample_group_color_unique)), 1, 2, 1)
    }
    plot.new()
    legend("right", legend = legend_text, pt.cex = 1, pch = legend_pch, lty = 0, lwd = 2, col = legend_col, cex = cex_axis, text.font = legend_font)
}

prepare_metadata_table = function(dir_csca_input_table, selected_sample_groups, spp) {
    files = list.files(dir_csca_input_table, pattern = ".*metadata.*")
    metadata_paths = file.path(dir_csca_input_table, files)
    metadata_tables = lapply(metadata_paths, function(metadata_path) {
        read.table(metadata_path, header = TRUE, sep = '\t', quote = '', comment.char = '', check.names = FALSE)
    })
    if (length(metadata_tables) == 0) {
        stop('No metadata files found in csca input table directory.')
    }
    df_metadata = do.call(rbind, metadata_tables)
    df_metadata = df_metadata[(df_metadata[['sample_group']] %in% selected_sample_groups) & (df_metadata[['scientific_name']] %in% spp),]
    df_metadata = df_metadata[, !startsWith(colnames(df_metadata), 'Unnamed')]
    return(df_metadata)
}

get_label_orders = function(df_metadata) {
    order_cg = order(df_metadata[['sample_group']])
    label_orders = unique(paste(df_metadata[order_cg, 'scientific_name'], df_metadata[order_cg, 'sample_group'], sep = '_'))
    label_orders = sub(' ', '_', label_orders)
    return(label_orders)
}

extract_ortholog_expression_table_core = function(df_singleog, expression_tables,
                                                  transform_species_tc = NULL,
                                                  transform_combined_tc = NULL) {
    orthologs = list()
    rowname_order = rownames(df_singleog)
    for (d in correction_labels) {
        tc_slices = list()
        for (sp_filled in colnames(df_singleog)) {
            tc = expression_tables[[d]][[sp_filled]]
            if (is.null(tc)) {
                next
            }
            if (!is.null(transform_species_tc)) {
                tc = transform_species_tc(tc = tc, sp_filled = sp_filled, correction = d)
            }
            row_idx = match(df_singleog[[sp_filled]], rownames(tc))
            tc_slice = tc[row_idx, , drop = FALSE]
            rownames(tc_slice) = rowname_order
            tc_slices[[length(tc_slices) + 1]] = tc_slice
        }
        if (length(tc_slices) == 0) {
            combined_tc = data.frame(matrix(nrow = length(rowname_order), ncol = 0))
            rownames(combined_tc) = rowname_order
        } else {
            combined_tc = as.data.frame(do.call(cbind, tc_slices), check.names = FALSE)
            rownames(combined_tc) = rowname_order
        }
        if (!is.null(transform_combined_tc)) {
            combined_tc = transform_combined_tc(tc = combined_tc, correction = d)
        }
        orthologs[[d]] = combined_tc
    }
    return(orthologs)
}

extract_ortholog_mean_expression_table = function(df_singleog, averaged_tcs, label_orders) {
    averaged_orthologs = extract_ortholog_expression_table_core(
        df_singleog = df_singleog,
        expression_tables = averaged_tcs,
        transform_combined_tc = function(tc, correction) {
            tc = sort_averaged_tc(tc)
            available_label_orders = label_orders[label_orders %in% colnames(tc)]
            tc = tc[, available_label_orders, drop = FALSE]
            tc
        }
    )
    cat(nrow(averaged_orthologs[['corrected']]), 'orthologs were found before filtering.\n')
    return(averaged_orthologs)
}

load_unaveraged_expression_tables = function(dir_csca_input_table, spp_filled, batch_effect_alg) {
    unaveraged_tcs = list()
    unaveraged_tcs[['uncorrected']] = list()
    unaveraged_tcs[['corrected']] = list()
    all_files = list.files(dir_csca_input_table, pattern = "*.tc.tsv")
    uncorrected_files = all_files[grepl("uncorrected", all_files)]
    corrected_files = all_files[((!grepl("uncorrected", all_files)) & (grepl(batch_effect_alg, all_files)))]
    for (sp in spp_filled) {
        uncorrected_file = uncorrected_files[startsWith(uncorrected_files, sp)]
        uncorrected_path = file.path(dir_csca_input_table, uncorrected_file)
        corrected_file = corrected_files[startsWith(corrected_files, sp)]
        corrected_path = file.path(dir_csca_input_table, corrected_file)
        if ((length(uncorrected_path) == 0) | (length(corrected_path) == 0)) {
            cat(paste("Skipping. `amalgkit curate` output(s) not found:", sp, "\n"), file = stderr())
            next
        }
        unaveraged_tcs[['uncorrected']][[sp]] = read.delim(uncorrected_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
        unaveraged_tcs[['corrected']][[sp]] = read.delim(corrected_path, header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
    }
    return(unaveraged_tcs)
}

extract_ortholog_unaveraged_expression_table = function(df_singleog, unaveraged_tcs) {
    unaveraged_orthologs = extract_ortholog_expression_table_core(
        df_singleog = df_singleog,
        expression_tables = unaveraged_tcs,
        transform_species_tc = function(tc, sp_filled, correction) {
            colnames(tc) = paste(sp_filled, colnames(tc), sep = '_')
            tc
        },
        transform_combined_tc = function(tc, correction) {
            tc_order = order(sub('.*_', '', colnames(tc)))
            tc = tc[, tc_order, drop = FALSE]
            tc
        }
    )
    return(unaveraged_orthologs)
}

get_df_labels_averaged = function(df_metadata, label_orders, selected_sample_groups, sample_group_colors) {
    metadata_tmp = df_metadata[(df_metadata[['exclusion']] == 'no'),]
    df_label = unique(metadata_tmp[, c('scientific_name', 'sample_group')])
    metadata_key = paste(metadata_tmp[['scientific_name']], metadata_tmp[['sample_group']], sep = '|||')
    label_key = paste(df_label[['scientific_name']], df_label[['sample_group']], sep = '|||')
    num_bp_by_key = tapply(metadata_tmp[['bioproject']], metadata_key, function(x) { length(unique(x)) })
    num_run_by_key = tapply(metadata_tmp[['run']], metadata_key, function(x) { length(unique(x)) })
    df_label[['num_bp']] = as.numeric(num_bp_by_key[label_key])
    df_label[['num_run']] = as.numeric(num_run_by_key[label_key])
    df_label = df_label[order(df_label[['sample_group']], df_label[['scientific_name']]),]
    df_label = sort_labels(df_label, label_orders)
    df_label = add_color_to_metadata(df_label, selected_sample_groups, sample_group_colors)
    rownames(df_label) = NULL
    write.table(df_label, paste0('csca_color_averaged.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
    return(df_label)
}

get_df_labels_unaveraged = function(df_metadata, selected_sample_groups, sample_group_colors) {
    cols = c('run', 'bioproject', 'sample_group', 'scientific_name', 'sp_color', 'sample_group_color', 'bp_color')
    metadata_tmp = df_metadata[(df_metadata[['exclusion']] == 'no'),]
    df_color = add_color_to_metadata(metadata_tmp, selected_sample_groups, sample_group_colors)
    df_color = df_color[, cols]
    label_order = order(df_color[['run']])
    df_color = df_color[label_order,]
    write.table(df_color, paste0('csca_color_unaveraged.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
    return(df_color)
}

add_color_to_metadata = function(df, selected_sample_groups, sample_group_colors) {
    df = df[, (!colnames(df) %in% c('bp_color', 'sp_color', 'sample_group_color'))]
    scientific_name = as.character(df[['scientific_name']])
    sample_group = as.character(df[['sample_group']])
    scientific_name_unique = sort(scientific_name[!duplicated(scientific_name)])
    if (length(sample_group_colors) == 1 && sample_group_colors == 'DEFAULT') {
        if (length(selected_sample_groups) <= 8) {
            sample_group_color = brewer.pal(length(unique(sample_group)), "Dark2")
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
        } else if (length(selected_sample_groups) <= 12) {
            sample_group_color = brewer.pal(length(unique(sample_group)), "Paired")
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
        } else {
            sample_group_color = rainbow_hcl(length(selected_sample_groups), c = 100)
            sp_color = rainbow_hcl(length(unique(scientific_name)), c = 150)
        }
    } else {
        if (length(sample_group_colors) != length(selected_sample_groups)) {
            stop("Length of sample_group_colors must match length of selected_sample_groups")
        }
        sample_group_color = sample_group_colors
        sp_color = rainbow_hcl(length(unique(scientific_name)), c = 100)
    }
    sample_group_unique = selected_sample_groups
    sample_group_palette = sample_group_color[1:length(sample_group_unique)]
    names(sample_group_palette) = sample_group_unique
    sp_palette = sp_color[1:length(scientific_name_unique)]
    names(sp_palette) = scientific_name_unique
    is_known_sample_group = sample_group %in% sample_group_unique
    df = df[is_known_sample_group, , drop = FALSE]
    scientific_name = as.character(df[['scientific_name']])
    sample_group = as.character(df[['sample_group']])
    df[['sp_color']] = unname(sp_palette[scientific_name])
    df[['sample_group_color']] = unname(sample_group_palette[sample_group])
    if ('bioproject' %in% colnames(df)) {
        bioproject = as.character(df[['bioproject']])
        bp_unique = unique(bioproject)
        bp_color = rainbow_hcl(length(bp_unique), c = 50)
        bp_palette = bp_color[1:length(bp_unique)]
        names(bp_palette) = bp_unique
        df[['bp_color']] = unname(bp_palette[bioproject])
    }
    return(df)
}

save_averaged_tsne_plot = function(tc, df_label) {
    cat('Generating averaged t-SNE plot.\n')
    tsne_result = compute_tsne_embedding(tc = tc)
    if (!tsne_result[['ok']]) {
        if (tsne_result[['reason']] == 'insufficient_samples') {
            cat('Skipping averaged t-SNE: insufficient number of samples.\n')
        } else {
            cat("t-SNE failed and will be skipped.\n")
        }
        return()
    }
    try_out = tryCatch(
    {
        draw_species_scatter_base(
            x = tsne_result[['Y']][, 1],
            y = tsne_result[['Y']][, 2],
            df_label = df_label,
            xlab = "t-SNE dimension 1",
            ylab = "t-SNE dimension 2",
            cex = 2,
            lwd = 1,
            las = 1
        )
    },
        error = function(a) { return("t-SNE plot failed.") }
    )
    if (mode(try_out) == "character") {
        cat('t-SNE failed.\n')
        draw_empty_panel()
    }
}

get_pca_coordinates = function(tc, df_label, by = 'species_sample_group') {
    tc_dist_matrix = cor(tc, method = 'pearson')
    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
    #set.seed(1)
    pca = prcomp(tc_dist_matrix)
    pca_summary = summary(pca)
    num_pc_available = ncol(pca[['x']])
    labels = c()
    for (i in 1:5) {
        if (i <= num_pc_available) {
            labels = c(labels, paste0("Principal component ", i, " (", round(pca_summary$importance[2, i] * 100, digits = 1), "%)"))
        } else {
            labels = c(labels, paste0("Principal component ", i, " (NA%)"))
        }
    }
    pca_x = pca[['x']]
    empty_pc = rep(NA_real_, nrow(pca_x))
    PC1 = if ('PC1' %in% colnames(pca_x)) pca_x[, 'PC1'] else empty_pc
    PC2 = if ('PC2' %in% colnames(pca_x)) pca_x[, 'PC2'] else empty_pc
    PC3 = if ('PC3' %in% colnames(pca_x)) pca_x[, 'PC3'] else empty_pc
    PC4 = if ('PC4' %in% colnames(pca_x)) pca_x[, 'PC4'] else empty_pc
    PC5 = if ('PC5' %in% colnames(pca_x)) pca_x[, 'PC5'] else empty_pc
    tmp = data.frame(PC1, PC2, PC3, PC4, PC5)
    if (by == 'species_sample_group') {
        df_label[by] = paste0(sub(' ', '_', df_label[['scientific_name']]), '_', df_label[['sample_group']])
        tmp[by] = rownames(tmp)
    } else if (by == 'run') {
        tmp[by] = sub('.*_', '', rownames(tmp))
    } else {
        tmp[by] = rownames(tmp)
    }
    tmp = merge(df_label, tmp, by = by)
    return(list(tmp, labels))
}

save_unaveraged_pca_plot = function(unaveraged_orthologs, df_color_unaveraged, df_metadata) {
    cat('Generating unaveraged PCA plot.\n')
    sample_group_colors = get_sample_group_palette(df_color_unaveraged)
    for (d in correction_labels) {
        out = get_pca_coordinates(tc = unaveraged_orthologs[[d]], df_label = df_color_unaveraged, by = 'run')
        tmp = out[[1]]
        pc_contributions = out[[2]]
        pc_cols = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')
        pc_cols2 = paste(pc_cols, d, sep = '_')
        sorted_cols = c(colnames(df_metadata), pc_cols2)
        tmp2 = tmp[, c('run', pc_cols)]
        colnames(tmp2) = c('run', pc_cols2)
        df_metadata = append_columns_by_key(left_df = df_metadata, right_df = tmp2, key_col = 'run')
        df_metadata = df_metadata[, sorted_cols]

        for (pcxy in list(c(1, 2), c(3, 4))) {
            pcx = pcxy[1]
            pcy = pcxy[2]
            colx = paste0('PC', pcx)
            coly = paste0('PC', pcy)

            axis_limits = compute_scatter_axis_limits(
                xvals = tmp[[colx]],
                yvals = tmp[[coly]],
                require_variation = TRUE,
                pad_ratio = 0.01,
                flat_half_span = 0.5
            )
            if (!axis_limits[['ok']]) {
                if (axis_limits[['reason']] == 'no_finite_data') {
                    message(sprintf("Skipping PC%d vs PC%d plot for '%s' - no finite data.", pcx, pcy, d))
                } else if (axis_limits[['reason']] == 'no_variation') {
                    message(sprintf("Skipping PC%d vs PC%d plot for '%s' - no variation in data.", pcx, pcy, d))
                }
                next
            }

            g = build_unaveraged_scatter_plot(
                df = tmp,
                x_col = colx,
                y_col = coly,
                x_label = pc_contributions[pcx],
                y_label = pc_contributions[pcy],
                axis_limits = axis_limits,
                sample_group_colors = sample_group_colors,
                font_size = font_size,
                point_size = 0.7,
                point_alpha = 0.5,
                na_rm = TRUE
            )
            filename = paste0('csca_unaveraged_pca_PC', pcx, pcy, '.', d, '.pdf')
            save_scatter_plot_pdf(
                g = g,
                filename = filename,
                fail_prefix = "PCA could not be computed for file",
                width = 4.25,
                height = 2.15
            )
        }
    }
    return(df_metadata)
}

get_tsne_coordinates = function(tc, df_label, by = 'run') {
    tsne_result = compute_tsne_embedding(tc = tc)
    if (!tsne_result[['ok']]) {
        return(NULL)
    }
    tmp = data.frame(tsne1 = tsne_result[['Y']][, 1], tsne2 = tsne_result[['Y']][, 2])
    tmp[[by]] = sub('.*_', '', colnames(tc))
    tmp = merge(df_label, tmp, by = by)
    return(tmp)
}

save_unaveraged_tsne_plot = function(unaveraged_orthologs, df_color_unaveraged) {
    cat('Generating unaveraged t-SNE plot.\n')
    sample_group_colors = get_sample_group_palette(df_color_unaveraged)
    for (d in correction_labels) {
        tmp = get_tsne_coordinates(tc = unaveraged_orthologs[[d]], df_label = df_color_unaveraged)
        if (is.null(tmp)) {
            cat(sprintf('Skipping unaveraged t-SNE (%s): insufficient number of samples.\n', d))
            next
        }
        pcx = 1
        pcy = 2
        colx = paste0('tsne', pcx)
        coly = paste0('tsne', pcy)
        axis_limits = compute_scatter_axis_limits(
            xvals = tmp[[colx]],
            yvals = tmp[[coly]],
            require_variation = FALSE,
            pad_ratio = 0.01,
            flat_half_span = 0.5
        )
        if (!axis_limits[['ok']]) {
            cat(sprintf('Skipping unaveraged t-SNE (%s): no finite coordinates.\n', d))
            next
        }

        g = build_unaveraged_scatter_plot(
            df = tmp,
            x_col = colx,
            y_col = coly,
            x_label = 't-SNE dimension 1',
            y_label = 't-SNE dimension 2',
            axis_limits = axis_limits,
            sample_group_colors = sample_group_colors,
            font_size = font_size,
            point_size = 0.7,
            point_alpha = 1,
            na_rm = TRUE
        )
        filename = paste0('csca_unaveraged_tsne.', d, '.pdf')
        save_scatter_plot_pdf(
            g = g,
            filename = filename,
            fail_prefix = "t-SNE could not be computed for file",
            width = 4.25,
            height = 2.15
        )
    }
}

save_averaged_heatmap_plot = function(averaged_orthologs, df_color_averaged, averaged_plot_cache = NULL) {
    cat('Generating averaged heatmap.\n')
    file_name = 'csca_SVA_heatmap.pdf'
    with_pdf_defaults(
        file = file_name,
        width = 7.2,
        height = 3.3,
        font_size = font_size,
        font_family = font_family,
        plot_fn = function(local_font_size, local_font_family) {
            layout_matrix = matrix(c(
                1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                1, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                2, 10, byrow = TRUE
            )
            layout(t(layout_matrix))
            par(mar = c(0, 0, 0, 0))
            draw_empty_panel()
            text(0.27, 0.5, 'Uncorrected', srt = 0, cex = 1)
            text(0.80, 0.5, 'Corrected', srt = 0, cex = 1)
            for_each_averaged_correction(averaged_orthologs, df_color_averaged, function(correction, tc, df_label) {
                tc_dist_matrix = get_cached_plot_data(averaged_plot_cache, correction, 'tc_cor_plot', fallback = NULL)
                par(mar = c(0, 0, 0, 0))
                draw_multisp_heatmap(
                    tc = tc,
                    df_label = df_label,
                    tc_dist_matrix = tc_dist_matrix,
                    font_size = local_font_size,
                    show_colorbar = (correction == 'corrected')
                )
            })
        }
    )
}

save_averaged_dendrogram_plot = function(averaged_orthologs, df_color_averaged, averaged_plot_cache = NULL) {
    cat('Generating averaged dendrogram.\n')
    file_name = 'csca_SVA_dendrogram.pdf'
    with_pdf_defaults(
        file = file_name,
        width = 7.2,
        height = 2.5,
        font_size = font_size,
        font_family = font_family,
        plot_fn = function(local_font_size, local_font_family) {
            layout_matrix = matrix(c(1, 2), 2, 1, byrow = TRUE)
            layout(t(layout_matrix))
            for_each_averaged_correction(averaged_orthologs, df_color_averaged, function(correction, tc, df_label) {
                tc_dist_dist = get_cached_plot_data(averaged_plot_cache, correction, 'tc_dist_dendrogram', fallback = NULL)
                label_cex = compute_dense_label_cex(
                    num_labels = ncol(tc),
                    base_pt = local_font_size,
                    min_pt = 3.5,
                    soft_limit = 18
                )
                par(mar = c(10, 5.5, 0, 0), mgp = c(4, 0.7, 0))
                draw_multisp_dendrogram(
                    tc = tc,
                    df_label = df_label,
                    df_metadata = df_metadata,
                    label_cex = label_cex,
                    axis_cex = 1,
                    tc_dist_dist = tc_dist_dist
                )
            })
        }
    )
}

save_averaged_dimensionality_reduction_summary = function(averaged_orthologs, df_color_averaged, averaged_plot_cache = NULL) {
    cat('Generating averaged dimensionality reduction summary.\n')
    file_name = 'csca_averaged_summary.pdf'
    with_pdf_defaults(
        file = file_name,
        width = 7.2,
        height = 7.2,
        font_size = font_size,
        font_family = font_family,
        plot_fn = function(local_font_size, local_font_family) {
            layout_matrix = matrix(c(1, 1, 1, 4, 4, 4, 7, 7, 2, 2, 2, 5, 5, 5, 7, 7, 3, 3, 3, 6, 6, 6, 7, 7), 3, 8, byrow = TRUE)
            layout(layout_matrix)
            for_each_averaged_correction(averaged_orthologs, df_color_averaged, function(correction, tc, df_label) {
                tc_dist_matrix = get_cached_plot_data(averaged_plot_cache, correction, 'tc_cor_plot', fallback = NULL)
                tc_dist_dist = get_cached_plot_data(averaged_plot_cache, correction, 'tc_dist_mds', fallback = NULL)
                par(mar = c(4, 4, 0.1, 1)); draw_multisp_pca(tc = tc, df_label = df_label, tc_dist_matrix = tc_dist_matrix)
                par(mar = c(4, 4, 0.1, 1)); draw_multisp_tsne(tc = tc, df_label = df_label)
                par(mar = c(4, 4, 0.1, 1)); draw_multisp_mds(tc = tc, df_label = df_label, tc_dist_dist = tc_dist_dist)
            })
            par(mar = c(0, 0, 0, 0)); draw_multisp_legend(df_label)
        }
    )
}

draw_multisp_boxplot = function(df_metadata, tc_dist_matrix, fontsize = 8) {
    is_same_sp = outer(df_metadata[['scientific_name']], df_metadata[['scientific_name']], function(x, y) { x == y })
    is_same_sample_group = outer(df_metadata[['sample_group']], df_metadata[['sample_group']], function(x, y) { x == y })
    plot(c(0.5, 4.5), c(0, 1), type = 'n', xlab = '', ylab = "Pearson's correlation\ncoefficient", las = 1, xaxt = 'n')
    safe_boxplot = function(values, at) {
        values = values[is.finite(values)]
        if (length(values) == 0) {
            return(invisible(FALSE))
        }
        boxplot(values, at = at, add = TRUE, col = 'gray', yaxt = 'n')
        return(invisible(TRUE))
    }
    safe_boxplot(tc_dist_matrix[(!is_same_sp) & (!is_same_sample_group)], at = 1)
    safe_boxplot(tc_dist_matrix[(is_same_sp) & (!is_same_sample_group)], at = 2)
    safe_boxplot(tc_dist_matrix[(!is_same_sp) & (is_same_sample_group)], at = 3)
    safe_boxplot(tc_dist_matrix[(is_same_sp) & (is_same_sample_group)], at = 4)
    labels = c('bw\nbw', 'bw\nwi', 'wi\nbw', 'wi\nwi')
    axis(side = 1, at = c(1, 2, 3, 4), labels = labels, padj = 0.5)
    axis(side = 1, at = 0.35, labels = 'Sample group\nSpecies', padj = 0.5, hadj = 1, tick = FALSE)
}

save_averaged_box_plot = function(averaged_orthologs, df_color_averaged, averaged_plot_cache = NULL) {
    cat('Generating averaged boxplot.\n')
    file_name = 'csca_boxplot.pdf'
    with_pdf_defaults(
        file = file_name,
        width = 7.2,
        height = 3.6,
        font_size = font_size,
        font_family = font_family,
        plot_fn = function(local_font_size, local_font_family) {
            par(mfrow = c(1, 2))
            for_each_averaged_correction(averaged_orthologs, df_color_averaged, function(correction, tc, df_label) {
                tc[tc < 0] = 0
                tc_dist_matrix = get_cached_plot_data(averaged_plot_cache, correction, 'tc_cor_boxplot', fallback = NULL)
                if (is.null(tc_dist_matrix)) {
                    tc_dist_matrix = as.matrix(suppressWarnings(cor(tc, method = 'pearson')))
                    tc_dist_matrix[is.na(tc_dist_matrix)] = 0
                }
                draw_multisp_boxplot(df_label, tc_dist_matrix, fontsize = local_font_size)
            })
        }
    )
}

calculate_correlation_within_group = function(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_sample_groups, dist_method = 'pearson') {
    median_fun = function(x) { median(x, na.rm = TRUE) }
    metadata_keys = paste(df_metadata[['scientific_name']], df_metadata[['run']], sep = '|||')
    metadata_rows_by_key = split(seq_len(nrow(df_metadata)), metadata_keys)
    metadata_first_row = vapply(metadata_rows_by_key, function(idx) idx[[1]], integer(1))
    metadata_sample_group_by_key = as.character(df_metadata[['sample_group']][metadata_first_row])
    sample_index_cache = list()

    build_sample_index = function(sample_ids) {
        sample_species = sub('_', ' ', sub('^([^_]+_[^_]+)_.*$', '\\1', sample_ids))
        sample_runs = sub('^[^_]+_[^_]+_', '', sample_ids)
        sample_keys = paste(sample_species, sample_runs, sep = '|||')
        rows_by_sample = unname(metadata_rows_by_key[sample_keys])
        has_metadata = !vapply(rows_by_sample, is.null, logical(1))
        sample_groups = as.character(metadata_sample_group_by_key[sample_keys])
        list(
            rows_by_sample = rows_by_sample,
            has_metadata = has_metadata,
            sample_groups = sample_groups
        )
    }

    for (d in correction_labels) {
        averaged_tc = averaged_orthologs[[d]]
        unaveraged_tc = unaveraged_orthologs[[d]]
        ortholog_med = data.frame(matrix(NA_real_, nrow(averaged_tc), length(selected_sample_groups)))
        colnames(ortholog_med) = selected_sample_groups
        rownames(ortholog_med) = rownames(averaged_tc)

        group_col_idx = lapply(selected_sample_groups, function(sample_group) {
            which(endsWith(colnames(averaged_tc), sample_group))
        })
        names(group_col_idx) = selected_sample_groups
        for (sample_group in selected_sample_groups) {
            idx = group_col_idx[[sample_group]]
            if (length(idx) == 0) {
                ortholog_med[, sample_group] = NA_real_
                next
            }
            ortholog_med[, sample_group] = apply(as.data.frame(averaged_tc[, idx, drop = FALSE]), 1, median_fun)
        }

        stopifnot(all(rownames(unaveraged_tc) == rownames(ortholog_med)))
        target_col = paste0('within_group_cor_', d)
        nongroup_col = paste0('max_nongroup_cor_', d)
        df_metadata[, target_col] = NA_real_
        df_metadata[, nongroup_col] = NA_real_

        cor_by_group = as.matrix(suppressWarnings(cor(unaveraged_tc, ortholog_med, method = dist_method, use = 'pairwise.complete.obs')))
        sample_ids = colnames(unaveraged_tc)
        group_ids = colnames(ortholog_med)

        is_expected_orientation = identical(rownames(cor_by_group), sample_ids) && identical(colnames(cor_by_group), group_ids)
        is_transposed_orientation = identical(rownames(cor_by_group), group_ids) && identical(colnames(cor_by_group), sample_ids)
        if (is_transposed_orientation) {
            cor_by_group = t(cor_by_group)
        }
        expected_nrow = length(sample_ids)
        expected_ncol = length(group_ids)
        if ((!is_expected_orientation) && (!is_transposed_orientation)) {
            if ((nrow(cor_by_group) == expected_ncol) && (ncol(cor_by_group) == expected_nrow)) {
                cor_by_group = t(cor_by_group)
            } else if ((nrow(cor_by_group) != expected_nrow) || (ncol(cor_by_group) != expected_ncol)) {
                stop('Unexpected correlation matrix shape in calculate_correlation_within_group().')
            }
        }
        rownames(cor_by_group) = sample_ids
        colnames(cor_by_group) = group_ids

        sample_cache_key = paste(sample_ids, collapse = '\r')
        if (!sample_cache_key %in% names(sample_index_cache)) {
            sample_index_cache[[sample_cache_key]] = build_sample_index(sample_ids)
        }
        sample_index = sample_index_cache[[sample_cache_key]]

        missing_metadata = !sample_index[['has_metadata']]
        if (any(missing_metadata)) {
            warning(sprintf('Samples skipped (missing metadata): %s', paste(sample_ids[missing_metadata], collapse = ', ')))
        }
        group_idx = match(sample_index[['sample_groups']], group_ids)
        missing_group = sample_index[['has_metadata']] & is.na(group_idx)
        if (any(missing_group)) {
            warning(sprintf(
                'Sample groups not found in ortholog medians: %s',
                paste(unique(sample_index[['sample_groups']][missing_group]), collapse = ', ')
            ))
        }

        valid_sample_idx = which(sample_index[['has_metadata']] & !is.na(group_idx))
        if (length(valid_sample_idx) == 0) {
            next
        }
        within_cor = rep(NA_real_, length(sample_ids))
        within_cor[valid_sample_idx] = cor_by_group[cbind(valid_sample_idx, group_idx[valid_sample_idx])]
        max_nongroup_cor = rep(NA_real_, length(sample_ids))
        cor_valid = cor_by_group[valid_sample_idx, , drop = FALSE]
        if (ncol(cor_valid) > 1) {
            valid_group_idx = group_idx[valid_sample_idx]
            cor_valid[cbind(seq_len(nrow(cor_valid)), valid_group_idx)] = NA_real_
            max_nongroup_valid = apply(cor_valid, 1, function(row_vals) {
                nongroup_vals = row_vals[is.finite(row_vals)]
                if (length(nongroup_vals) == 0) {
                    return(NA_real_)
                }
                max(nongroup_vals)
            })
            max_nongroup_cor[valid_sample_idx] = as.numeric(max_nongroup_valid)
        }

        metadata_rows_for_samples = sample_index[['rows_by_sample']][valid_sample_idx]
        repeat_counts = lengths(metadata_rows_for_samples)
        metadata_rows_flat = unlist(metadata_rows_for_samples, use.names = FALSE)
        df_metadata[metadata_rows_flat, target_col] = rep(within_cor[valid_sample_idx], repeat_counts)
        df_metadata[metadata_rows_flat, nongroup_col] = rep(max_nongroup_cor[valid_sample_idx], repeat_counts)
    }
    return(df_metadata)
}

save_group_cor_histogram = function(df_metadata, df_color_unaveraged, font_size = 8) {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family('Helvetica')
    cat('Generating unaveraged group correlation histogram.\n')

    cor_cols = c('within_group_cor_uncorrected', 'within_group_cor_corrected')
    fill_by_vars = c('sample_group', 'scientific_name')
    valid_rows = Reduce('&', lapply(cor_cols, function(col) {
        is.finite(df_metadata[[col]]) & !is.na(df_metadata[[col]]) & (df_metadata[[col]] >= 0) & (df_metadata[[col]] <= 1)
    }))
    df_clean = df_metadata[valid_rows, , drop = FALSE]

    build_palette = function(df, label_col, color_col) {
        tmp = unique(df[, c(label_col, color_col)])
        tmp = tmp[!is.na(tmp[[label_col]]) & !is.na(tmp[[color_col]]), , drop = FALSE]
        out = as.character(tmp[[color_col]])
        names(out) = as.character(tmp[[label_col]])
        out
    }
    species_palette = build_palette(df_color_unaveraged, 'scientific_name', 'sp_color')
    sample_group_palette = build_palette(df_color_unaveraged, 'sample_group', 'sample_group_color')

    make_legend_panel = function(palette, legend_title) {
        if (length(palette) == 0) {
            return(ggplot2::ggplot() + theme_void(base_size = font_size, base_family = font_family))
        }
        legend_df = data.frame(
            x = rep(1, length(palette)),
            y = rep(1, length(palette)),
            group = factor(names(palette), levels = names(palette)),
            stringsAsFactors = FALSE
        )
        ggplot2::ggplot(legend_df, ggplot2::aes(x = x, y = y, color = group)) +
            ggplot2::geom_point(alpha = 0, show.legend = TRUE) +
            ggplot2::scale_color_manual(values = palette, breaks = names(palette), drop = FALSE, name = legend_title) +
            ggplot2::theme_void(base_size = font_size, base_family = font_family) +
            ggplot2::guides(color = ggplot2::guide_legend(
                ncol = 1,
                byrow = TRUE,
                order = 1,
                override.aes = list(shape = 15, size = 4, alpha = 1)
            )) +
            ggplot2::theme(
                legend.position = 'center',
                legend.direction = 'vertical',
                legend.box = 'vertical',
                legend.title = ggplot2::element_text(size = font_size, family = font_family, color = 'black'),
                legend.text = ggplot2::element_text(size = font_size, family = font_family, color = 'black'),
                legend.key.height = grid::unit(0.12, 'in'),
                legend.key.width = grid::unit(0.18, 'in'),
                plot.margin = ggplot2::margin(0, 0, 0, 0, unit = 'pt')
            )
    }

    plot_list <- list()
    for (col in cor_cols) {
        for (fill_by in fill_by_vars) {
            tmp = df_clean[!is.na(df_clean[[col]]), ]
            if (fill_by == 'scientific_name') {
                palette = species_palette
                legend_title = 'Species'
            } else {
                palette = sample_group_palette
                legend_title = 'Sample group'
            }
            if (nrow(tmp) > 0) {
                tmp[[fill_by]] = factor(as.character(tmp[[fill_by]]), levels = names(palette))
            }
            if (nrow(tmp) == 0) {
                g = ggplot2::ggplot() +
                    theme_void(base_size = font_size, base_family = font_family) +
                    annotate('text', x = 0.5, y = 0.5, label = paste('No finite data:', col), family = font_family, size = font_size / ggplot2::.pt)
            } else {
                g = ggplot2::ggplot(tmp) +
                    geom_histogram(aes(x = !!rlang::sym(col), fill = !!rlang::sym(fill_by)),
                                   position = "stack", alpha = 0.7, bins = 40, na.rm = TRUE, show.legend = FALSE) +
                    scale_fill_manual(values = palette, breaks = names(palette), drop = FALSE, name = legend_title) +
                    theme_bw(base_size = font_size, base_family = font_family) +
                    coord_cartesian(xlim = c(0, 1)) +
                    labs(x = col, y = 'Sample count') +
                    guides(fill = 'none') +
                    build_standard_ggplot_theme(font_size = font_size, font_family = font_family)
            }
            plot_list[[paste0(col, "_", fill_by)]] <- g
        }
    }

    plot_order = c(
        'within_group_cor_uncorrected_scientific_name',
        'within_group_cor_corrected_scientific_name',
        'within_group_cor_uncorrected_sample_group',
        'within_group_cor_corrected_sample_group'
    )
    plots = lapply(plot_order, function(name) {
        if (name %in% names(plot_list)) {
            return(plot_list[[name]])
        }
        ggplot2::ggplot() + theme_void(base_size = font_size, base_family = font_family) + annotate('text', x = 0.5, y = 0.5, label = paste('Missing panel:', name), family = font_family, size = font_size / ggplot2::.pt)
    })

    legend_species_plot = make_legend_panel(species_palette, 'Species')
    legend_sample_group_plot = make_legend_panel(sample_group_palette, 'Sample group')
    max_legend_items = max(length(species_palette), length(sample_group_palette))
    legend_height_in = max(1.0, 0.13 * max_legend_items + 0.35)
    pdf_height = 5.2 + legend_height_in

    with_pdf_defaults(
        file = "csca_within_group_cor.pdf",
        width = 7.2,
        height = pdf_height,
        font_size = font_size,
        font_family = font_family,
        plot_fn = function(local_font_size, local_font_family) {
            grid::grid.newpage()
            grid::pushViewport(
                grid::viewport(
                    layout = grid::grid.layout(
                        nrow = 3,
                        ncol = 2,
                        heights = grid::unit.c(
                            grid::unit(1, "null"),
                            grid::unit(1, "null"),
                            grid::unit(legend_height_in, "in")
                        )
                    )
                )
            )
            print(plots[[1]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
            print(plots[[2]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
            print(plots[[3]], vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
            print(plots[[4]], vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
            print(legend_species_plot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
            print(legend_sample_group_plot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
        }
    )
}

extract_selected_tc_only = function(unaveraged_tcs, df_metadata) {
    selected_runs = unique(df_metadata[(df_metadata[['exclusion']] == 'no'), 'run'])
    for (d in correction_labels) {
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            is_selected = colnames(unaveraged_tcs[[d]][[sci_name]]) %in% selected_runs
            if (sum(is_selected) == 0) {
                warning(paste('No', d, 'samples were selected:', sci_name))
                next
            }
            unaveraged_tcs[[d]][[sci_name]] = unaveraged_tcs[[d]][[sci_name]][, is_selected]
        }
    }
    return(unaveraged_tcs)
}

unaveraged2averaged = function(unaveraged_tcs, df_metadata, selected_sample_groups) {
    metadata_selected = df_metadata[(df_metadata[['exclusion']] == 'no'), c('scientific_name', 'sample_group', 'run')]
    metadata_selected[['sci_name_ub']] = sub(' ', '_', metadata_selected[['scientific_name']])
    metadata_selected[['sci_group_key']] = paste(metadata_selected[['sci_name_ub']], metadata_selected[['sample_group']], sep = '|||')
    runs_by_sci_group = split(metadata_selected[['run']], metadata_selected[['sci_group_key']])

    averaged_tcs = list()
    for (d in correction_labels) {
        averaged_tcs[[d]] = list()
        scientific_names = names(unaveraged_tcs[[d]])
        for (sci_name in scientific_names) {
            tc_df = unaveraged_tcs[[d]][[sci_name]]
            tc_mat = as.matrix(tc_df)
            tc_colnames = colnames(tc_mat)
            averaged_cols = list()
            averaged_colnames = c()
            for (sample_group in selected_sample_groups) {
                key = paste(sci_name, sample_group, sep = '|||')
                candidate_runs = runs_by_sci_group[[key]]
                if (is.null(candidate_runs) || (length(candidate_runs) == 0)) {
                    next
                }
                run_idx = match(candidate_runs, tc_colnames)
                run_idx = run_idx[!is.na(run_idx)]
                if (length(run_idx) == 0) {
                    next
                }
                label = paste(sci_name, sample_group, sep = '_')
                if (length(run_idx) == 1) {
                    averaged_col = tc_mat[, run_idx]
                } else {
                    averaged_col = rowMeans(tc_mat[, run_idx, drop = FALSE], na.rm = FALSE)
                }
                averaged_cols[[length(averaged_cols) + 1]] = averaged_col
                averaged_colnames = c(averaged_colnames, label)
            }
            if (length(averaged_cols) == 0) {
                averaged_tcs[[d]][[sci_name]] = NULL
            } else {
                averaged_df = as.data.frame(do.call(cbind, averaged_cols), check.names = FALSE)
                if (length(averaged_cols) == 1) {
                    colnames(averaged_df) = averaged_colnames[1]
                } else {
                    colnames(averaged_df) = averaged_colnames
                }
                rownames(averaged_df) = rownames(tc_df)
                averaged_tcs[[d]][[sci_name]] = averaged_df
            }
        }
    }
    return(averaged_tcs)
}

save_group_cor_scatter = function(df_metadata, font_size = 8) {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family('Helvetica')
    cat('Generating unaveraged group correlation scatter plot.\n')
    alpha_value = 0.2
    improvement_xymin = 0.5
    improvement_xymax = 2.0

    # Compute new columns
    df_metadata[['corrected_per_uncorrected_group_cor']] = df_metadata[['within_group_cor_corrected']] / df_metadata[['within_group_cor_uncorrected']]
    df_metadata[['corrected_per_uncorrected_max_nongroup_cor']] = df_metadata[['max_nongroup_cor_corrected']] / df_metadata[['max_nongroup_cor_uncorrected']]

    # Limit values to the specified range
    df_metadata[['corrected_per_uncorrected_group_cor']] = pmax(pmin(df_metadata[['corrected_per_uncorrected_group_cor']], improvement_xymax), improvement_xymin)
    df_metadata[['corrected_per_uncorrected_max_nongroup_cor']] = pmax(pmin(df_metadata[['corrected_per_uncorrected_max_nongroup_cor']], improvement_xymax), improvement_xymin)

    # Remove rows with NA, NaN, or Inf
    df_clean <- df_metadata[is.finite(df_metadata$within_group_cor_uncorrected) &
                            is.finite(df_metadata$within_group_cor_corrected) &
                            is.finite(df_metadata$max_nongroup_cor_uncorrected) &
                            is.finite(df_metadata$max_nongroup_cor_corrected) &
                            is.finite(df_metadata$corrected_per_uncorrected_group_cor) &
                            is.finite(df_metadata$corrected_per_uncorrected_max_nongroup_cor), ]
    if (nrow(df_clean) == 0) {
        p = ggplot() + theme_void(base_size = font_size, base_family = font_family) + annotate('text', x = 0.5, y = 0.5, label = 'No finite data for group-correlation scatter', family = font_family, size = font_size / ggplot2::.pt)
        save_ggplot_grid(c(list(p), rep(list(ggplot() + theme_void(base_size = font_size, base_family = font_family)), 5)), filename = "csca_group_cor_scatter.pdf", width = 7.2, height = 4.8, nrow = 2, ncol = 3)
        return(invisible(NULL))
    }

    make_scatter_panel = function(x_col, y_col, x_limits, y_limits) {
        finite_x = df_clean[[x_col]][is.finite(df_clean[[x_col]])]
        finite_y = df_clean[[y_col]][is.finite(df_clean[[y_col]])]
        has_density_signal = (length(unique(finite_x)) > 1) && (length(unique(finite_y)) > 1)
        p = ggplot(df_clean, aes(x = !!rlang::sym(x_col), y = !!rlang::sym(y_col))) +
            coord_cartesian(xlim = x_limits, ylim = y_limits) +
            geom_point(alpha = alpha_value, na.rm = TRUE) +
            geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'blue') +
            theme_bw(base_size = font_size, base_family = font_family) +
            build_standard_ggplot_theme(font_size = font_size, font_family = font_family)
        if (has_density_signal) {
            p = p + stat_density_2d(bins = 12, linewidth = 0.25, color = 'gray', na.rm = TRUE)
        }
        return(p)
    }

    ps = list(
        make_scatter_panel('max_nongroup_cor_uncorrected', 'within_group_cor_uncorrected', c(0, 1), c(0, 1)),
        make_scatter_panel('max_nongroup_cor_corrected', 'within_group_cor_corrected', c(0, 1), c(0, 1)),
        make_scatter_panel('within_group_cor_uncorrected', 'within_group_cor_corrected', c(0, 1), c(0, 1)),
        make_scatter_panel('max_nongroup_cor_uncorrected', 'max_nongroup_cor_corrected', c(0, 1), c(0, 1)),
        make_scatter_panel('corrected_per_uncorrected_max_nongroup_cor', 'corrected_per_uncorrected_group_cor',
                           c(improvement_xymin, improvement_xymax), c(improvement_xymin, improvement_xymax)),
        ggplot() + theme_void(base_size = font_size, base_family = font_family)
    )
    save_ggplot_grid(ps, filename = "csca_group_cor_scatter.pdf", width = 7.2, height = 4.8, nrow = 2, ncol = 3)
}

write_pivot_table = function(df_metadata, unaveraged_tcs, selected_sample_groups) {
    d = 'corrected'
    is_selected_cg = (df_metadata[['sample_group']] %in% selected_sample_groups)
    is_loaded_run = FALSE
    for (sci_name_ub in names(unaveraged_tcs[[d]])) {
        sci_name = sub('_', ' ', sci_name_ub)
        is_loaded_run_sp = ((df_metadata[['run']] %in% colnames(unaveraged_tcs[[d]][[sci_name_ub]])) & (df_metadata[['scientific_name']] == sci_name))
        is_loaded_run = (is_loaded_run | is_loaded_run_sp)
    }
    is_selected = (is_selected_cg & is_loaded_run)
    tmp = df_metadata[is_selected, c('scientific_name', 'sample_group')]
    pivot_table = as.data.frame.matrix(table(tmp[['scientific_name']], tmp[['sample_group']]))
    pivot_table = cbind(data.frame(scientific_name = rownames(pivot_table)), pivot_table)
    write.table(pivot_table, 'csca_pivot_selected_samples.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
}

save_delta_pcc_plot = function(directory, plot_title) {
    local_font_size = resolve_csca_font_size(font_size)
    local_font_family = resolve_csca_font_family('Helvetica')

    first_row_means <- list()
    last_row_means <- list()
    all_data_points <- list()
    mean_cols = c('bwbw_mean', 'bwwi_mean', 'wiwi_mean', 'wibw_mean')

    file_list <- list.files(directory, pattern = "\\.correlation_statistics\\.tsv$")
    if (length(file_list) == 0) {
        cat('No correlation statistics files found. Skipping delta PCC plot.\n')
        return(invisible(NULL))
    }
    for (filename in file_list) {
        filepath <- file.path(directory, filename)
        df <- read.table(filepath, sep = '\t', header = TRUE)
        if ((nrow(df) == 0) || (!all(mean_cols %in% colnames(df)))) {
            next
        }
        first_row_means[[filename]] <- as.numeric(df[1, mean_cols, drop = TRUE])
        last_row_means[[filename]] <- as.numeric(df[nrow(df), mean_cols, drop = TRUE])
        all_data_points[[filename]] <- as.numeric(unlist(df[mean_cols]))
    }
    if (length(first_row_means) == 0) {
        cat('No valid correlation statistics rows were found. Skipping delta PCC plot.\n')
        return(invisible(NULL))
    }

    first_row_means_df <- as.data.frame(do.call(rbind, first_row_means), stringsAsFactors = FALSE)
    colnames(first_row_means_df) <- c('bwbw_uncorrected', 'bwwi_uncorrected', 'wiwi_uncorrected', 'wibw_uncorrected')
    last_row_means_df <- as.data.frame(do.call(rbind, last_row_means), stringsAsFactors = FALSE)
    colnames(last_row_means_df) <- c('bwbw_corrected', 'bwwi_corrected', 'wiwi_corrected', 'wibw_corrected')

    combined_means_df <- cbind(first_row_means_df, last_row_means_df)

    #cat("Combined Means DataFrame:\n")
    #print(combined_means_df)

    # Calculate deltas
    combined_means_df$delta_bwbw_wibw_uncorrected <- combined_means_df$bwbw_uncorrected - combined_means_df$wibw_uncorrected
    combined_means_df$delta_bwbw_wibw_corrected <- combined_means_df$bwbw_corrected - combined_means_df$wibw_corrected
    combined_means_df$delta_wiwi_bwwi_uncorrected <- combined_means_df$bwwi_uncorrected - combined_means_df$wiwi_uncorrected
    combined_means_df$delta_wiwi_bwwi_corrected <- combined_means_df$bwwi_corrected - combined_means_df$wiwi_corrected

    # Select only the delta columns and take the absolute values
    delta_means_df <- combined_means_df[, c('delta_bwbw_wibw_uncorrected', 'delta_bwbw_wibw_corrected',
                                            'delta_wiwi_bwwi_uncorrected', 'delta_wiwi_bwwi_corrected')]
    delta_means_df <- abs(delta_means_df)
    delta_means_df <- na.omit(delta_means_df)
    if (nrow(delta_means_df) == 0) {
        cat('No finite delta values were found. Skipping delta PCC plot.\n')
        return(invisible(NULL))
    }

    #cat("\nDelta Means DataFrame:\n")
    #print(delta_means_df)

    if (nrow(delta_means_df) > 2) {
        cat("\nShapiro-Wilk Test Results:\n")
        cat("delta_bwbw_wibw_uncorrected: ", shapiro.test(delta_means_df$delta_bwbw_wibw_uncorrected)$p.value, "\n")
        cat("delta_bwbw_wibw_corrected: ", shapiro.test(delta_means_df$delta_bwbw_wibw_corrected)$p.value, "\n")
        cat("delta_wiwi_bwwi_uncorrected: ", shapiro.test(delta_means_df$delta_wiwi_bwwi_uncorrected)$p.value, "\n")
        cat("delta_wiwi_bwwi_corrected: ", shapiro.test(delta_means_df$delta_wiwi_bwwi_corrected)$p.value, "\n")

        cat("\nDependent T-Test Results:\n")
        t_test_result_bw = t.test(delta_means_df$delta_bwbw_wibw_uncorrected, delta_means_df$delta_bwbw_wibw_corrected, paired = TRUE)$p.value
        t_test_result_wi = t.test(delta_means_df$delta_wiwi_bwwi_uncorrected, delta_means_df$delta_wiwi_bwwi_corrected, paired = TRUE)$p.value
        cat("delta_bwbw_wibw_uncorrected vs delta_bwbw_wibw_corrected: ", t_test_result_bw, "\n")
        cat("delta_wiwi_bwwi_uncorrected vs delta_wiwi_bwwi_corrected: ", t_test_result_wi, "\n")

        p_label1 <- paste("p =", signif(t_test_result_bw, 3))
        p_label2 <- paste("p =", signif(t_test_result_wi, 3))
    } else {
        cat("Not enough data points to perform statistical tests. P values will not be shown.\n")
    }


    with_pdf_defaults(
        file = plot_title,
        width = 4.5,
        height = 4.5,
        font_size = local_font_size,
        font_family = local_font_family,
        plot_fn = function(local_font_size_inner, local_font_family_inner) {
            plot(c(0.5, 4.5), c(0, 0.45), type = 'n', xlab = '', ylab = expression(Delta ~ "mean PCC"), las = 1, xaxt = 'n')
            boxplot(delta_means_df$delta_bwbw_wibw_uncorrected, at = 1, add = TRUE, col = 'gray', yaxt = 'n')
            boxplot(delta_means_df$delta_bwbw_wibw_corrected, at = 2, add = TRUE, col = 'gray', yaxt = 'n')
            boxplot(delta_means_df$delta_wiwi_bwwi_uncorrected, at = 3, add = TRUE, col = 'gray', yaxt = 'n')
            boxplot(delta_means_df$delta_wiwi_bwwi_corrected, at = 4, add = TRUE, col = 'gray', yaxt = 'n')

            if (nrow(delta_means_df) > 2) {
                segments(x0 = 1, y0 = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.2, x1 = 2, y1 = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.2, col = "black", lwd = 0.5)
                segments(x0 = 3, y0 = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.2, x1 = 4, y1 = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.2, col = "black", lwd = 0.5)
                text(x = 1.5, y = mean(delta_means_df$delta_bwbw_wibw_uncorrected) + 0.22, labels = p_label1, xpd = TRUE, srt = 0)
                text(x = 3.5, y = mean(delta_means_df$delta_wiwi_bwwi_uncorrected) + 0.22, labels = p_label2, xpd = TRUE, srt = 0)
            }

            axis(side = 1, at = c(1, 2, 3, 4), labels = c("uncorr.\n", "corr.\n", "uncorr.\n", "corr.\n"), padj = 0.5, tick = FALSE)
            axis(side = 1, at = 0.35, labels = 'Correction\nSample group', padj = 0.5, hadj = 1, tick = FALSE)
            axis(side = 1, at = c(1.5, 3.5), labels = c("\nbetween sample group", "\nwithin sample group"), padj = 0.5, tick = FALSE)
        }
    )
}

save_sample_number_heatmap <- function(df_metadata, font_size = 8, dpi = 300) {
    font_size = resolve_csca_font_size(font_size)
    font_family = resolve_csca_font_family('Helvetica')
    sampled_data <- df_metadata[df_metadata$exclusion == 'no', c('scientific_name', 'sample_group')]
    freq_table <- table(sampled_data$scientific_name, sampled_data$sample_group)
    df_sample_count <- as.data.frame(freq_table)
    names(df_sample_count) <- c('scientific_name', 'sample_group', 'Freq')
    all_scientific_names <- unique(df_metadata$scientific_name)
    all_sample_groups <- unique(df_metadata$sample_group)
    complete_grid <- expand.grid(scientific_name = all_scientific_names, sample_group = all_sample_groups, stringsAsFactors = FALSE)
    df_sample_count <- merge(complete_grid, df_sample_count, by = c("scientific_name", "sample_group"), all.x = TRUE)
    df_sample_count$Freq[is.na(df_sample_count$Freq)] <- 0
    df_sample_count$log2_Freq <- ifelse(df_sample_count$Freq > 0, log2(df_sample_count$Freq + 1), 0)
    fill_max <- max(df_sample_count$log2_Freq)
    freq_breaks <- c(0, 1, 2, 4, 8, 16, 32, 64, 128)
    freq_breaks <- freq_breaks[freq_breaks <= max(df_sample_count$Freq)]
    if (!0 %in% freq_breaks) { freq_breaks <- c(0, freq_breaks) }
    if (!1 %in% freq_breaks) { freq_breaks <- c(freq_breaks, 1) }
    freq_breaks <- sort(unique(freq_breaks))
    if (!is.finite(fill_max) || (fill_max <= 0)) {
        log2_breaks <- c(0)
        fill_max <- 1
    } else {
        log2_breaks <- log2(freq_breaks + 1)
        log2_breaks <- log2_breaks[log2_breaks <= fill_max]
    }
    n_sample_groups <- length(all_sample_groups)
    n_scientific_names <- length(all_scientific_names)
    base_width <- 5
    base_height <- 5
    per_sample_group_width <- 0.1
    per_scientific_name_height <- 0.1
    width <- base_width + (per_sample_group_width * n_sample_groups)
    height <- base_height + (per_scientific_name_height * n_scientific_names)
    p <- ggplot(data = df_sample_count, aes(x = sample_group, y = scientific_name, fill = log2_Freq)) +
        geom_tile(color = "grey80") +
        geom_text(aes(label = Freq), size = font_size / ggplot2::.pt, family = font_family, color = "black") +
        scale_fill_gradientn(
            colors = c("white", "#440154", "#482878", "#3E4989", "#31688E", "#35B779", "#4DCD63", "#76D730", "#B8DE29", "#FDE725"),
            values = sort(unique(c(0, log2_breaks / fill_max, 1))),
            limits = c(0, fill_max),
            name = "# of samples",
            breaks = log2_breaks,
            labels = freq_breaks
        ) +
        xlab('') +
        ylab('') +
        scale_y_discrete(
            limits = rev(unique(df_sample_count$scientific_name))
        ) +
        theme_minimal(base_size = font_size, base_family = font_family) +
        build_standard_ggplot_theme(
            font_size = font_size,
            font_family = font_family,
            x_angle = 90
        )
    ggsave('csca_sample_number_heatmap.pdf', plot = p, width = width, height = height, dpi = dpi)
}

write_ortholog_tables_for_correction = function(correction,
                                                 averaged_orthologs,
                                                 unaveraged_orthologs,
                                                 imputed_averaged_orthologs,
                                                 imputed_unaveraged_orthologs) {
    table_specs = list(
        list(df = averaged_orthologs[[correction]], prefix = 'csca_ortholog_averaged'),
        list(df = unaveraged_orthologs[[correction]], prefix = 'csca_ortholog_unaveraged'),
        list(df = imputed_averaged_orthologs[[correction]], prefix = 'csca_ortholog_averaged.imputed'),
        list(df = imputed_unaveraged_orthologs[[correction]], prefix = 'csca_ortholog_unaveraged.imputed')
    )
    for (spec in table_specs) {
        write_table_with_index_name(
            df = spec[['df']],
            file_path = paste0(spec[['prefix']], '.', correction, '.tsv'),
            index_name = 'target_id',
            sort = FALSE
        )
    }
}

df_og = read.table(file_orthogroup, header = TRUE, sep = '\t', row.names = 1, quote = '', check.names = FALSE)
df_gc = read.table(file_genecount, header = TRUE, sep = '\t', quote = '', check.names = FALSE)
rownames(df_gc) = df_gc[['orthogroup_id']]
df_gc[, 'orthogroup_id'] = NULL
spp_filled = colnames(df_gc)

is_singlecopy = get_singlecopy_bool_index(df_gc, spp_filled)
df_singleog = df_og[is_singlecopy, spp_filled, drop = FALSE]
spp = sub('_', ' ', spp_filled)
df_metadata = prepare_metadata_table(dir_csca_input_table, selected_sample_groups, spp)
label_orders = get_label_orders(df_metadata)
df_color_averaged = get_df_labels_averaged(df_metadata, label_orders, selected_sample_groups, sample_group_colors)
df_color_unaveraged = get_df_labels_unaveraged(df_metadata, selected_sample_groups, sample_group_colors)
cat('Number of orthologs in input table:', nrow(df_og), '\n')
cat('Number of selected single-copy orthologs:', nrow(df_singleog), '\n')
cat('Number of selected species:', length(spp), '\n')
save_sample_number_heatmap(df_metadata, font_size = font_size, dpi = 300)

unaveraged_tcs = load_unaveraged_expression_tables(dir_csca_input_table, spp_filled, batch_effect_alg)
unaveraged_tcs = extract_selected_tc_only(unaveraged_tcs, df_metadata)

# if a species was skipped during load_unaveraged_expression_tables(), it will cause indexing issues down the line, due to df_singleog having more species entries than the tcs
if (length(unaveraged_tcs[['corrected']]) < length(colnames(df_singleog))) {
    df_singleog = df_singleog[, names(unaveraged_tcs[['corrected']]), drop = FALSE]
}

unaveraged_orthologs = extract_ortholog_unaveraged_expression_table(df_singleog, unaveraged_tcs)
averaged_tcs = unaveraged2averaged(unaveraged_tcs, df_metadata, selected_sample_groups)
averaged_orthologs = extract_ortholog_mean_expression_table(df_singleog, averaged_tcs, label_orders)
write_pivot_table(df_metadata, unaveraged_tcs, selected_sample_groups)

cat('Applying expression level imputation for missing orthologs.\n')
quiet_impute_expression = function(dat) {
    withCallingHandlers(
        impute_expression(dat, strategy = missing_strategy),
        warning = function(w) {
            if (grepl('Precision for components .*below \\.Machine\\$double\\.eps', conditionMessage(w))) {
                invokeRestart('muffleWarning')
            }
        }
    )
}
imputed_averaged_orthologs = list()
imputed_unaveraged_orthologs = list()
for (d in correction_labels) {
    imputed_averaged_orthologs[[d]] = quiet_impute_expression(averaged_orthologs[[d]])
    imputed_unaveraged_orthologs[[d]] = quiet_impute_expression(unaveraged_orthologs[[d]])
    write_ortholog_tables_for_correction(
        correction = d,
        averaged_orthologs = averaged_orthologs,
        unaveraged_orthologs = unaveraged_orthologs,
        imputed_averaged_orthologs = imputed_averaged_orthologs,
        imputed_unaveraged_orthologs = imputed_unaveraged_orthologs
    )
}
cat(nrow(imputed_unaveraged_orthologs[['corrected']]), 'orthologs were found after filtering and imputation.\n')
df_metadata = calculate_correlation_within_group(unaveraged_orthologs, averaged_orthologs, df_metadata, selected_sample_groups)
averaged_plot_cache = build_averaged_plot_cache(imputed_averaged_orthologs)
tryCatch(save_group_cor_scatter(df_metadata, font_size = font_size), error = function(e) cat("Warning: save_group_cor_scatter skipped:", conditionMessage(e), "\n"))
tryCatch(save_group_cor_histogram(df_metadata, df_color_unaveraged, font_size = font_size), error = function(e) cat("Warning: save_group_cor_histogram skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_tsne_plot(tc = imputed_unaveraged_orthologs[['corrected']], df_label = df_color_unaveraged), error = function(e) cat("Warning: save_averaged_tsne_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_heatmap_plot(imputed_averaged_orthologs, df_color_averaged, averaged_plot_cache = averaged_plot_cache), error = function(e) cat("Warning: save_averaged_heatmap_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_dendrogram_plot(imputed_averaged_orthologs, df_color_averaged, averaged_plot_cache = averaged_plot_cache), error = function(e) cat("Warning: save_averaged_dendrogram_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_dimensionality_reduction_summary(imputed_averaged_orthologs, df_color_averaged, averaged_plot_cache = averaged_plot_cache), error = function(e) cat("Warning: save_averaged_dimensionality_reduction_summary skipped:", conditionMessage(e), "\n"))
tryCatch(save_averaged_box_plot(imputed_averaged_orthologs, df_color_averaged, averaged_plot_cache = averaged_plot_cache), error = function(e) cat("Warning: save_averaged_box_plot skipped:", conditionMessage(e), "\n"))
tryCatch({df_metadata = save_unaveraged_pca_plot(imputed_unaveraged_orthologs, df_color_unaveraged, df_metadata)}, error = function(e) cat("Warning: save_unaveraged_pca_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_unaveraged_tsne_plot(imputed_unaveraged_orthologs, df_color_unaveraged), error = function(e) cat("Warning: save_unaveraged_tsne_plot skipped:", conditionMessage(e), "\n"))
tryCatch(save_delta_pcc_plot(directory = dir_csca_input_table, plot_title = 'csca_delta_pcc_boxplot.pdf'), error = function(e) cat("Warning: save_delta_pcc_plot skipped:", conditionMessage(e), "\n"))

file_metadata_out = file.path(dir_csca, 'metadata.tsv')
write.table(df_metadata, file_metadata_out, row.names = FALSE, sep = '\t', quote = FALSE)

cat(sprintf('Number of SRA samples for exclusion potting: %s\n', formatC(nrow(df_metadata), format = 'd', big.mark = ',')))
out_path = file.path(dir_csca, 'csca_exclusion.pdf')
save_exclusion_plot(df = df_metadata, out_path = out_path, font_size = font_size, y_label = "Sample count")

if (file.exists('Rplots.pdf')) {
    file.remove('Rplots.pdf')
}
cat('csca.r completed!\n')
