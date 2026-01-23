library(dplyr)
library(ggplot2)
library(ggpubr)

#' Plot a gene by cell type with group comparison, optionally with pseudobulk per sample
#'
#' This function produces:
#' 1) A cell level plot: per cell gene values grouped by cell type and split by a comparison factor
#' 2) An optional pseudobulk plot: per sample mean gene values within each (cell type, comparison) group
#'
#' Statistical testing:
#' Performs a within cell type comparison between levels of compare_col using ggpubr::stat_compare_means.
#' The p value labels are displayed as custom significance symbols via symnum.args.
#'
#' @param df Data frame containing gene expression and metadata columns.
#' @param gene String. Column name in df containing the gene expression values.
#' @param celltype_col String. Column name in df defining cell type labels for the x axis.
#' @param sample_col String. Column name in df defining biological sample identifiers.
#' @param compare_col String. Column name in df defining the grouping variable to compare
#'   within each cell type (for example location, genotype, treatment).
#' @param compare_levels Optional character vector. If provided, restrict to these levels and
#'   set factor order accordingly. Useful to enforce comparisons like c("lung","tibia").
#' @param do_pseudobulk Logical. If TRUE, compute per sample mean gene values and produce a pseudobulk plot.
#' @param stat_method String. Statistical test passed to ggpubr::stat_compare_means (default wilcox.test).
#' @param sym_cutpoints Numeric vector. Cutpoints for mapping p values to custom symbols.
#'   Must be increasing and include 0 and Inf.
#' @param sym_symbols Character vector. Symbols corresponding to sym_cutpoints intervals.
#'   Length must be length(sym_cutpoints) - 1.
#' @param jitter_width Numeric. Horizontal jitter used for cell level points.
#' @param dodge_width Numeric. Dodge width used to separate comparison groups within each cell type.
#'
#' @return A named list containing:
#'   - cell_level_plot: ggplot object for cell level visualization
#'   - cell_level_df: cleaned and standardized data used for the cell level plot
#'   - pseudobulk_plot: ggplot object for pseudobulk visualization (only if do_pseudobulk = TRUE)
#'   - pseudobulk_df: summarized pseudobulk data (only if do_pseudobulk = TRUE)
#'
#' @examples
#' res <- plot_gene_by_celltype_compare(
#'   df = spp1_df,
#'   gene = "Spp1",
#'   celltype_col = "Ann_Level3",
#'   sample_col = "sample_name",
#'   compare_col = "location",
#'   compare_levels = c("lung","tibia")
#' )
#' res$cell_level_plot
#' res$pseudobulk_plot
plot_gene_by_celltype_compare <- function(
  df,
  gene,
  celltype_col,
  sample_col,
  compare_col,
  compare_levels = NULL,
  do_pseudobulk = TRUE,
  stat_method = "wilcox.test",
  sym_cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  sym_symbols   = c("''", "'", "**", "*", "ns"),
  jitter_width = 0.25,
  dodge_width = 0.8
) {
  # ---- Input validation ----------------------------------------------------
  # Ensure required columns exist before doing any work so errors are immediate and clear.
  needed <- c(gene, celltype_col, sample_col, compare_col)
  missing_cols <- setdiff(needed, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in df: ", paste(missing_cols, collapse = ", "))
  }

  # Validate significance mapping configuration.
  # sym_cutpoints should define bins and sym_symbols should define labels for those bins.
  if (!is.numeric(sym_cutpoints) || length(sym_cutpoints) < 2) {
    stop("sym_cutpoints must be a numeric vector of length >= 2.")
  }
  if (length(sym_symbols) != length(sym_cutpoints) - 1) {
    stop("sym_symbols must have length length(sym_cutpoints) - 1.")
  }
  if (any(diff(sym_cutpoints) <= 0)) {
    stop("sym_cutpoints must be strictly increasing.")
  }

  # ---- Standardize the input data -----------------------------------------
  # Rename user supplied columns to stable internal names so downstream code stays readable.
  # This also isolates the plotting code from changes in the original column names.
  df2 <- df %>%
    dplyr::select(
      gene_val = all_of(gene),
      celltype = all_of(celltype_col),
      sample   = all_of(sample_col),
      compare  = all_of(compare_col)
    ) %>%
    # Drop rows with missing values in any required field.
    dplyr::filter(!is.na(gene_val), !is.na(celltype), !is.na(sample), !is.na(compare))

  # If the user specified compare_levels, restrict to those levels and enforce their order.
  # This avoids accidental issues from typos or inconsistent factor ordering.
  if (!is.null(compare_levels)) {
    df2 <- df2 %>% dplyr::filter(compare %in% compare_levels)
    df2$compare <- factor(df2$compare, levels = compare_levels)
  } else {
    df2$compare <- as.factor(df2$compare)
  }

  # Build the configuration list for ggpubr significance symbol mapping.
  sym_args <- list(cutpoints = sym_cutpoints, symbols = sym_symbols)

  # Create a human readable caption from the cutpoints and symbols.
  # This assumes the cutpoints follow the common pattern 0 < ... < Inf.
  caption_txt <- paste0(
    "Significance codes: ",
    sym_symbols[1], " p \u2264 ", sym_cutpoints[2],
    "; ", sym_symbols[2], " p \u2264 ", sym_cutpoints[3],
    "; ", sym_symbols[3], " p \u2264 ", sym_cutpoints[4],
    "; ", sym_symbols[4], " p \u2264 ", sym_cutpoints[5],
    "; ", sym_symbols[5], " p > ",  sym_cutpoints[5]
  )

  # ---- Cell level plot -----------------------------------------------------
  # Plot per cell values:
  # - X axis is cell type
  # - Groups are compare levels (dodged within each cell type)
  # - Points are drawn first so they appear behind boxplots
  p_cells <- ggplot(
    df2,
    aes(x = celltype, y = gene_val, fill = compare, color = compare)
  ) +
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = jitter_width,
        dodge.width = dodge_width
      ),
      alpha = 0.6,
      size = 1.2
    ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.35,
      width = 0.7,
      position = position_dodge(width = dodge_width),
      # Force a consistent outline color for boxplots while preserving fill by group.
      color = "black"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      x = "Cell type",
      y = paste0(gene, " expression"),
      fill = compare_col,
      color = compare_col
    ) +
    # Run within cell type comparisons between compare levels.
    # This tests compare groups separately for each x category (celltype).
    stat_compare_means(
      aes(group = compare),
      method = stat_method,
      label = "p.signif",
      hide.ns = FALSE,
      symnum.args = sym_args
    ) +
    labs(caption = caption_txt)

  # Package outputs in a named list for convenient access.
  out <- list(cell_level_plot = p_cells, cell_level_df = df2)

  # ---- Pseudobulk plot (optional) -----------------------------------------
  # Pseudobulk here means mean gene value per sample within each (celltype, compare) group.
  # This reduces pseudoreplication and makes each point represent a biological sample.
  if (isTRUE(do_pseudobulk)) {
    pseudo_df <- df2 %>%
      group_by(sample, celltype, compare) %>%
      summarise(
        gene_val = mean(gene_val, na.rm = TRUE),
        .groups = "drop"
      )

    # Plot pseudobulk per sample values with the same grouping structure as the cell level plot.
    p_pseudo <- ggplot(
      pseudo_df,
      aes(x = celltype, y = gene_val, fill = compare, color = compare)
    ) +
      geom_jitter(
        position = position_jitterdodge(
          jitter.width = 0.15,
          dodge.width = dodge_width
        ),
        size = 2,
        alpha = 0.8
      ) +
      geom_boxplot(
        outlier.shape = NA,
        alpha = 0.35,
        width = 0.7,
        position = position_dodge(width = dodge_width),
        color = "black"
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(
        x = "Cell type",
        y = paste0("Pseudobulk ", gene, " expression (mean per sample)"),
        fill = compare_col,
        color = compare_col
      ) +
      stat_compare_means(
        aes(group = compare),
        method = stat_method,
        label = "p.signif",
        hide.ns = FALSE,
        symnum.args = sym_args
      ) +
      labs(caption = caption_txt)

    out$pseudobulk_plot <- p_pseudo
    out$pseudobulk_df <- pseudo_df
  }

  return(out)
}

