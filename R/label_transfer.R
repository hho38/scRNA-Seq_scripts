
library(Seurat)

#' Transfer reference labels to a query Seurat object
#'
#' Assign cell labels from a reference Seurat object to a query Seurat object
#' using Seurat's anchor-based label transfer workflow. This function identifies
#' anchors between the reference and query with `FindTransferAnchors()`,
#' predicts labels for each query cell with `TransferData()`, and writes the
#' predicted labels directly into the query object's metadata.
#'
#' The reference object must already contain a PCA reduction named `"pca"`,
#' and `label_col` must be present in the reference object's metadata.
#'
#' @param obj_of_interest A Seurat object containing the query cells to annotate.
#'   The returned object will contain the transferred labels in its metadata.
#' @param reference_obj A Seurat object containing the annotated reference cells.
#'   This object should already be processed for label transfer and must include
#'   a PCA reduction named `"pca"`.
#' @param label_col A character string giving the name of the metadata column in
#'   `reference_obj` that contains the labels to transfer (for example,
#'   `"celltype"` or `"seurat_annotations"`).
#' @param metadata_col A character string giving the name of the metadata column
#'   to create in `obj_of_interest` for the transferred labels.
#' @param dims An integer vector of dimensions to use for anchor finding and
#'   label transfer. Defaults to `1:30`.
#' @param include_scores Logical. If `TRUE`, also adds a metadata column named
#'   `paste0(metadata_col, "_score")` containing the maximum prediction score
#'   for each transferred label. If `FALSE`, only the transferred labels are
#'   added. Defaults to `FALSE`.
#'
#' @return A Seurat object identical to `obj_of_interest` but with a new metadata
#'   column containing transferred labels. If `include_scores = TRUE`, an
#'   additional score column is also added.
#'
#' @details
#' Internally, this function performs three steps:
#' 1. Finds transfer anchors between `reference_obj` and `obj_of_interest`.
#' 2. Transfers labels from `reference_obj[[label_col]]` to the query cells.
#' 3. Adds the transferred labels to the query object's metadata.
#'
#' When `include_scores = TRUE`, only the top prediction score
#' (`prediction.score.max`) is retained. Per-class prediction score columns
#' returned by `TransferData()` are discarded.
#'
#' @examples
#' query_obj <- label_transfer_simple(
#'   obj_of_interest = query_obj,
#'   reference_obj = ref_obj,
#'   label_col = "celltype",
#'   metadata_col = "predicted_celltype"
#' )
#'
#' query_obj <- label_transfer_simple(
#'   obj_of_interest = query_obj,
#'   reference_obj = ref_obj,
#'   label_col = "celltype",
#'   metadata_col = "predicted_celltype",
#'   include_scores = TRUE
#' )
#'
#' head(query_obj[[]][, c("predicted_celltype", "predicted_celltype_score")])
label_transfer_simple <- function(
  obj_of_interest,
  reference_obj,
  label_col,
  metadata_col,
  dims = 1:30,
  include_scores = FALSE
) {
  anchors <- FindTransferAnchors(
    reference = reference_obj,
    query = obj_of_interest,
    dims = dims,
    reference.reduction = "pca"
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_obj[[label_col, drop = TRUE]],
    dims = dims
  )

  colnames(predictions)[colnames(predictions) == "predicted.id"] <- metadata_col

  if (include_scores) {
    colnames(predictions)[colnames(predictions) == "prediction.score.max"] <- paste0(metadata_col, "_score")
    keep_cols <- c(metadata_col, paste0(metadata_col, "_score"))
  } else {
    keep_cols <- metadata_col
  }

  predictions <- predictions[, keep_cols, drop = FALSE]

  obj_of_interest <- AddMetaData(obj_of_interest, metadata = predictions)
  obj_of_interest
}

