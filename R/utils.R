#' Construct a Combined SNN Graph Integrating Gene Expression and Spatial Coordinates
#'
#' This function computes a combined shared nearest neighbor (SNN) graph
#' by integrating gene expression-based SNN from Seurat with spatial
#' coordinate-based SNN. The final adjacency matrix is a convex combination
#' of the two, controlled by parameter `alpha`.
#'
#' @param seurat_obj A Seurat object containing single-cell or spatial data.
#' @param alpha Numeric scalar in [0,1]. Weight of spatial coordinate-based
#'        SNN in the final graph. (1 - alpha) is the weight of the gene expression SNN.
#' @param position A matrix or data frame of spatial coordinates, with rownames
#'        corresponding to cell/spot names and columns typically being `x, y`
#'        (or higher-dimensional embeddings).
#' @param snn Gene expression-based SNN graph. Default: `seurat_obj@graphs$SCT_snn`.
#'
#' @return A sparse adjacency matrix representing the combined SNN graph.
#'
#' @export
getCOSNN <- function(seurat_obj, alpha, position, snn = seurat_obj@graphs$SCT_snn){
  sp_obj <- seurat_obj
  spatial_coords <- position[colnames(sp_obj), ]

  coord_snn <- as.matrix(FindNeighbors(as.matrix(spatial_coords), compute.SNN = TRUE)$snn)
  ge_snn <- as.matrix(snn)


  alpha <- alpha
  combined_adj <- Matrix::Matrix((1 - alpha) * ge_snn + alpha * coord_snn)
  return(combined_adj)

}


#' Generate Data Matrix from Phenotype-associated Peaks
#'
#' This function quantifies accessibility of differential accessible regions (DARs)
#' in single-cell ATAC-seq data and constructs a Seurat object containing a ChromatinAssay.
#' Peaks can be provided as character strings in the format "chr:start-end".
#'
#' @param peaks Character vector of peak coordinates.
#' @param seq Character vector of length 2 specifying separators, e.g., c("-", "-").
#' @param fragmentsPath Path to a tabix-indexed single-cell fragment file.
#' @param cells Character vector of cell barcodes to include.
#' @param min.cells Minimum number of cells required for each peak (default: 3).
#' @return A Seurat object containing the ATAC assay with peak x cell counts.
#'
#' @export
runscCountsDARs <- function(peaks, seq, fragmentsPath, cells, min.cells = 3){

  DARs_GR <- makeGRangesFromDataFrame(data.frame(t(matrix(unlist(strsplit(unlist(strsplit(peaks,seq[1])),seq[2])),3))),
                                      seqnames.field="X1",
                                      start.field="X2",end.field="X3"
  )

  sc_bulk_peak_counts <- FeatureMatrix(
    fragments = CreateFragmentObject(
      path = fragmentsPath
    ),
    features = DARs_GR,
    cells = cells)

  sc_bulk_seurat <- CreateSeuratObject(
    counts = CreateChromatinAssay(
      counts = sc_bulk_peak_counts
    ),
    assay = "ATAC",min.cells = min.cells)

  return(sc_bulk_seurat)
}


#' Sigmoid-based feature compression
#'
#' This function applies a standardized sigmoid transformation to a numeric vector.
#' The input is first z-score normalized (centered by mean and scaled by standard deviation),
#' then passed through the logistic sigmoid function to map values into the (0, 1) range.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector of the same length as `x`, with values compressed into (0, 1).
#'
#'
sigmoid_compress <- function(x) {

  1 / (1 + exp(-(x-mean(x))/sd(x)))

}

