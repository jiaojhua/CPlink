
#' Run MCA on Phenotype-specific Feature to Construct Latent Phenotype Space
#'
#' This function performs Multiple Correspondence Analysis (MCA) on a given
#' set of phenotype-specific features to generate low-dimensional embeddings of cells/spots.
#'
#' @param Mtx A numeric data matrix (features x cells/spots).
#' @param markers Character vector of phenotype-specific feature (genes or peaks) to be used.
#' @param reduction_name Character, the name assigned to the MCA reduction (default: "mca").
#'
#' @return A list containing:
#' \itemize{
#'   \item embeddings: Cell/spot coordinates in MCA space.
#'   \item features: Features used in the MCA.
#' }
#'
#' @export
runMCA <- function(Mtx = Mtx,
                   markers = markers,
                   reduction_name = "mca"){


  mat <- as.matrix(Mtx)
  markers <- markers
  reduction_name <- reduction_name

  # intersect markers with rownames
  markers.use <- intersect(markers, rownames(mat))
  data_matrix <- mat[markers.use,]

  # These functions are based on https://github.com/RausellLab/CelliD
  MCA_output <- RunMCA(X = data_matrix,
                       features = markers.use)

  return(MCA_output)
}


#' Compute marker–cell/spot association matrix in the latent phenotype space
#'
#' Compute similarity scores between marker features and cells/spots using coordinates
#' in the latent phenotype space (MCA). The resulting marker × cell/spot matrix quantifies
#' feature–cell associations.
#'
#' @param MCA_obj A list-like object representing the latent phenotype space with at least:
#'        \itemize{
#'          \item featuresCoordinates: numeric matrix (features × dims) of co-embedded marker coordinates
#'          \item cellsCoordinates: numeric matrix (cells/spots × dims) of cell/spot embeddings in the latent space
#'          Row and column names should be present and consistent with marker and cell/spot identifiers.
#'        }
#' @param markers Character vector of marker features ( genes / peaks) used to define the phenotype-specific signatures.

#' @param ues_dim Integer scalar; number of latent dimensions to use from the MCA embedding (default: 50).
#' @param method Character; similarity metric to quantify feature–cell association. Supported options:
#'        \itemize{
#'          \item "cosine": Cosine similarity (default)
#'          \item "pearson": Pearson correlation
#'          \item "spearman": Spearman correlation
#'          \item "euclidean": Inverse Euclidean distance (1/(1+distance))
#'        }
#'
#' @return A numeric matrix of association scores with markers as rows and cells/spots as columns.
#'    Values reflect the degree of association between each marker feature and each cell/spot
#'   in the latent phenotype space; this matrix can be used as input for diffusion-based
#'   smoothing and permutation testing to derive phenotype risk scores.
#'
#' @export
GetCellMarkerDist <- function(MCA_obj = MCA_obj,
                              markers = markers,
                              ues_dim = 50,
                              method = c("cosine", "pearson", "spearman", "euclidean")[1]){

  # Validate inputs (latent phenotype space must provide co-embedded coordinates)
  if (is.null(MCA_obj) || !is.list(MCA_obj)) {

    stop("MCA_obj must be a list-like latent phenotype object containing 'featuresCoordinates' and 'cellsCoordinates'.")

  }
  if (is.null(MCA_obj$featuresCoordinates) || is.null(MCA_obj$cellsCoordinates)) {

    stop("MCA_obj must include 'featuresCoordinates' and 'cellsCoordinates'.")

  }

  markers.update <- intersect(markers,
                              rownames(MCA_obj$featuresCoordinates))

  markersCoor <- MCA_obj$featuresCoordinates[markers.update, 1:ues_dim]
  cellsCoor <- MCA_obj$cellsCoordinates[,1:ues_dim]

  if (ncol(markersCoor) != ncol(cellsCoor)) {

    stop("featuresCoordinates and cellsCoordinates must have the same latent dimensionality (columns).")

  }

  if(method == "cosine"){

    mtx <- proxyC::simil(markersCoor, cellsCoor,
                         method = "cosine", drop0=TRUE)

    return(mtx)
  }

  if(method == "pearson"){

    mtx <- cor(t(markersCoor), t(cellsCoor),
               method = c("pearson", "spearman")[1])

    return(mtx)
  }

  if(method == "spearman"){

    mtx <- cor(t(markersCoor), t(cellsCoor),
               method = c("pearson", "spearman")[2])

    return(mtx)
  }

  if(method == "euclidean"){

    mtx <- 1/(1+rdist::cdist(markersCoor, cellsCoor))
    colnames(mtx) <- rownames(cellsCoor)
    rownames(mtx) <- rownames(markersCoor)

    return(mtx)
  }

}


#' Compute Weighted Cell/Spot Scores Based on Phenotype-Associated Features
#'
#' This function calculates initial cell/spot-phenotype association scores by
#' aggregating phenotype-associated feature-cell associations derived from
#' single-cell or spatial omics data.
#'
#' @param dist_mat A matrix (features x cells/spots) where rows are features
#'        (genes/peaks) and columns are cells/spots. Higher values indicate
#'        stronger feature-cell/spot associations.
#' @param marker_mat Data frame or matrix of phenotype-associated markers
#' identified  from reference (bulk) data. Must include:
#'        \itemize{
#'          \item Rows: Marker features (must overlap with dist_mat rownames)
#'          \item Column specified by direction_col: signed values representing
#'                effect sizes or directionality (e.g., log fold change)
#'        }
#' @param up_weight Numeric scalar specifying the contribution of upregulated markers
#'        to the final cell score (default: 0.5)
#' @param down_weight Numeric scalar specifying the contribution of downregulated markers
#'        to the final cell score (default: 0.5)
#' @param direction_col Column name in marker_mat containing directional
#'        information (default: "logFC")
#'
#' @return Numeric vector of length equal to the number of cells or spots,
#'         representing phenotype-association scores for each cell or spot.
#'
#' @export
runCellScores <- function(dist_mat,
                          marker_mat,
                          up_weight = 0.5,
                          down_weight = 0.5,
                          direction_col = "logFC"){

  common_features <- intersect(rownames(dist_mat), rownames(marker_mat))
  if (length(common_features) == 0) {

    stop("No common features between distance matrix and marker matrix")

  }

  up_markers <- common_features[marker_mat[common_features, direction_col] > 0]
  down_markers <- common_features[marker_mat[common_features, direction_col] < 0]

  if (length(up_markers) == 0) warning("No upregulated markers found")
  if (length(down_markers) == 0) warning("No downregulated markers found")

  up_weights <- abs(marker_mat[up_markers, direction_col])
  down_weights <- abs(marker_mat[down_markers, direction_col])

  if (length(up_markers) > 0) up_weights <- up_weights/sum(up_weights)
  if (length(down_markers) > 0) down_weights <- down_weights/sum(down_weights)


  up_scores <- if (length(up_markers) > 0) {

    colSums(dist_mat[up_markers, , drop = FALSE] * up_weights, na.rm = TRUE)

  } else { 0 }

  down_scores <- if (length(down_markers) > 0) {

    colSums(dist_mat[down_markers, , drop = FALSE] * down_weights, na.rm = TRUE)

  } else { 0 }

  cell_scores <- up_weight * up_scores - down_weight * down_scores

  return(cell_scores)
}


#' Transfer Phenotype Association Scores from Spatial Transcriptomics to CODEX Cells
#'
#' This function computes initial phenotype association scores from spatial
#' transcriptomics (ST) data and transfers these scores to spatially aligned
#' CODEX cells. It extends the CPlink framework to spatial proteomics by
#' leveraging phenotype-associated molecular features derived from reference
#' bulk omics data and their spatial correspondence with CODEX cells.
#'
#' Specifically, the function first constructs a latent phenotype space using
#' Multiple Correspondence Analysis (MCA) restricted to phenotype-associated
#' features. It then quantifies feature-to-cell associations in this space
#' and calculates weighted phenotype association scores for ST cells/spots.
#' Finally, these scores are transferred to CODEX cells through spatial
#' nearest-neighbor mapping with Gaussian kernel weighting.
#'
#' @param STMtx A numeric count matrix (features x cells/spots) from spatial
#'   transcriptomics data. Rows represent genes and columns represent
#'   cells/spots. Must contain row names and column names.
#' @param DEGs_mtx A data.frame containing phenotype-associated molecular
#'   features and their corresponding effect weights. Must contain:
#'   \itemize{
#'     \item markers: Character vector of phenotype-associated feature IDs.
#'     \item logFC: Signed effect weights indicating the direction and
#'       magnitude of phenotype associations.
#'   }
#'   Feature identifiers must also be provided as row names for compatibility
#'   with the downstream scoring function.
#' @param MCA_dim Integer scalar. Number of MCA dimensions retained for
#'   calculating feature-to-cell associations (default: 50).
#' @param sim_methods Character string specifying the similarity metric
#'   used to quantify associations between molecular features and ST cells/spots.
#'   Supported methods include "cosine", "pearson", "spearman", and
#'   "euclidean" (default: "cosine").
#' @param up_weight Numeric scalar specifying the contribution of positively
#'   associated features to the initial phenotype score (default: 0.5).
#' @param down_weight Numeric scalar specifying the contribution of negatively
#'   associated features to the initial phenotype score (default: 0.5).
#' @param k Integer scalar. Number of spatially nearest ST cells/spots used
#'   to estimate the phenotype association score of each CODEX cell
#'   (default: 20).
#' @param ST_coords A matrix or data.frame containing spatial coordinates
#'   of ST cells/spots. Must contain columns "x" and "y".
#'   Row order must correspond to the ST cell/spot score vector.
#' @param SP_coords A matrix or data.frame containing spatial coordinates
#'   of CODEX cells. Must contain columns "x" and "y", with row names
#'   representing CODEX cell identifiers.
#'
#' @return A named numeric vector containing the transferred initial
#'   phenotype association score for each CODEX cell. Vector names
#'   correspond to CODEX cell identifiers.
#'
#' @seealso runMCA, GetCellMarkerDist, runCellScores,
#'   mapSPTScoreToCODEX, runCPlink_CODEX
#'
#' @export
getCodexInitCellScore <- function(STMtx,
                                  DEGs_mtx,
                                  MCA_dim = 50,
                                  sim_methods = c("cosine", "pearson", "spearman", "euclidean")[1],
                                  up_weight = 0.5,
                                  down_weight = 0.5,
                                  k = 20,
                                  ST_coords,
                                  SP_coords


){
  if (is.null(DEGs_mtx$markers)) stop("`DEGs_mtx` must contain element `markers` (character vector of feature IDs).")
  if (!sim_methods %in% c("cosine","pearson","spearman","euclidean")) stop("`sim_methods` must be one of 'cosine','pearson','spearman','euclidean'.")

  message("Starting CPlink pipeline")

  # ---- MCA embedding restricted to phenotype-specific features ----
  message("Running MCA on phenotype-specific features")
  MCA_obj <- runMCA(Mtx = STMtx, markers = DEGs_mtx$markers)

  # ---- Compute feature-to-cell association matrix (features x cells) ----
  cell_feature_mtx <- GetCellMarkerDist(MCA_obj = MCA_obj, markers = DEGs_mtx$markers, ues_dim = MCA_dim,
                                        method = sim_methods)

  # ---- Compute initial per-cell/spot score using marker weights ----
  message("Calculating initial phenotype association scores")
  ST_CellScore <- runCellScores(cell_feature_mtx,  DEGs_mtx, up_weight = up_weight, down_weight = down_weight, direction_col = "logFC")

  CODEX_CellScore <- mapSPTScoreToCODEX(ST_coords, SP_coords, ST_CellScore, k = k)

  return(CODEX_CellScore)
}


#' Transfer Spatial Transcriptomics Phenotype Scores to CODEX Cells
#'
#' This function transfers phenotype association scores from spatial
#' transcriptomics (ST) cells/spots to spatially aligned CODEX cells using
#' a spatial k-nearest-neighbor (kNN) mapping strategy with Gaussian
#' kernel weighting.
#'
#' For each CODEX cell, the function identifies its k nearest ST cells/spots
#' based on their spatial coordinates. The phenotype association scores of
#' these neighboring ST cells/spots are then aggregated using distance-dependent
#' Gaussian weights, assigning greater contributions to spatially closer
#' cells/spots.
#'
#' This approach enables the transfer of transcriptomics-derived phenotype
#' information to spatial proteomics data through shared spatial organization,
#' without requiring direct correspondence between transcriptomic and
#' proteomic molecular features.
#'
#' @param STcoord A matrix or data.frame containing spatial coordinates
#'   of ST cells/spots. Rows represent ST cells/spots and columns contain
#'   spatial coordinates (e.g., "x" and "y").
#' @param CODEXcoord A matrix or data.frame containing spatial coordinates
#'   of CODEX cells. Rows represent CODEX cells and columns contain
#'   spatial coordinates (e.g., "x" and "y").
#' @param STCellScore A numeric vector containing initial phenotype
#'   association scores for ST cells/spots. The order of scores must
#'   correspond to the row order of STcoord.
#' @param k Integer scalar. Number of nearest ST cells/spots used for
#'   estimating each CODEX cell score (default: 20).
#' @param sigma Numeric scalar or NULL. Bandwidth parameter controlling
#'   the spatial decay of Gaussian kernel weights. If NULL, sigma is
#'   estimated as the median of nonzero distances among the identified
#'   ST-to-CODEX nearest-neighbor pairs (default: NULL).
#' @param x_col Character string specifying the column containing
#'   x coordinates (default: "x").
#' @param y_col Character string specifying the column containing
#'   y coordinates (default: "y").
#'
#' @details
#' For each CODEX cell i, the transferred phenotype score is computed as:
#'
#'   Score_i = sum_j(w_ij * Score_j) / sum_j(w_ij)
#'
#' where Score_j denotes the initial phenotype association score of
#' neighboring ST cell/spot j, and w_ij is the Gaussian kernel weight:
#'
#'   w_ij = exp(-d_ij^2 / (2 * sigma^2))
#'
#' Here, d_ij represents the Euclidean distance between CODEX cell i
#' and ST cell/spot j in the shared spatial coordinate system.
#'
#' If sigma cannot be estimated from nonzero neighbor distances,
#' it is set to 1. Missing or undefined transferred scores are replaced
#' with the median of the computed CODEX scores.
#'
#' Spatial coordinates from both modalities must be registered to a
#' common coordinate system before applying this function.
#'
#' @return A named numeric vector containing the spatially transferred
#'   phenotype association scores for CODEX cells. The vector names
#'   correspond to the row names of CODEXcoord.
#'
#' @seealso getCodexInitCellScore, runCPlink_CODEX
#'
#' @export
mapSPTScoreToCODEX <- function(STcoord,
                               CODEXcoord,
                               STCellScore,
                               k = 20,
                               sigma = NULL,
                               x_col = "x",
                               y_col = "y") {

  # 坐标矩阵
  st_xy <- (STcoord[, c(x_col, y_col), drop = FALSE])
  codex_xy <- (CODEXcoord[, c(x_col, y_col), drop = FALSE])

  coor_ST_SP_nn <- FNN::get.knnx(
    data = st_xy,
    query = codex_xy,
    k = k
  )

  if (is.null(sigma)) {
    sigma <- median(coor_ST_SP_nn$nn.dist[coor_ST_SP_nn$nn.dist > 0], na.rm = TRUE)
  }

  if (is.na(sigma) || sigma <= 0) {
    sigma <- 1
  }

  w <- exp(-(coor_ST_SP_nn$nn.dist^2) / (2 * sigma^2))

  score_mat <- matrix(
    STCellScore[coor_ST_SP_nn$nn.index],
    nrow = nrow(codex_xy),
    ncol = k
  )

  denom <- rowSums(w)
  codex_score <- rowSums(w * score_mat) / denom

  codex_score[is.na(codex_score) | denom == 0] <- median(codex_score, na.rm = TRUE)

  names(codex_score) <- rownames(codex_xy)

  return(codex_score)
}
