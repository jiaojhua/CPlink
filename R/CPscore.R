
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
                              method = c("cosine", "pearcon", "spearman", "euclidean")[1]){

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


