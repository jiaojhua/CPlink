#' CPlink: identify phenotype-associated cells/spots from target single-cell / spatial data
#'
#' This function executes the core CPlink pipeline for transferring phenotype
#' knowledge (features and weights) learned from reference bulk data to target
#' single-cell or spatial data and then identifying phenotype-associated
#' cells/spots via diffusion (Random Walk with Restart) in permutation testing
#' framework.
#'
#' @param Mtx A count matrix (features x cells). Rows are feature identifiers
#'   (genes/peaks) and columns are cell/spot identifiers. Must have rownames and colnames.
#' @param DEGs_mtx A data.frame containing phenotype-specific features and effect weights with:
#'        \itemize{
#'          \item markers: Vector of phenotype-specific features
#'          \item logFC: Column name containing log fold changes (effect weights)
#'        }
#' @param snn Sparse similarity matrix (e.g., SNN graph) for cell-cell (spot-spot) interactions.
#'    Typically the SNN graph produced by Seurat or a custom nearest-neighbor graph.
#'   Row/column names must match the column names of `Mtx`. Accepts dense or dgCMatrix.
#' @param MCA_dim Integer scalar. Number of MCA components (dimensions) to retain when
#'   computing feature-cell associations (default: 50).
#' @param sim_methods Similarity calculation method, one of:
#'        c("cosine", "pearson", "spearman", "euclidean") (default: "cosine")
#' @param up_weight Weight for positive markers when computing initial cell
#' scores (default: 0.5).
#' @param down_weight Weight for negative markers when computing initial cell
#' scores (default: 0.5).
#' @param seed_ratio Numeric in (0,1) or NULL. Quantile threshold to select seed
#' nodes from the CellScore distribution (e.g., 0.95). If NULL, the function
#' chooses a value adaptively  (default: NULL).
#' @param restart_prob umeric in [0,1] or NULL. Restart probability for RWR. If NULL,
#'   is chosen adaptively (default: NULL).
#' @param conv_threshold Numeric. Convergence threshold for RWR iterations (default: 1e-05)
#' @param nPermutations Integer. Number of permutations used to construct null distributions (default: 100)
#' @param nCores Integer. Number of CPU cores to use for permutation parallelization (default: 1).
#' @param SigThreshold Numeric in (0,1). Significance cutoff for empirical permutation p-values (default: 0.05).
#'
#' @return A list containing:
#' \itemize{
#'   \item pos_neg_score: Matrix with columns:
#'     \itemize{
#'       \item CellScore: Initial phenotype association score (before diffusion)
#'       \item pos_score: RWR diffusion score seeded from high CellScore cells
#'       \item neg_score: RWR diffusion score seeded from low CellScore cells
#'       \item CPlink_pos_score: Sigmoid-compressed positive diffusion score
#'       \item CPlink_neg_score: Sigmoid-compressed negative diffusion score
#'       \item CPlinkScore: Balanced CPlink phenotype score
#'     }
#'   \item res_permu: Permutation test results object
#' }
#'
#' @export
runCPlink <- function(Mtx,
                      DEGs_mtx,
                      snn,
                      MCA_dim = 50,
                      sim_methods = c("cosine", "pearson", "spearman", "euclidean")[1],
                      up_weight = 0.5,
                      down_weight = 0.5,
                      seed_ratio =NULL,
                      restart_prob = NULL,
                      conv_threshold = 1e-05,
                      nPermutations = 100,
                      nCores = 1,
                      SigThreshold = 0.05){

  if (is.null(DEGs_mtx$markers)) stop("`DEGs_mtx` must contain element `markers` (character vector of feature IDs).")
  if (nrow(snn) != ncol(snn)) stop("`snn` must be square (nodes x nodes).")
  if (!sim_methods %in% c("cosine","pearson","spearman","euclidean")) stop("`sim_methods` must be one of 'cosine','pearson','spearman','euclidean'.")

  message("Starting CPlink pipeline")

  # ---- MCA embedding restricted to phenotype-specific features ----
  message("Running MCA on phenotype-specific features")
    MCA_obj <- runMCA(Mtx = Mtx, markers = DEGs_mtx$markers)

    # ---- Compute feature-to-cell association matrix (features x cells) ----
    cell_feature_mtx <- GetCellMarkerDist(MCA_obj = MCA_obj, markers = DEGs_mtx$markers, ues_dim = MCA_dim,
                                      method = sim_methods)

    # ---- Compute initial per-cell/spot score using marker weights ----
    message("Calculating initial phenotype association scores")
    CellScore <- runCellScores(cell_feature_mtx,  DEGs_mtx, up_weight = up_weight, down_weight = down_weight, direction_col = "logFC")


    # ---- Adaptive seed_ratio & restart_prob (if requested) ----
    if(is.null(seed_ratio) || is.null(restart_prob)) {

      dens_score <- density(CellScore)

      find_maximal_points <- function(x) {

        m <- rle(x)
        maximal_points <- which(rep(diff(sign(diff(c(-Inf, m$values, -Inf)))) == -2,
                                    times = m$lengths))

        list(maximal_points = maximal_points, n = length(maximal_points))
      }

      maxpoints_info <- find_maximal_points(dens_score$y)
      maxpoints_x <- dens_score$x[maxpoints_info$maximal_points]
      maxpoints_y <- dens_score$y[maxpoints_info$maximal_points]

      pos_maxima <- maxpoints_x[maxpoints_x > max(CellScore)/4]
      neg_maxima <- maxpoints_x[maxpoints_x < min(CellScore)/4]
      mid_maxima <- maxpoints_x[(maxpoints_x <= max(CellScore)/4) & (maxpoints_x >= min(CellScore)/4)]

      pos_maxima_max <- ifelse(length(pos_maxima) > 0, max(pos_maxima), max(dens_score$x))
      neg_maxima_min <- ifelse(length(neg_maxima) > 0, min(neg_maxima), min(dens_score$x))
      mid_maxima_height <- ifelse(length(mid_maxima) > 0,
                                  max(maxpoints_y[maxpoints_x %in% mid_maxima]), 0)

      if(is.null(seed_ratio)) {

        if(pos_maxima_max > max(CellScore)/2 && neg_maxima_min < min(CellScore)/2) {

          seed_ratio <- 0.99

        } else {

          seed_ratio <- 0.95

        }
      }

      if(is.null(restart_prob)) {

        pos_maxima_height <- ifelse(length(pos_maxima) > 0, max(maxpoints_y[maxpoints_x %in% pos_maxima]), 0.00001)
        neg_maxima_height <- ifelse(length(neg_maxima) > 0, max(maxpoints_y[maxpoints_x %in% neg_maxima]), 0.00001)

        if(mid_maxima_height/pos_maxima_height > 50 || mid_maxima_height/neg_maxima_height > 50) {

          restart_prob <- 0.05

        } else {

          restart_prob <- 0.01

        }
      }
    }

    message("Restart Probability: ", restart_prob, "; Seeds Ratio: ", seed_ratio)

    # ---- select seeds based on quantiles of CellScore ----
    pos_seed_idx <- CellScore >= quantile(CellScore, seed_ratio)
    neg_seed_idx <- (-CellScore) >= quantile((-CellScore), seed_ratio)

    pos_queryCells <- rownames(snn)[pos_seed_idx]
    neg_queryCells <- rownames(snn)[neg_seed_idx]

    if (length(pos_queryCells) == 0) warning("No positive seed cells selected at the chosen seed_ratio.")
    if (length(neg_queryCells) == 0) warning("No negative seed cells selected at the chosen seed_ratio.")

    # ---- diffusion via Random Walk with Restart from positive & negative seeds ----
    pos_score <- runRWR(mtx = snn, seedcells = pos_queryCells,
                        restart_prob = restart_prob, conv_threshold = conv_threshold)
    neg_score <- runRWR(mtx = snn, seedcells = neg_queryCells,
                        restart_prob = restart_prob, conv_threshold = conv_threshold)

    # ---- sigmoid compression and final score normalization ----
    sigmoid_compress <- function(x) {

        1 / (1 + exp(-(x-mean(x))/sd(x)))

        }

    CPlink_pos_score <- sigmoid_compress(pos_score)
    CPlink_neg_score <- sigmoid_compress(neg_score)
    CPlinkScore <- (CPlink_pos_score-CPlink_neg_score)/max(abs(CPlink_pos_score-CPlink_neg_score))

    pos_neg_score <- cbind(CellScore, pos_score, neg_score, CPlink_pos_score, CPlink_neg_score, CPlinkScore)

    # ---- Permutation test to assign significance labels ----
    message("Computing empirical p-values and assign CPlink+ / CPlink- labels")
    res_permu <- runPermutationTest(nnSparseMtx = snn,
                               posSeedIdx = pos_seed_idx,
                               posSeedRWRScores = pos_score,
                               negSeedIdx = neg_seed_idx,
                               negSeedRWRScores = neg_score,
                               nPermutations=nPermutations,
                               SigThreshold=SigThreshold,
                               nCores=nCores,
                               restart_prob=restart_prob,
                               conv_threshold = conv_threshold)

    message("CPlink run completed successfully")
    return(list(pos_neg_score, res_permu))
}


#' CPlink-CODEX: Identify Phenotype-Associated Cells in Spatial Proteomics Data
#'
#' This function extends the CPlink framework to CODEX spatial proteomics
#' data by identifying phenotype-associated cells using transferred
#' phenotype association scores and cell-cell similarity networks.
#'
#' Unlike the standard CPlink pipeline, which directly computes initial
#' phenotype association scores from molecular features in target data,
#' CPlink-CODEX uses externally derived initial scores, typically transferred
#' from spatially aligned transcriptomics data.
#'
#' The function identifies positive and negative seed cells based on the
#' distribution of initial phenotype association scores. Phenotype signals
#' are subsequently propagated through the CODEX cell-cell graph using
#' Random Walk with Restart (RWR). Positive and negative diffusion scores
#' are combined into a balanced phenotype association score.
#'
#' Finally, a permutation-based statistical framework evaluates the
#' significance of diffusion scores and classifies CODEX cells into
#' CPlink+, CPlink-, and background populations.
#'
#' @param InitCellScore A numeric vector containing initial phenotype
#'   association scores for CODEX cells. These scores are typically
#'   obtained using getCodexInitCellScore(). The vector must correspond
#'   to the cell order in snn.
#' @param snn A square cell-cell similarity or adjacency matrix
#'   representing interactions among CODEX cells. Rows and columns
#'   correspond to CODEX cells and must follow the same ordering.
#'   Accepts dense or sparse matrices (e.g., dgCMatrix).
#' @param seed_ratio Numeric scalar in (0,1) or NULL. Quantile threshold
#'   for selecting positive and negative seed cells based on initial
#'   phenotype association scores. If NULL, the threshold is selected
#'   adaptively from the score distribution (default: NULL).
#' @param restart_prob Numeric scalar in [0,1] or NULL. Restart probability
#'   controlling the balance between network diffusion and return
#'   to seed nodes during RWR. If NULL, it is selected adaptively
#'   based on the score distribution (default: NULL).
#' @param conv_threshold Numeric scalar. Convergence threshold for
#'   iterative RWR diffusion (default: 1e-05).
#' @param nPermutations Integer scalar. Number of random permutations
#'   used to construct empirical null distributions for evaluating
#'   the statistical significance of diffusion scores (default: 100).
#' @param nCores Integer scalar. Number of CPU cores used for
#'   permutation parallelization (default: 1).
#' @param SigThreshold Numeric scalar in (0,1). Significance threshold
#'   applied to empirical permutation p-values for identifying
#'   phenotype-associated CODEX cells (default: 0.05).
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item pos_neg_score: A matrix containing per-cell phenotype
#'     association and diffusion scores, including:
#'     \itemize{
#'       \item CellScore: Initial phenotype association score.
#'       \item pos_score: RWR diffusion score from positive seed cells.
#'       \item neg_score: RWR diffusion score from negative seed cells.
#'       \item CPlink_pos_score: Sigmoid-compressed positive diffusion score.
#'       \item CPlink_neg_score: Sigmoid-compressed negative diffusion score.
#'       \item CPlinkScore: Balanced phenotype association score.
#'     }
#'   \item res_permu: A data.frame containing positive and negative
#'     seed indicators, permutation-based significance indicators,
#'     empirical p-values, and final cell labels.
#' }
#'
#' @seealso getCodexInitCellScore, mapSPTScoreToCODEX,
#'   runRWR, runPermutationTest
#'
#' @export
runCPlink_CODEX <- function(#Mtx,
  InitCellScore,
  snn,
  seed_ratio =NULL,
  restart_prob = NULL,
  conv_threshold = 1e-05,
  nPermutations = 100,
  nCores = 1,
  SigThreshold = 0.05){

  #if (is.null(DEGs_mtx$markers)) stop("`DEGs_mtx` must contain element `markers` (character vector of feature IDs).")
  if (nrow(snn) != ncol(snn)) stop("`snn` must be square (nodes x nodes).")
  # if (!sim_methods %in% c("cosine","pearson","spearman","euclidean")) stop("`sim_methods` must be one of 'cosine','pearson','spearman','euclidean'.")

  #  message("Starting CPlink pipeline")

  # ---- MCA embedding restricted to phenotype-specific features ----
  #  message("Running MCA on phenotype-specific features")
  #   MCA_obj <- runMCA(Mtx = Mtx, markers = DEGs_mtx$markers)

  # ---- Compute feature-to-cell association matrix (features x cells) ----
  #   cell_feature_mtx <- GetCellMarkerDist(MCA_obj = MCA_obj, markers = DEGs_mtx$markers, ues_dim = MCA_dim,
  #                                     method = sim_methods)

  # ---- Compute initial per-cell/spot score using marker weights ----
  #    message("Calculating initial phenotype association scores")
  CellScore <- InitCellScore


  # ---- Adaptive seed_ratio & restart_prob (if requested) ----
  if(is.null(seed_ratio) || is.null(restart_prob)) {

    dens_score <- density(CellScore)

    find_maximal_points <- function(x) {

      m <- rle(x)
      maximal_points <- which(rep(diff(sign(diff(c(-Inf, m$values, -Inf)))) == -2,
                                  times = m$lengths))

      list(maximal_points = maximal_points, n = length(maximal_points))
    }

    maxpoints_info <- find_maximal_points(dens_score$y)
    maxpoints_x <- dens_score$x[maxpoints_info$maximal_points]
    maxpoints_y <- dens_score$y[maxpoints_info$maximal_points]

    pos_maxima <- maxpoints_x[maxpoints_x > max(CellScore)/4]
    neg_maxima <- maxpoints_x[maxpoints_x < min(CellScore)/4]
    mid_maxima <- maxpoints_x[(maxpoints_x <= max(CellScore)/4) & (maxpoints_x >= min(CellScore)/4)]

    pos_maxima_max <- ifelse(length(pos_maxima) > 0, max(pos_maxima), max(dens_score$x))
    neg_maxima_min <- ifelse(length(neg_maxima) > 0, min(neg_maxima), min(dens_score$x))
    mid_maxima_height <- ifelse(length(mid_maxima) > 0,
                                max(maxpoints_y[maxpoints_x %in% mid_maxima]), 0)

    if(is.null(seed_ratio)) {

      if(pos_maxima_max > max(CellScore)/2 && neg_maxima_min < min(CellScore)/2) {

        seed_ratio <- 0.99

      } else {

        seed_ratio <- 0.95

      }
    }

    if(is.null(restart_prob)) {

      pos_maxima_height <- ifelse(length(pos_maxima) > 0, max(maxpoints_y[maxpoints_x %in% pos_maxima]), 0.00001)
      neg_maxima_height <- ifelse(length(neg_maxima) > 0, max(maxpoints_y[maxpoints_x %in% neg_maxima]), 0.00001)

      if(mid_maxima_height/pos_maxima_height > 50 || mid_maxima_height/neg_maxima_height > 50) {

        restart_prob <- 0.05

      } else {

        restart_prob <- 0.01

      }
    }
  }

  message("Restart Probability: ", restart_prob, "; Seeds Ratio: ", seed_ratio)

  # ---- select seeds based on quantiles of CellScore ----
  pos_seed_idx <- CellScore >= quantile(CellScore, seed_ratio)
  neg_seed_idx <- (-CellScore) >= quantile((-CellScore), seed_ratio)

  pos_queryCells <- rownames(snn)[pos_seed_idx]
  neg_queryCells <- rownames(snn)[neg_seed_idx]

  if (length(pos_queryCells) == 0) warning("No positive seed cells selected at the chosen seed_ratio.")
  if (length(neg_queryCells) == 0) warning("No negative seed cells selected at the chosen seed_ratio.")

  # ---- diffusion via Random Walk with Restart from positive & negative seeds ----
  pos_score <- runRWR(mtx = snn, seedcells = pos_queryCells,
                      restart_prob = restart_prob, conv_threshold = conv_threshold)
  neg_score <- runRWR(mtx = snn, seedcells = neg_queryCells,
                      restart_prob = restart_prob, conv_threshold = conv_threshold)

  # ---- sigmoid compression and final score normalization ----
  sigmoid_compress <- function(x) {

    1 / (1 + exp(-(x-mean(x))/sd(x)))

  }

  CPlink_pos_score <- sigmoid_compress(pos_score)
  CPlink_neg_score <- sigmoid_compress(neg_score)
  CPlinkScore <- (CPlink_pos_score-CPlink_neg_score)/max(abs(CPlink_pos_score-CPlink_neg_score))

  pos_neg_score <- cbind(CellScore, pos_score, neg_score, CPlink_pos_score, CPlink_neg_score, CPlinkScore)

  # ---- Permutation test to assign significance labels ----
  message("Computing empirical p-values and assign CPlink+ / CPlink- labels")
  res_permu <- runPermutationTest(nnSparseMtx = snn,
                                  posSeedIdx = pos_seed_idx,
                                  posSeedRWRScores = pos_score,
                                  negSeedIdx = neg_seed_idx,
                                  negSeedRWRScores = neg_score,
                                  nPermutations=nPermutations,
                                  SigThreshold=SigThreshold,
                                  nCores=nCores,
                                  restart_prob=restart_prob,
                                  conv_threshold = conv_threshold)

  message("CPlink run completed successfully")
  return(list(pos_neg_score, res_permu))
}
