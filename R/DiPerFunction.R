
#' Random Walk with Restart (RWR) diffusion on a cell/spot graph
#'
#' Implements a random-walk-with-restart diffusion process to propagate
#' seed signals over a cell/spot similarity/adjacency matrix. This routine
#' follows the formulation commonly used in CPlink: starting from a
#' probability vector concentrated on seed nodes, the score is iteratively
#' diffused along edges with periodic restarting to the seed distribution.
#'
#' @param mtx Numeric matrix representing the adjacency / transition counts
#'        between nodes. Rows/columns correspond to the same set of nodes.
#'        The function will convert this to a column-stochastic transition
#'        matrix internally.
#' @param seedcells Character vector or integer indices specifying seed nodes.
#' @param restart_prob Restart probability (scalar in [0,1]). Typical
#'        values (e.g. 0.05) preserve local propagation while allowing return
#'        to seeds.
#' @param conv_threshold Convergence threshold on L1 change between iterations.
#' @param max_iter Maximum number of iterations to run (guards against non-convergence).
#'
#' @return A named numeric vector of diffusion scores for each node.
#'
#' @export
runRWR <- function(mtx,
                   seedcells,
                   restart_prob=0.05,
                   conv_threshold=1e-5){

  mtx <- t(mtx) / colSums(mtx)

  prob_vector <- numeric(ncol(mtx))
  names(prob_vector) <- colnames(mtx)
  prob_vector[seedcells] <- 1
  prob_vector <- t(prob_vector/sum(prob_vector))

  RWR <- function(t_mtx, p_init, p_rest, cutoff) {
    Wt <- t(t_mtx)
    p_curr <- p_init
    iter <- 0
    delta <- 1

    while  (delta > cutoff) {
      p_next <- t(((1-p_rest)*Wt) %*% t(p_curr)) + (p_rest*p_init)
      delta <- sum(abs(p_next - p_curr))
      p_curr <- p_next
      iter <- iter + 1
    }

    #    message("Step: ", iter, "; Delta: ", delta)
    return(drop(p_curr))
  }

  output <- RWR(t_mtx = mtx, p_init = prob_vector, p_rest = restart_prob, cutoff = conv_threshold)

  return(output)
}


#' Permutation Test for Identifying Phenotype-Associated Cells or spots (CPlink)
#'
#' This function implements the permutation-based statistical framework of CPlink
#' to assess the significance of Random Walk with Restart (RWR) diffusion scores.
#' It evaluates whether cells/spots are significantly associated with positive
#' or negative seed sets compared to a null distribution constructed by random seeds.
#'
#' @param nnSparseMtx Sparse transition matrix
#' (cells x cells or spots x spots) used for RWR.
#' @param posSeedIdx Integer indices of positive seed nodes.
#' @param posSeedRWRScores Numeric vector of diffusion scores from positive seeds.
#' @param negSeedIdx Integer indices of negative seed nodes.
#' @param negSeedRWRScores Numeric vector of diffusion scores from negative seeds.
#' @param nPermutations Number of random permutations to generate the null distribution (default: 100).
#' @param SigThreshold Significance threshold for adjusted p-values (default: 0.05).
#' @param nCores Number of parallel cores for permutation (default: 1).
#' @param restart_prob Restart probability for RWR (default: 0.05).
#' @param conv_threshold Convergence threshold for RWR iterations (default: 1e-5).
#'
#' @return A data frame with per-cell or per-spot results.
#'
#'
runPermutationTest <- function(nnSparseMtx,
                               posSeedIdx,
                               posSeedRWRScores,
                               negSeedIdx,
                               negSeedRWRScores,
                               nPermutations=100,
                               SigThreshold=0.05,
                               nCores=1,
                               restart_prob=0.05,
                               conv_threshold = 1e-05){

  # ----- Generate null distribution by random seeds -----
  permutationScores <- mclapply(seq_len(nPermutations), mc.cores = nCores, function(i){

    sampleIdx <- sample(1:nrow(nnSparseMtx), round((nrow(nnSparseMtx[posSeedIdx,]) + nrow(nnSparseMtx[negSeedIdx,]))/2))
    sampleScores <- runRWR(mtx=nnSparseMtx, seedcells=rownames(nnSparseMtx)[as.numeric(sampleIdx)],
                           restart_prob=restart_prob, conv_threshold = conv_threshold)

    return(sampleScores)
  })
  names(permutationScores) <- paste0("permutation_", seq_len(nPermutations))
  permutationScores <- data.frame(sapply(permutationScores, c))

  posPermutationLogical <- data.frame(matrix(nrow=nrow(nnSparseMtx), ncol=nPermutations))
  posPermutationLogical <- apply(permutationScores, 2, function(x) {

    temp <- x > posSeedRWRScores

    return(temp)
  })

  posCellsPval <- (1+rowSums(posPermutationLogical)) / (1+nPermutations)
  names(posCellsPval) <- rownames(nnSparseMtx)
  posTrueCellIdx <- posCellsPval <= SigThreshold


  negPermutationLogical <- data.frame(matrix(nrow=nrow(nnSparseMtx), ncol=nPermutations))
  negPermutationLogical <- apply(permutationScores, 2, function(x) {

    temp <- x > negSeedRWRScores

    return(temp)
  })


  negCellsPval <- (1+rowSums(negPermutationLogical)) / (1+nPermutations)
  names(negCellsPval) <- rownames(nnSparseMtx)
  negTrueCellIdx <- negCellsPval <= SigThreshold

  overlapCells <- intersect(which(posTrueCellIdx), which(negTrueCellIdx))
  if (length(overlapCells) > 0) {

    for (cell in overlapCells) {

      if (posCellsPval[cell] < negCellsPval[cell]) {

        negTrueCellIdx[cell] <- FALSE

      } else {

        posTrueCellIdx[cell] <- FALSE

      }
    }
  }

  cellLabels <- rep("Background", nrow(nnSparseMtx))
  cellLabels[posTrueCellIdx] <- "CPlink+ cells"
  cellLabels[negTrueCellIdx] <- "CPlink- cells"

  TureCellDF <- data.frame(posSeedIdx, negSeedIdx, posTrueCellIdx, negTrueCellIdx, cellLabels, posCellsPval, negCellsPval)

  return(TureCellDF)
}

