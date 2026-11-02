#' Identify prognosis-associated Features with Cox Regression
#'
#' For survival phenotypes, performs univariate Cox proportional hazards regression on normalized bulk data.
#' Features with lowest 5% coefficient of variation are excluded to reduce noise. Log2 hazard ratios (HR)
#' are used as phenotype-specific weights
#'
#' @param bulk_count A count matrix with features as rows and samples as columns
#' @param bulk_meta A data frame containing sample metadata with matching row names to count columns
#' @param design_formula design_formula A formula specifying the experimental design (default: "Surv(time, status) ~")
#' @param time_event_name the column name for survival time and event status
#' @param cv_quantile Filter out a certain proportion of low-variance features by using the coefficient of variation
#' @param p.adj Adjusted p-value cutoff for significance (default: 0.05)
#' @param norm Logical. Whether to perform DESeq2 size-factor normalization on bulk count data (default: FALSE).
#' @param min_prop Numeric scalar in [0,1]. Minimum proportion of samples with feature expression greater than 1 required for retaining a feature (default: 0.05).
#' @param p_cox Numeric scalar. Nominal p-value threshold for Cox regression (default: 0.05). Currently not used in the function.
#' @param n_cores Number of cores for parallel processing. default: parallel::detectCores() (all available cores)
#'
#' @return A list containing:
#' \itemize{
#'   \item Cox results for all features
#'   \item Cox results for significant genes only
#' }
#'
#' @export
runCox <- function(bulk_count,
                   bulk_meta,
                   design_formula = "1",
                   time_event_name = c("time", "status"),
                   norm = FALSE,
                   min_prop = 0.05,
                   p.adj = 0.05,p_cox=0.05,
                   n_cores = parallel::detectCores()){

  if (!all(colnames(bulk_count) %in% rownames(bulk_meta))) {
    stop("Column names of bulk_count must be present in rownames(bulk_meta)")
  }
  bulk_meta <- bulk_meta[colnames(bulk_count), , drop = FALSE]

  min_prop <- min_prop
  bulk_count <- bulk_count[rowMeans(bulk_count > 1) >= min_prop, ]
  registerDoParallel(cores = n_cores)
  if ( norm == TRUE){
    dds <- DESeqDataSetFromMatrix(bulk_count, colData = bulk_meta, design = ~1)
    dds <- estimateSizeFactors(dds)

    bulk_norm <- t(counts(dds, normalized=TRUE))
  } else {
    bulk_norm <- t(bulk_count)
  }

  bulk_norm_mat <- as.matrix(bulk_norm)

  genes <- colnames(bulk_norm_mat)
  cox_results <- foreach(gene = genes, .combine = rbind) %dopar% {


    dat <- data.frame(
      bulk_meta[rownames(bulk_norm_mat), time_event_name],
      gene_expr = bulk_norm_mat[, gene]
    )
    dat$gene_expr <- as.numeric(scale(dat$gene_expr))

    rhs <- if (design_formula == "1") {
      "gene_expr"
    } else {
      paste(design_formula, "gene_expr", sep = " + ")
    }

    formula_text <- paste0(
      "Surv(", time_event_name[1], ", ", time_event_name[2], ") ~ ", rhs
    )

    formula <- as.formula(formula_text)


    cox_fit <- coxph(formula, data = dat)
    cox_sum <- summary(cox_fit)


    cat(gene)
    data.frame(
      Gene = gene,
      beta = cox_sum$coef["gene_expr", "coef"],
      HR = cox_sum$coef["gene_expr", "exp(coef)"],
      logHR = log2(cox_sum$coef["gene_expr",  "exp(coef)"]),

      Pvalue = cox_sum$coef["gene_expr",  "Pr(>|z|)"]
    )
  }

  rownames(cox_results) <- cox_results$Gene
  cox_results$FDR <- p.adjust(cox_results$Pvalue, method = "fdr")
  significant_genes <- subset(cox_results, FDR < p.adj)
  significant_genes <- significant_genes[order(significant_genes$FDR), ]

  return(list(cox_results, significant_genes))

}


#' Identify Phenotype-specific Features with DESeq2
#'
#' For binary group phenotypes (e.g., disease vs. control), this function performs differential
#' feature analysis on bulk count data using DESeq2. Features significantly associated
#' with the phenotype are ranked, and top features are selected as potential phenotype-specific markers.
#' The log2 fold-change values are used as phenotype-specific importance weights.
#'
#' @param bulk.count A count matrix with features as rows and samples as columns
#' @param bulk.meta A data frame containing sample metadata with matching row names to count columns
#' @param design_formula A formula specifying the experimental design (default: ~ phenotype)
#' @param pAdjustMethod Method for p-value adjustment ("fdr", "BH", "bonferroni", etc.)
#' @param p.adj Adjusted p-value cutoff for significance (default: 0.05)
#' @param logFC Minimum absolute log2 fold-change threshold (default: 0)
#' @param parallel Logical indicating whether to use parallel processing (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item Full DESeq2 results object (all features)
#'   \item Filtered results containing only significant features meeting p.adj and logFC thresholds
#' }
#'
#' @export
runDESeq2 <- function(bulk.count,
                      bulk.meta,
                      design_formula = ~ phenotype,
                      pAdjustMethod = "fdr",
                      p.adj = 0.05,
                      logFC = 0,
                      parallel = TRUE){

  if (!all(colnames(bulk.count) %in% rownames(bulk.meta))) {
    stop("Column names of bulk.count must be present as rownames of bulk.meta")
  }

  dds <- DESeqDataSetFromMatrix(countData = as.matrix(bulk.count),
                                colData = as.data.frame(bulk.meta),
                                design = design_formula)
  dds <- DESeq(dds, parallel = parallel)
  res <- results(dds, pAdjustMethod=pAdjustMethod)

  sig_markers <- res[which(res$padj < p.adj), ]
  sig_markers <- sig_markers[which(abs(sig_markers$log2FoldChange) > logFC),]

  return(list(res, sig_markers))
}



#' Identify Phenotype-Associated Proteins Using DEP2
#'
#' This function performs differential protein expression analysis on
#' bulk proteomics data using the DEP2 framework. It identifies proteins
#' associated with a phenotype of interest by comparing protein abundance
#' profiles between experimental conditions.
#'
#' The workflow includes construction of a SummarizedExperiment object,
#' optional missing-value filtering, variance stabilizing normalization
#' (VSN), missing-value imputation, and differential protein expression
#' analysis using linear models and empirical Bayes statistics.
#'
#' Significant proteins are identified based on statistical significance
#' and absolute log2 fold-change thresholds. The resulting differential
#' protein profiles can be used to derive phenotype-associated molecular
#' features and their corresponding weights for downstream CPlink analysis.
#'
#' @param bulk.data A numeric protein abundance matrix with proteins
#'   as rows and samples as columns. Row names represent protein
#'   identifiers and column names represent sample identifiers.
#'   Missing abundance values should be represented as NA.
#' @param phenotype A character vector or factor specifying the
#'   experimental condition of each sample. Its length and order
#'   must correspond to the columns of bulk.data.
#' @param replicate A vector specifying replicate identifiers for
#'   individual samples. Its length and order must correspond
#'   to the columns of bulk.data.
#' @param filter Logical. Whether to perform missing-value filtering
#'   before normalization (default: TRUE).
#' @param thr Integer scalar. Maximum number of missing values
#'   allowed in at least one experimental condition for a protein
#'   to pass the corresponding filter (default: 1).
#' @param fraction Numeric scalar. Maximum proportion of missing
#'   values allowed per protein when missing-value filtering is
#'   enabled (default: 0.5).
#' @param normalization Logical. Whether to apply variance
#'   stabilizing normalization using VSN (default: TRUE).
#' @param imputation Logical. Whether to impute missing protein
#'   abundance values (default: TRUE).
#' @param imp_fun Character string specifying the missing-value
#'   imputation method supported by DEP2 (default: "MinProb").
#' @param imp_q Numeric scalar. Quantile parameter used by the
#'   MinProb imputation method (default: 0.01).
#' @param ctrl Character string specifying the reference or
#'   control condition (default: "normal").
#' @param treat Character string specifying the experimental
#'   condition of interest (default: "glioblastoma").
#' @param p_adj Logical. Whether to use adjusted p-values
#'   instead of raw p-values for selecting significant proteins
#'   (default: TRUE).
#' @param p_val Numeric scalar. Significance threshold for
#'   differential protein expression analysis (default: 0.05).
#' @param lfc Numeric scalar. Minimum absolute log2 fold-change
#'   required to select significant proteins (default: 0).
#'
#' @return A named list containing:
#' \itemize{
#'   \item full: A data.frame containing differential protein
#'     expression results for all retained proteins, including
#'     protein identifiers, log2 fold changes, raw p-values,
#'     adjusted p-values, and significance annotations.
#'   \item sig: A data.frame containing proteins satisfying the
#'     specified statistical significance and absolute log2
#'     fold-change thresholds for the treat versus ctrl comparison.
#' }
#'
#' @seealso DEP2::make_unique, DEP2::make_se, DEP2::filter_se,
#'   DEP2::normalize_vsn, DEP2::impute, DEP2::test_diff,
#'   DEP2::add_rejections, DEP2::get_results
#'
#' @export
run_DEP2 <- function(bulk.data,
                     phenotype,
                     replicate,
                     filter = TRUE,
                     thr = 1,
                     fraction = 0.5,
                     normalization = TRUE,
                     imputation = TRUE,
                     imp_fun = "MinProb",
                     imp_q = 0.01,
                     ctrl = "normal",
                     treat = "glioblastoma",
                     p_adj = FALSE,
                     p_val = 0.05,
                     lfc = 0) {

  bulk_design <- data.frame(label = colnames(bulk.data),  condition = phenotype, replicate = replicate)

  bulk_exp <- data.frame(bulk.data, name = rownames(bulk.data), ID = rownames(bulk.data))

  bulk_exp <- make_unique(bulk_exp, "name", "ID", delim = ";")

  num_cols <- which(sapply(bulk_exp, is.numeric))
  se <- make_se(bulk_exp, num_cols, bulk_design)

  if(filter){
    se <- filter_se(se,  thr = thr,
                    fraction = fraction
    )
  }

  if(normalization){
    se <- normalize_vsn(se)
  }

  if(imputation){
    se <- impute(se, fun = imp_fun, q = imp_q)
  }

  diff_res <- test_diff(se, type = "control", control = ctrl, fdr.type = "BH")

  dep_res <- DEP2::add_rejections(
    diff_res,
    alpha = p_val,
    lfc = lfc
  )

  res_table <- DEP2::get_results(dep_res)

  fc_col <- paste0(treat, "_vs_", ctrl, "_ratio")
  if(p_adj){
    p_col <- paste0(treat, "_vs_", ctrl, "_p.adj")
    sig <- res_table[which(res_table[[p_col]] < p_val & abs(res_table[[fc_col]]) > lfc), ]
  } else{
    p_col <- paste0(treat, "_vs_", ctrl, "_p.val")
    sig <- res_table[which(res_table[[p_col]] < p_val & abs(res_table[[fc_col]]) > lfc), ]
  }

  return(list(full = res_table, sig = sig))

}
