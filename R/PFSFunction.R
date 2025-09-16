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
                   design_formula = "Surv(time, status) ~",
                   time_event_name = c("time", "status"),
                   cv_quantile = 0.05,
                   p.adj = 0.05,
                   n_cores = parallel::detectCores()){

  if (!all(colnames(bulk_count) %in% rownames(bulk_meta))) {
    stop("Column names of bulk_count must be present in rownames(bulk_meta)")
  }

  registerDoParallel(cores = n_cores)
  dds <- DESeqDataSetFromMatrix(bulk_count, colData = bulk_meta, design = ~1)
  dds <- estimateSizeFactors(dds)
  bulk_norm <- t(counts(dds, normalized=TRUE))
  gene_cv <- foreach(g=1:nrow(bulk_norm), .combine=c) %dopar% {
    x <- as.numeric(bulk_norm[g, ])
    sd(x)/mean(x)
  }
  bulk_norm_filter <- bulk_norm[, gene_cv > quantile(gene_cv, 0.05)]
  surv_df <- data.frame(
    bulk_meta[rownames(bulk_norm_filter), time_event_name],
    bulk_norm_filter
  )
  genes <- colnames(surv_df)[3:ncol(surv_df)]
  cox_results <- foreach(gene = genes, .combine = rbind) %dopar% {
    formula <- as.formula(paste(design_formula, gene))
    cox_fit <- coxph(formula, data = surv_df)
    cox_sum <- summary(cox_fit)
    data.frame(
      Gene = gene,
      HR = cox_sum$coef[1, "exp(coef)"],
      logHR = log2(cox_sum$coef[1, "exp(coef)"]),

      Pvalue = cox_sum$coef[1, "Pr(>|z|)"]
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


#' Identify Phenotype-specific Features with edgeR
#'
#' Similar to DESeq2, uses edgeR for differential analysis on bulk counts for binary phenotypes.
#' Features significantly associated with the target phenotype are selected as markers, and
#' log2 fold-change values are used as importance weights for downstream CPlink analysis.
#'
#' @param bulk.count A count matrix with features as rows and samples as columns
#' @param group A vector or factor specifying the experimental groups
#' @param design_formula A formula specifying the experimental design (default: ~ed_phenotype)
#' @param p.adj Adjusted p-value cutoff for significance (default: 0.05)
#' @param logFC Minimum absolute log2 fold-change threshold (default: 0)
#'
#' @return A list containing:
#' \itemize{
#'   \item Full edgeR results table (all features)
#'   \item Filtered results containing only significant features meeting p.adj and logFC thresholds
#' }
#'
#'
runedgeR <- function(bulk.count, group,
                     design_formula = ~ed_phenotype,
                     p.adj = 0.05,
                     logFC = 0){

  if (length(group) != ncol(bulk.count)) stop("group must have length equal to ncol(bulk.count)")

  dge <- DGEList(counts = bulk.count, group = group)
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge)
  fit <- glmFit(dge, design = model.matrix(design_formula))
  lrt <- glmLRT(fit)
  edgeR_res <- topTags(lrt, n = Inf)$table
  edgeR_res_sig <- subset(edgeR_res, FDR < p.adj)
  edgeR_sig_markers <- edgeR_res_sig[which(abs(edgeR_res_sig$logFC) > logFC),]

  return(list(edgeR_res, edgeR_sig_markers))
}


#' Identify Phenotype-specific Features with limma-voom
#'
#' Performs differential analysis using limma-voom for bulk RNA-seq data with binary phenotypes.
#' Top features significantly associated with the phenotype are selected and assigned log2 fold-change
#' weights.
#'
#' @param bulk.count A count matrix with features as rows and samples as columns
#' @param conditions A factor vector specifying the experimental conditions
#' @param design_formula A formula specifying the experimental design (default: ~conditions)
#' @param p.adj Adjusted p-value cutoff for significance (default: 0.05)
#' @param logFC Minimum absolute log2 fold-change threshold (default: 0)
#'
#' @return A list containing:
#' \itemize{
#'   \item Full limma results table (all features)
#'   \item Filtered results containing only significant features meeting p.adj and logFC thresholds
#' }
#'
#'
runlimma <- function(bulk.count,
                     conditions,
                     design_formula = ~conditions,
                     p.adj = 0.05,
                     logFC = 0){

  if (length(conditions) != ncol(bulk.count)) stop("conditions must match ncol(bulk.count)")

  design <- model.matrix(design_formula)
  colnames(design) <- levels(conditions)
  rownames(design) <- colnames(bulk.count)
  fit1 <- lmFit(bulk.count, design)
  efit <- eBayes(fit1)
  limma_res <- topTable(efit, coef=ncol(design), n=nrow(bulk.count), adjust.method = "fdr")
  limma_res_sig <- subset(limma_res, P.Value < p.adj)
  limma_sig_markers <- limma_res_sig[which(abs(limma_res_sig$logFC) > logFC),]

  return(list(limma_res, limma_sig_markers))
}

