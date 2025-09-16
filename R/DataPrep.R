#' Run standard scRNA-seq preprocessing and visualization
#'
#' This function performs an end-to-end single-cell RNA-seq analysis workflow, including
#' quality control, normalization, highly variable feature selection, dimensionality reduction,
#' optional batch effect correction using Harmony, clustering, and UMAP-based visualization.
#' The resulting Seurat object can be used for downstream analyses.
#'
#' @param counts A raw or filtered count matrix with genes as rows and cells as columns.
#' @param meta.data A data frame containing cell-level metadata. Row names must match the column names of `counts`.
#' @param batch_correct Logical, whether to perform batch effect correction (default: TRUE).
#' @param batch_var Character string specifying the metadata column containing batch information (default: "patient").
#' @param min.cells Minimum number of cells in a feature for it to be retained (default: 3).
#' @param nfeatures Number of highly variable features to select (default: 2000).
#' @param resolution Cluster resolution parameter (default: 0.8).
#' @param dims_Neighbors Dimensions to use for FindNeighbors (default: 1:50).
#' @param dims_UMAP Dimensions to use for UMAP projection (default: 1:50).
#'
#' @return A Seurat object containing all analysis results.
#'
#' @export
runRNA <- function(counts,
                   meta.data,
                   batch_correct = TRUE,
                   batch_var = "batch",
                   min.cells = 3,
                   nfeatures = 2000,
                   resolution = 0.8,
                   dims_Neighbors = 1:50,
                   dims_UMAP = 1:50) {

  # basic checks
  if (is.null(colnames(counts)) || is.null(rownames(meta.data))) {
    stop("counts must have column names and meta.data must have row names matching column names of counts")
  }
  if (!all(colnames(counts) %in% rownames(meta.data))) {
    stop("Not all cell names in counts are present in rownames(meta.data)")
  }

  # ---- Seurat object creation and preprocessing ----
  sc_seurat <- CreateSeuratObject(counts = counts,  min.cells = min.cells, meta.data = meta.data)
  sc_seurat <- NormalizeData(sc_seurat, verbose = FALSE)
  sc_seurat <- FindVariableFeatures(sc_seurat, nfeatures = nfeatures, verbose = FALSE)
  sc_seurat <- ScaleData(sc_seurat, verbose = FALSE)
  sc_seurat <- RunPCA(sc_seurat, verbose = FALSE)

  # ---- Batch correction with Harmony if requested ----
  if (batch_correct) {
    if (!batch_var %in% colnames(meta.data)) {
      stop(paste("Batch variable", batch_var, "not found in metadata"))
    }
    sc_seurat <- RunHarmony(sc_seurat, group.by.vars = batch_var, verbose = FALSE)
    reduction_used <- "harmony"
  } else {
    reduction_used <- "pca"
  }

  # ---- Graph construction, clustering, and UMAP embedding ----
  sc_seurat <- FindNeighbors(sc_seurat, reduction = reduction_used, dims = dims_Neighbors, verbose = FALSE)
  sc_seurat <- FindClusters(sc_seurat, resolution = resolution, verbose = FALSE)
  sc_seurat <- RunUMAP(sc_seurat, reduction = reduction_used, dims = dims_UMAP, verbose = FALSE)

  return(sc_seurat)
}


#' Run standard scATAC-seq preprocessing and visualization
#'
#' This function performs a complete single-cell ATAC-seq analysis workflow, including
#' chromatin assay construction, TF-IDF normalization, highly informative feature selection,
#' Latent Semantic Indexing (LSI) dimensionality reduction, optional batch effect correction via Harmony,
#' graph-based clustering, and UMAP visualization. The resulting Seurat object can be used for
#' downstream analyses.
#'
#' @param matrix A count matrix with peaks/bins as rows and cells as columns.
#' @param meta.data A data frame containing cell-level metadata. Row names must match the column names of `matrix`.
#' @param ranges A GRanges object specifying genomic ranges for peaks. If NULL, ranges will be inferred from `matrix` row names.
#' @param genome Genome identifier (e.g., "hg38"). Required if ranges not provided.
#' @param sep Separator characters in peak names (default c("-", "-") for "chr-start-end" format)
#' @param batch_correct Logical indicating whether to perform Harmony batch correction (default: FALSE)
#' @param batch_var Metadata column name containing batch information (default: "tech")
#' @param min.cells Minimum cells required to keep a peak (default: 3)
#' @param min.cutoff Minimum cutoff for FindTopFeatures ("q0" uses all features, see ?FindTopFeatures)
#' @param resolution Cluster resolution parameter (default: 0.8)
#' @param dims_Neighbors Dimensions for FindNeighbors (default: 2:50)
#' @param dims_UMAP Dimensions for UMAP projection (default: 2:50)
#' @param fragments Path to fragments file (optional, for fragment counting)
#'
#' @return A Seurat object containing all analysis results.
#'
#' @export

runATAC <- function(matrix,
                    meta.data,
                    ranges = NULL,
                    genome = NULL,
                    sep = c("-", "-"),
                    batch_correct = FALSE,
                    batch_var = "batch",
                    min.cells = 3,
                    min.cutoff = "q0",
                    resolution = 0.8,
                    dims_Neighbors = 2:50,
                    dims_UMAP = 2:50,
                    fragments = NULL) {

  # input checks
  if (is.null(colnames(matrix)) || is.null(rownames(meta.data))) {
    stop("matrix must have column names and meta.data must have rownames matching cell names")
  }
  if (!all(colnames(matrix) %in% rownames(meta.data))) {
    stop("Not all cell names in matrix are present in rownames(meta.data)")
  }
  # try to infer ranges if NULL
  if (is.null(ranges) && is.null(genome)) {
    warning("ranges not supplied. If peak names are not in chr-start-end format, please provide 'ranges' or 'genome'")
  }

  # ---- Create Seurat object with ChromatinAssay ----
  sc_seurat <- CreateSeuratObject(counts =
                                    CreateChromatinAssay(counts = matrix,
                                                         ranges = ranges,
                                                         genome = genome,
                                                         sep = sep,
                                                         fragments = fragments,
                                                         min.cells = min.cells,
                                                         min.features = 0),
                                  assay = "ATAC", meta.data = meta.data
  )

  # ---- TF-IDF normalization and feature selection ----
  sc_seurat <- RunTFIDF(sc_seurat, verbose = FALSE)
  sc_seurat <- FindTopFeatures(sc_seurat, min.cutoff = min.cutoff, verbose = FALSE)
  sc_seurat <- RunSVD(sc_seurat, verbose = FALSE)

  # ---- Optional batch correction using Harmony ----
  if (batch_correct) {
    if (!batch_var %in% colnames(meta.data)) {
      stop(paste("Batch variable", batch_var, "not found in metadata"))
    }
    sc_seurat <- RunHarmony(
      object = sc_seurat,
      group.by.vars = batch_var,
      reduction = 'lsi',
      assay.use = 'peaks',
      project.dim = FALSE,
      verbose = FALSE
    )
    reduction_used <- "harmony"
  } else {
    reduction_used <- "lsi"
  }

  # ---- Graph construction, clustering, and UMAP visualization ----
  sc_seurat <- FindNeighbors(object = sc_seurat, reduction = reduction_used, dims = dims_Neighbors, verbose = FALSE)
  sc_seurat <- FindClusters(object = sc_seurat, algorithm = 3, resolution = resolution, verbose = FALSE)
  sc_seurat <- RunUMAP(object = sc_seurat, reduction = reduction_used, dims = dims_UMAP, verbose = FALSE)

  return(sc_seurat)
}


#' Run standard ST preprocessing and visualization
#'
#' This function provides a complete workflow for processing spatial transcriptomics data.
#' It reads count matrices and spatial image information, constructs a Seurat object,
#' performs normalization using SCTransform, reduces dimensionality, identifies clusters,
#' and generates UMAP embeddings for visualization. The output Seurat object integrates both
#' gene expression and spatial image information for downstream analyses.
#'
#' @param data_dir Path to the base directory containing the spatial transcriptomics dataset.
#' @param image_name Name of the spatial image file (default: "tissue_lowres_image.png").
#' @param min.cells Minimum number of cells required to keep a feature (default: 3).
#' @param min.features Minimum number of features required to keep a cell (default: 0).
#' @param dims Dimensions to use for PCA/Neighbors/UMAP (default: 1:50).
#' @param resolution Resolution parameter for clustering (default: 0.8).
#' @param verbose Logical; whether to print progress messages (default: FALSE).
#'
#' @return A Seurat object with spatial image and analysis results.
#'
#' @export
runSpatial_10X <- function(data_dir,
                       image_name = "tissue_lowres_image.png",
                       min.cells = 3,
                       min.features = 0,
                       dims_Neighbors = 1:50,
                       dims_UMAP = 1:50,
                       resolution = 0.8,
                       verbose = FALSE) {

  # ---- Read spatial transcriptomics counts ----
  counts <- Read10X(file.path(data_dir, "filtered_feature_bc_matrix"))

  # ---- Create Seurat object with spatial assay ----
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
    assay = "Spatial"
  )

  # ---- Load and attach spatial image ----
  spatial_image <- Read10X_Image(
    image.dir = file.path(data_dir, "spatial"),
    image.name = image_name,
    filter.matrix = TRUE
  )
  spatial_image <- spatial_image[Cells(x = seurat_obj)]
  DefaultAssay(object = spatial_image) <- "Spatial"
  seurat_obj[["slice1"]] <- spatial_image

  # ---- Preprocessing: normalization, dimensionality reduction, Graph construction, clustering, and UMAP visualization ----
  seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = verbose)
  seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = verbose)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = dims_Neighbors, verbose = verbose)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = verbose)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = dims_UMAP, verbose = verbose)

  return(seurat_obj)
}


#' Run standard ST preprocessing and visualization
#'
#' This function provides a complete workflow for processing spatial transcriptomics data.
#' It reads count matrices and spatial image information, constructs a Seurat object,
#' performs normalization using SCTransform, reduces dimensionality, identifies clusters,
#' and generates UMAP embeddings for visualization. The output Seurat object integrates both
#' gene expression and spatial image information for downstream analyses.
#'
#' @param counts A count matrix with peaks/bins as rows and spots as columns.
#' @param counts A data frame containing spot-level metadata, including coordinates.
#' @param min.cells Minimum number of cells required to keep a feature (default: 3).
#' @param min.features Minimum number of features required to keep a cell (default: 0).
#' @param dims Dimensions to use for PCA/Neighbors/UMAP (default: 1:50).
#' @param resolution Resolution parameter for clustering (default: 0.8).
#' @param verbose Logical; whether to print progress messages (default: FALSE).
#'
#' @return A Seurat object with spatial image and analysis results.
#'
#' @export
runSpatial <- function(counts,
                           position,
                           min.cells = 3,
                           min.features = 0,
                           dims_Neighbors = 1:50,
                           dims_UMAP = 1:50,
                           resolution = 0.8,
                           verbose = FALSE) {


  # ---- Create Seurat object with spatial assay ----
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features, meta.data = position,
    assay = "Spatial"
  )

  # ---- Preprocessing: normalization, dimensionality reduction, Graph construction, clustering, and UMAP visualization ----
  seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = verbose)
  seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = verbose)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = dims_Neighbors, verbose = verbose)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = verbose)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = dims_UMAP, verbose = verbose)

  return(seurat_obj)
}
