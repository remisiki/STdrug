source("/home/wanxingz/project/STDrug/STDrug/drug_score_intg/LoadLibHdf5.R")

mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  } else {
    message(paste("Path", path, "already exists"))
  }
}

subsetSeurat <- function(object, key = NULL, value = NULL, operator = "eq") {
  if (!is.null(key)) {
    expr <- Seurat::FetchData(object, key, clean = "none")[, 1]
    object <- if (operator == "eq") {
      object[, which(expr == value)]
    } else if (operator == "in") {
      object[, which(expr %in% value)]
    } else if (operator == "notin") {
      object[, which(!expr %in% value)]
    } else if (operator == "notna") {
      object[, which(!is.na(expr))]
    } else if (operator == "neq") {
      object[, which(expr != value)]
    } else if (operator == "lt") {
      object[, which(expr < value)]
    } else if (operator == "gt") {
      object[, which(expr > value)]
    } else if (operator == "leq") {
      object[, which(expr <= value)]
    } else if (operator == "geq") {
      object[, which(expr >= value)]
    } else {
      object
    }
  }
  for (image in names(object@images)) {
    if (
      ("coordinates" %in% slotNames(object@images[[image]])) &&
      (nrow(object@images[[image]]@coordinates) == 0)
    ) {
      object@images[[image]] <- NULL
    }
  }
  object
}

mergeSeuratList <- function(objects, merge_reduction = T) {
  if (length(objects) > 1) {
    object <- merge(objects[[1]], objects[-1])
    if (merge_reduction) {
      for (reduction in names(objects[[1]]@reductions)) {
        object@reductions[[reduction]] <- merge(objects[[1]]@reductions[[reduction]], lapply(objects[-1], function(x) x@reductions[[reduction]]))
      }
    }
  } else {
    object <- objects[[1]]
  }
  object
}

progressbar <- function() {
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  for (i in 1:100) {
    Sys.sleep(0.1)  # Simulate work
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

loadVisium <- function(
  visium_path,
  batch = "visium",
  patient = "patient",
  condition = "N",
  ...
) {
  object <- Seurat::Load10X_Spatial(visium_path, slice = batch, ...)
  object$batch <- batch
  object$patient <- patient
  object$condition <- condition
  # Read the correct spot diameter for fullres
  # https://github.com/satijalab/seurat/issues/7706
  # scalefactors <- jsonlite::read_json(file.path(visium_path, "spatial", "scalefactors_json.json"))
  # object@images[[1]]@scale.factors$spot <- scalefactors$spot_diameter_fullres
  return(object)
}

preprocessQualityControl <- function(
  object,
  min_counts = NULL,
  min_genes = NULL,
  percent_mt = NULL
) {
  # Seurat qc for Visium data
  n_cell <- ncol(object)
  # Filter cells by min_genes
  if (!is.null(min_genes)) {
    object <- object[, object$nFeature_Spatial >= min_genes]
  }
  # Filter cells by min_counts
  if (!is.null(min_counts)) {
    object <- object[, object$nCount_Spatial >= min_counts]
  }
  # Filter cells by percent_mt
  if (!is.null(percent_mt)) {
    # Calculate MT gene percentage
    object$percent_mt <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
    object <- object[, object$percent_mt <= percent_mt]
  }
  n_cell_qc <- ncol(object)
  message(paste0(n_cell_qc, "/", n_cell, " cells remained after quality control"))
  return(object)
}

preprocessNormalize <- function(object) {
  # object <- Seurat::SCTransform(object, assay = "Spatial")
  object <- Seurat::NormalizeData(object, assay = "Spatial")
  Seurat::DefaultAssay(object) <- "Spatial"
  return(object)
}

lsplit <- function(...) {
  return(as.list(unlist(strsplit(...))))
}

dropUnusedLevels <- function(obj, col) {
  if (is.factor(obj@meta.data[[col]])) {
    obj@meta.data[[col]] <- droplevels(obj@meta.data[[col]])
  }
  return(obj)
}

minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

renameSeuratAssay <- function(obj, from_assay = NULL, to_assay = NULL) {
  if (is.null(from_assay)) {
    from_assay <- obj@active.assay
  }
  if (!is.null(to_assay)) {
    obj <- SeuratObject::RenameAssays(obj, setNames(to_assay, from_assay))
  }
  obj
}

toTitleCase <- function(s) {
  sapply(s, function(x) {
    l <- nchar(x)
    if (l >= 2) {
      paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, l)))
    } else if (l >= 1) {
      toupper(x)
    } else {
      x
    }
  }, USE.NAMES = F)
}

grepSeuratGenes <- function(object, pattern, ignore_case = T) {
  rownames(object)[grep(pattern, rownames(object), ignore.case = ignore_case)]
}

revFactorOrder <- function(x) {
  if (!is.factor(x)) {
    stop("Input vector is not a factor")
  }
  factor(x, levels = rev(levels(x)))
}

minmax <- function(x, a = 0, b = 1) {
  a + (x - min(x)) * (b - a) / (max(x) - min(x))
}