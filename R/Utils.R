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

lsplit <- function(...) {
  return(as.list(unlist(strsplit(...))))
}

dropUnusedLevels <- function(obj, col) {
  if (is.factor(obj@meta.data[[col]])) {
    obj@meta.data[[col]] <- droplevels(obj@meta.data[[col]])
  }
  return(obj)
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