source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

saveHccData <- function() {
  data_dir <- "/nfs/dcmb-lgarmire/shared/public/PMC8683021"
  output_dir <- data_dir
  patients <- c("HCC01", "HCC02", "HCC03", "HCC04")
  conditions <- c("N", "T")
  for (patient in patients) {
    for (condition in conditions) {
      batch <- paste0(patient, condition)
      message(paste0("Reading raw data for ", batch))
      visium_path <- file.path(data_dir, patient, condition)
      data <- loadVisium(visium_path, batch, patient, condition)
      data <- preprocessQualityControl(data, min_genes = 3)
      data <- preprocessNormalize(data)
      # Save data
      output_file <- file.path(output_dir, paste0(batch, ".rds"))
      message(paste0("Save to ", output_file))
      saveRDS(data, output_file)
    }
  }
}

saveProstateData <- function() {
  data_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/prostate/erickson")
  output_dir <- data_dir
  
  batch_map <- list(
    patient1 = list("H1_1", "H1_2", "H1_4", "H1_5", "H2_1", "H2_2", "H2_5", "V1_1", "V1_2"),
    patient2 = list("H2_1", "H2_2", "H3_1", "H3_2", "H3_4", "H3_5", "H3_6", "V1_1", "V1_2", "V1_3", "V1_4", "V1_5", "V1_6", "V2_1", "V2_2")
  )
  tumor_cases <- list("patient1_V2_1", "patient1_H2_1", "patient1_H2_2", "patient1_H2_5", "patient1_H1_1", "patient1_H1_2", "patient1_H1_4", "patient1_H1_5", "patient2_H2_1", "patient2_H2_2", "patient2_H3_1", "patient2_H3_6")
  for (patient in names(batch_map)) {
    integrated_data <- list()
    for (batch in batch_map[[patient]]) {
      batch <- paste0(patient, "_", batch)
      message(paste0("Reading raw data for ", batch))
      visium_path <- file.path(data_dir, batch)
      condition <- ifelse(batch %in% tumor_cases, "T", "N")
      data <- loadVisium(visium_path, batch, patient, condition)
      data <- preprocessQualityControl(data, min_genes = 3)
      integrated_data[[batch]] <- data
    }
    # Merge slices from one patient
    integrated_data <- mergeSeuratList(integrated_data)
    integrated_data <- preprocessNormalize(integrated_data)
    # Save data
    output_file <- file.path(output_dir, paste0(patient, ".rds"))
    message(paste0("Save to ", output_file))
    saveRDS(integrated_data, output_file)
  }
}

saveNsclcData <- function() {
  data_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/nsclc")
  output_dir <- data_dir
  
  batch_map <- list(
    P10 = list("B1", "B2", "T1", "T2", "T3", "T4"),
    P11 = list("B1", "B2", "T1", "T2", "T3", "T4"),
    P15 = list("B1", "B2", "T1", "T2"),
    P16 = list("B1", "B2", "T1", "T2"),
    P17 = list("B1", "B2", "T1", "T2"),
    P19 = list("B1", "B2", "T1", "T2"),
    P24 = list("B1", "B2", "T1", "T2"),
    P25 = list("B1", "B2", "T1", "T2")
  )
  tumor_cases <- list("P10_T1", "P10_T2", "P10_T3", "P10_T4", "P11_T1", "P11_T2", "P11_T3", "P11_T4", "P15_T1", "P15_T2", "P16_T1", "P16_T2", "P17_T1", "P17_T2", "P19_T1", "P19_T2", "P24_T1", "P24_T2", "P25_T1", "P25_T2")
  for (patient in names(batch_map)) {
    integrated_data <- list()
    for (batch in batch_map[[patient]]) {
      batch <- paste0(patient, "_", batch)
      message(paste0("Reading raw data for ", batch))
      visium_path <- file.path(data_dir, batch)
      condition <- ifelse(batch %in% tumor_cases, "T", "N")
      data <- loadVisium(visium_path, batch, patient, condition)
      data <- preprocessQualityControl(data, min_genes = 3)
      integrated_data[[batch]] <- data
    }
    # Merge slices from one patient
    integrated_data <- mergeSeuratList(integrated_data)
    integrated_data <- preprocessNormalize(integrated_data)
    # Save data
    output_file <- file.path(output_dir, paste0(patient, ".rds"))
    message(paste0("Save to ", output_file))
    saveRDS(integrated_data, output_file)
  }
}

removeCelltype <- function(
  data,
  ref,
  cell_abundance,
  drop_celltype,
  use_celltype = "celltype"
) {
  # Make a copy of original data to prevent direct writing
  data_cp <- data
  ref_cp <- ref
  # Subset data and ref by common genes
  common_genes <- intersect(rownames(ref_cp), rownames(data_cp))
  data_cp <- data_cp[common_genes, ]
  ref_cp <- ref_cp[common_genes, ]
  # Normalize target data and ref
  counts <- cbind(Seurat::GetAssayData(data_cp, slot = "counts"), Seurat::GetAssayData(ref_cp, slot = "counts"))
  normalized_data <- Seurat::NormalizeData(counts)
  data_cp <- Seurat::SetAssayData(data_cp, "data", normalized_data[, 1:ncol(data_cp)])
  ref_cp <- Seurat::SetAssayData(ref_cp, "data", normalized_data[, (ncol(data_cp) + 1):ncol(normalized_data)])
  # Calculate mean expression of cell types to drop in ref
  ref_celltype_mean <- rowMeans(Seurat::GetAssayData(
    subsetSeurat(ref_cp, use_celltype, drop_celltype, operator = "in")
  ))
  # Subset cell abundance by all cell types in ref
  celltypes <- unique(ref_cp@meta.data[[use_celltype]])
  cell_abundance <- cell_abundance[, celltypes]
  # Convert cell abundance to proportion
  cell_abundance <- t(apply(cell_abundance, 1, function(x)
    x / sum(x)))
  # Sum cell type proportions to drop
  cell_abundance <- cell_abundance[, drop_celltype, drop = F]
  cell_abundance <- rowSums(cell_abundance)
  # Calculate the amount of expression to subtract from target data using mean
  # expression in ref multiplied by cell type proportion
  pseudo_bulk <- outer(ref_celltype_mean, cell_abundance)
  # Correct the assay, only modify non zero entries
  assay <- Seurat::GetAssayData(data_cp)
  non_zero_indices <- Matrix::which(assay != 0, arr.ind = T)
  pseudo_bulk <- Matrix::sparseMatrix(
    i = non_zero_indices[, 1],
    j = non_zero_indices[, 2],
    x = pseudo_bulk[non_zero_indices],
    dims = c(nrow(pseudo_bulk), ncol(pseudo_bulk))
  )
  corrected_assay <- assay - pseudo_bulk
  # If corrected expression is negative, adjust to zero
  corrected_assay[corrected_assay < 0] <- 0
  # Assign the corrected assay to data
  data_cp <- Seurat::SetAssayData(object = data_cp, new.data = corrected_assay)
  data_cp
}

loadSampleData <- function(
  data_modality = "spatial",
  patients = NULL,
  cluster_output_path = NULL,
  domain_key = "stads_domain",
  data_name = "hcc",
  drop_celltype = NULL,
  deconv_result_path = NULL,
  sc_ref_path = NULL
) {
  data_list <- list()
  partition <- NULL
  domains <- NULL
  if (data_modality == "spatial") {
    if (!is.null(cluster_output_path)) {
      partition_file <- file.path(cluster_output_path, "partition.csv")
      if (file.exists(partition_file)) {
        message(paste0("Read STADS clustering result from ", partition_file))
        partition <- read.csv(partition_file, row.names = 1)
        domains <- unique(partition$cluster)
      } else {
        stop(paste0("Invalid STADS clustering result file ", partition_file))
      }
    }
    if (!is.null(sc_ref_path)) {
      sc_ref <- readRDS(sc_ref_path)
    }
    if (data_name == "hcc") {
      message("Load HCC dataset")
      data_dir <- "/nfs/dcmb-lgarmire/shared/public/PMC8683021"
      patients <- c("HCC01", "HCC02", "HCC03", "HCC04")
      for (patient in patients) {
        message(patient)
        normal_samples <- NULL
        tumor_samples <- NULL
        for (condition in c("N", "T")) {
          batch <- paste0(patient, condition)
          data_batch <- readRDS(file.path(data_dir, paste0(patient, condition, ".rds")))
          if (!is.null(drop_celltype)) {
            cell_abundance <- read.csv(
              file.path(deconv_result_path, batch, "cell_abundance.csv"),
              row.names = 1,
              check.names = F
            )
            data_batch <- removeCelltype(data_batch, sc_ref, cell_abundance, drop_celltype)
          }
          if (!is.null(partition)) {
            data_batch@meta.data[[domain_key]] <- partition[partition$batch == batch, "cluster"]
          }
          if (condition == "N") {
            normal_samples <- data_batch
          } else if (condition == "T") {
            tumor_samples <- data_batch
          }
        }
        data_list[[patient]] <- list(normal = normal_samples, tumor = tumor_samples)
      }
    } else if (data_name == "prostate") {
      message("Load Erickson prostate dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/prostate/erickson"
      patients <- c("patient1", "patient2")
      for (patient in patients) {
        message(patient)
        data <- readRDS(file.path(data_dir, paste0(patient, ".rds")))
        normal_batches <- unique(data@meta.data[data$condition == "N", "batch"])
        tumor_batches <- unique(data@meta.data[data$condition == "T", "batch"])
        normal_samples <- list()
        tumor_samples <- list()
        batches <- c(normal_batches, tumor_batches)
        for (batch in batches) {
          data_batch <- subsetSeurat(data, "batch", batch)
          if (!is.null(drop_celltype)) {
            cell_abundance <- read.csv(
              file.path(deconv_result_path, batch, "cell_abundance.csv"),
              row.names = 1,
              check.names = F
            )
            data_batch <- removeCelltype(data_batch, sc_ref, cell_abundance, drop_celltype)
          }
          if (!is.null(partition)) {
            data_batch@meta.data[[domain_key]] <- partition[partition$batch == batch, "cluster"]
          }
          if (batch %in% normal_batches) {
            normal_samples[[batch]] <- data_batch
          } else {
            tumor_samples[[batch]] <- data_batch
          }
        }
        normal_samples <- mergeSeuratList(normal_samples)
        tumor_samples <- mergeSeuratList(tumor_samples)
        data_list[[patient]] <- list(normal = normal_samples, tumor = tumor_samples)
      }
    } else if (data_name == "prostate_hirz") {
      message("Load Hirz prostate dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/prostate/hirz"
      patients <- c("patient1", "patient2")
      for (patient in patients) {
        message(patient)
        data <- readRDS(file.path(data_dir, paste0(patient, ".rds")))
        normal_samples <- subsetSeurat(data, "condition", "N")
        tumor_samples <- subsetSeurat(data, "condition", "T")
        if (!is.null(partition)) {
          batches <- unique(normal_samples$batch)
          normal_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
          batches <- unique(tumor_samples$batch)
          tumor_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
        }
        data_list[[patient]] <- list(normal = normal_samples, tumor = tumor_samples)
      }
    } else if (data_name == "nsclc") {
      message("Load nsclc dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/nsclc"
      patients <- c("P10", "P11", "P15", "P16", "P17", "P19", "P24", "P25")
      for (patient in patients) {
        message(patient)
        data <- readRDS(file.path(data_dir, paste0(patient, ".rds")))
        normal_samples <- subsetSeurat(data, "condition", "N")
        tumor_samples <- subsetSeurat(data, "condition", "T")
        if (!is.null(partition)) {
          batches <- unique(normal_samples$batch)
          normal_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
          batches <- unique(tumor_samples$batch)
          tumor_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
        }
        data_list[[patient]] <- list(normal = normal_samples, tumor = tumor_samples)
      }
    } else if (data_name == "luad") {
      message("Load luad dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/luad"
      patients <- c("P10", "P15", "P16", "P24", "P25")
      for (patient in patients) {
        message(patient)
        data <- readRDS(file.path(data_dir, paste0(patient, ".rds")))
        normal_samples <- subsetSeurat(data, "condition", "N")
        tumor_samples <- subsetSeurat(data, "condition", "T")
        if (!is.null(partition)) {
          batches <- unique(normal_samples$batch)
          normal_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
          batches <- unique(tumor_samples$batch)
          tumor_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
        }
        data_list[[patient]] <- list(normal = normal_samples, tumor = tumor_samples)
      }
    } else if (data_name == "lusc") {
      message("Load lusc dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/lusc"
      patients <- c("P11", "P17", "P19")
      for (patient in patients) {
        message(patient)
        data <- readRDS(file.path(data_dir, paste0(patient, ".rds")))
        normal_samples <- subsetSeurat(data, "condition", "N")
        tumor_samples <- subsetSeurat(data, "condition", "T")
        if (!is.null(partition)) {
          batches <- unique(normal_samples$batch)
          normal_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
          batches <- unique(tumor_samples$batch)
          tumor_samples[[domain_key]] <- partition[partition$batch %in% batches, ]$cluster
        }
        data_list[[patient]] <- list(normal = normal_samples, tumor = tumor_samples)
      }
    }
  } else if (data_modality == "scRNA") {
    if (data_name == "prostate") {
      message("Load prostate cell atlas dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/prostate/prostatecellatlas"
      data <- readRDS(file.path(data_dir, "prostate_ref.rds"))
      data <- renameSeuratAssay(data, to_assay = "RNA")
      patients <- unique(data$patient)
      data$batch <- data$sample
      data$condition <- ifelse(data$group == "normal", "N", "T")
      for (patient in patients) {
        data_patient <- subsetSeurat(data, "patient", patient)
        data_list[[patient]] <- list(
          normal = subsetSeurat(data_patient, "condition", "N"),
          tumor = subsetSeurat(data_patient, "condition", "T")
        )
      }
    } else if (data_name == "pancreas") {
      message("Load pancreas dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/pancreas/GSE212966"
      patients <- c("patient1", "patient2", "patient6")
      for (patient in patients) {
        message(patient)
        normal_samples <- readRDS(file.path(data_dir, paste0(patient, "n_ann.rds")))
        normal_samples$condition <- "N"
        tumor_samples <- readRDS(file.path(data_dir, paste0(patient, "t_ann.rds")))
        tumor_samples$condition <- "T"
        data_list[[patient]] <- list(normal = normal_samples, tumor = tumor_samples)
      }
    } else if (data_name == "hcc") {
      message("Load single cell HCC dataset")
      data_dir <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/hcc/GSE149614"
      data <- readRDS(file.path(data_dir, "data.rds"))
      patients <- unique(data$patient)
      for (patient in patients) {
        data_patient <- subsetSeurat(data, "patient", patient)
        data_list[[patient]] <- list(
          normal = subsetSeurat(data_patient, "condition", "N"),
          tumor = subsetSeurat(data_patient, "condition", "T")
        )
      }
    }
  }
  
  return(list(patients = patients, data = data_list, domains = domains))
}
