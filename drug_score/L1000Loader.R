loadLincsDrugResponse <- function(
  response_data_file,
  sig_info_file,
  cell_lines
) {
  message(paste0("Read signature info from ", sig_info_file))
  sig_info <- read.csv(sig_info_file, sep = "\t")
  # Subset sig info to cell lines of interest
  sig_ids <- sig_info[sig_info$cell_id %in% cell_lines &
                        sig_info$pert_type == "trt_cp" &
                        sig_info$pert_iname %in% Asgard::FDA.drug, "sig_id"]
  sig_info <- sig_info[, c("sig_id", "pert_iname")]
  rownames(sig_info) <- sig_info$sig_id
  sig_info <- sig_info[sig_ids, ]
  
  if (length(sig_ids) == 0) {
    message("No response data found for cell lines")
    return(NULL)
  } else {
    # Load drug response
    expr <- as.data.frame(cmapR::parse_gctx(response_data_file, cid = sig_ids)@mat)
    expr$gene_id <- rownames(expr)
    return(list(
      response = expr,
      sig_info = sig_info
    ))
  }
}

loadLincsTwoPhasesDrugResponse <- function(
  phase1_path,
  phase2_path,
  tissue = "liver"
) {
  # Read common cell line and gene info from phase I (GSE92742) folder
  cell_info_file <- normalizePath(file.path(phase1_path, "GSE92742_Broad_LINCS_cell_info.txt"))
  gene_info_file <- normalizePath(file.path(phase1_path, "GSE92742_Broad_LINCS_gene_info.txt"))
  message(paste0("Read cell info from ", cell_info_file))
  cell_info <- read.csv(cell_info_file, sep = "\t")
  cell_lines <- cell_info[cell_info$primary_site == tissue, ]$cell_id
  message(paste0("Read gene info from ", gene_info_file))
  gene_info <- read.csv(gene_info_file, sep = "\t")
  # Load response for phase I drugs
  message("Load LINCS phase I drug response data (GSE92742)")
  phase1_drugs <- loadLincsDrugResponse(
    response_data_file = file.path(phase1_path, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"),
    sig_info_file = file.path(phase1_path, "GSE92742_Broad_LINCS_sig_info.txt"),
    # gene_info = gene_info,
    cell_lines = cell_lines
  )
  # Load response for phase II drugs
  message("Load LINCS phase II drug response data (GSE70138)")
  phase2_drugs <- loadLincsDrugResponse(
    response_data_file = file.path(phase2_path, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"),
    sig_info_file = file.path(phase2_path, "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt"),
    # gene_info = gene_info,
    cell_lines = cell_lines
  )
  # Merge phase I and II drug responses
  if (is.null(phase1_drugs) && is.null(phase2_drugs)) {
    return(NULL)
  } else {
    if (is.null(phase1_drugs)) {
      drug_response <- phase2_drugs$response
      sig_info <- phase2_drugs$sig_info
    } else if (is.null(phase2_drugs)) {
      drug_response <- phase1_drugs$response
      sig_info <- phase1_drugs$sig_info
    } else {
      drug_response <- merge(phase1_drugs$response, phase2_drugs$response, by = "gene_id")
      sig_info <- rbind(phase1_drugs$sig_info, phase2_drugs$sig_info)
    }
    rownames(drug_response) <- drug_response[, 1]
    drug_response <- drug_response[, -1]
    return(list(
      response = drug_response,
      sig_info = sig_info,
      gene_info = gene_info
    ))
  }
}