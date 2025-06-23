source("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/stable/drug_score/DrugScore.R")
source("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/stable/drug_score/DataLoader.R")
source("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/stable/drug_score/L1000Loader.R")

spec <- matrix(
  c(
    "data-name", "n", 1, "character",
    "patients", "p", 1, "character",
    "tissue", "t", 1, "character",
    "output", "o", 1, "character",
    "lincs-phase1-path", "l1", 1, "character",
    "lincs-phase2-path", "l2", 1, "character",
    "cluster-output-path", "co", 1, "character",
    "checkpoint-path", "cp", 1, "character",
    "drug-annotation-path", "dp", 1, "character",
    "gdsc-path", "gp", 1, "character",
    "sider-path", "sp", 1, "character",
    "drop-celltype", "d", 1, "character",
    "deconv-result-path", "dvp", 1, "character",
    "sc-ref-path", "scp", 1, "character"
  ),
  byrow = T,
  ncol = 4
)

opt <- getopt::getopt(spec)

data_name <- opt[["data-name"]]
patients <- opt[["patients"]]
if (!is.null(patients)) {
  patients <- unlist(lsplit(patients, "\\|"))
}
tissue <- opt[["tissue"]]
output_dir <- opt[["output"]]
lincs_drug_response_phase1_path <- opt[["lincs-phase1-path"]]
lincs_drug_response_phase2_path <- opt[["lincs-phase2-path"]]
cluster_output_path <- opt[["cluster-output-path"]]
checkpoint_path <- opt[["checkpoint-path"]]
drug_annotation_path <- opt[["drug-annotation-path"]]
gdsc_path <- opt[["gdsc-path"]]
sider_path <- opt[["sider-path"]]
drop_celltype <- opt[["drop-celltype"]]
if (!is.null(drop_celltype)) {
  drop_celltype <- unlist(lsplit(drop_celltype, "\\|"))
}
deconv_result_path <- opt[["deconv-result-path"]]
sc_ref_path <- opt[["sc-ref-path"]]

library(Seurat)
set.seed(0)

mkdir(output_dir)
if (!is.null(checkpoint_path)) {
  mkdir(checkpoint_path)
}

# Load data
data_factory <- loadSampleData(
  cluster_output_path = cluster_output_path,
  data_name = data_name,
  drop_celltype = drop_celltype,
  deconv_result_path = deconv_result_path,
  sc_ref_path = sc_ref_path,
  patients = patients
)
patients <- data_factory$patients
data_list <- data_factory$data
domains <- data_factory$domains

if (!is.null(checkpoint_path)) {
  # Save stads clustered data into integrated seurat object
  cluster_rds_file <- file.path(checkpoint_path, "stads_cluster.rds")
  message(paste0("Save cluster result into seurat object ", cluster_rds_file))
  data_list_copy <- data_list
  for (patient in patients) {
    data_list_copy[[patient]] <- mergeSeuratList(data_list_copy[[patient]])
  }
  data_list_copy <- mergeSeuratList(data_list_copy)
  saveRDS(data_list_copy, cluster_rds_file)
}

# Load L1000, gpt-4o drug annotations, GDSC, Sider
l1000 <- loadLincsTwoPhasesDrugResponse(lincs_drug_response_phase1_path, lincs_drug_response_phase2_path, tissue)
drug_annotation <- read.csv(drug_annotation_path, row.names = 1)
gdsc <- read.csv(gdsc_path, row.names = 1)
sider <- read.csv(sider_path, row.names = 1)

# Calculate drug score
for (patient in patients) {
  message(paste0("Calculate drug score for ", patient))
  # Create normal and tumor samples
  tumor_samples <- data_list[[patient]]$tumor
  normal_samples <- data_list[[patient]]$normal
  drug_scores <-
    calcDrugScore(
      tissue,
      tumor_samples,
      normal_samples,
      patient,
      domains,
      l1000,
      drug_annotation,
      gdsc,
      sider,
      checkpoint_path = checkpoint_path
    )
}
