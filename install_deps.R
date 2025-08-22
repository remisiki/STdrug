# run in Rscript or R
source("renv/activate.R")

# 1) Pin Bioconductor for R 4.3.x (don’t set repositories; that setter doesn’t exist in renv 1.0.7)
renv::settings$bioconductor.version("3.18")

# 2) Use a reliable CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 3) Fix the broken cache symlinks
renv::repair(prompt = FALSE)

# 4) Restore ONLY what your pipeline needs (avoid the huge full restore)
needed <- c(
  "getopt",
  "Seurat","SeuratObject","Matrix","limma","xgboost",
  "SingleCellExperiment","SummarizedExperiment"
  # add "CellChat" ONLY if you truly need calcCciRatio(); otherwise skip it
  # "CellChat"
)
renv::restore(packages = needed, prompt = FALSE)

library(getopt); library(Seurat); library(SingleCellExperiment)
cat("renv minimal restore OK\n")