# STDrug: A Computational Method to Use Spatial Transcriptomics to Aid Personalized Drug-reposition Recommendation

Drug repurposing is a cost-effective strategy for accelerating therapeutic discovery, yet existing single-cell RNA-seq (scRNA-seq)-based methods often overlook the spatial context critical for capturing tissue-specific drug responses. We introduce STADS (Spatial Transcriptomics to Aid Drug-reposition Strategy), a personalized computational framework that leverages spatial transcriptomics data to improve drug repurposing.

![pipeline](figure/STADS%20pipeline%20cartoon_V8.jpeg)

**Illustration of STADS architecture.** STADS utilizes paired diseased and normal tissues from the same patients as input for its spatial domain identification module. This module first performs batch correction and sample alignment before applying a GCN combined with the coherent point drift (CPD) algorithm to identify corresponding spatial domains across conditions. These paired spatial domains then serve as inputs for the drug repurposing module, which identifies potentially reversible genes by comparing differentially expressed genes (DEGs) between spatial domains and integrating drug perturbation data from the L1000 dataset. To prioritize key reversible genes, STADS leverages weights extracted from an XGBoost model trained on potential drug information retrieved from GPT-4o. Additionally, STADS accounts for spatial domain interactions in its drug score calculation. The final drug score is computed by integrating spatial domain proportions and interactions, the significance and weighted influence of reversible genes from XGBoost, as well as drug side effect profiles and sensitivity data. The potential drugs are then validated using empirical evidence from literature and clinical trials, LLM-validated potential drug information, in-silico validation using EHR, and in-vitro validation using cell line experiments.

# Installation

## Docker

TBD

## Manual install

STDrug consists of a spatial domain matching module and a drug score calculation module, written in Python and R correspondingly. For a manual installation and custom package usage, please follow the directions here.

The package is tested under Python `3.12.10` and R `4.5.2`.

Install Python requirements using `pip`:

```bash
pip install 'git+https://github.com/remisiki/STdrug.git'
```

Install R environments using `devtools` and `BiocManager`:
```r
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("cmapR", "limma"))
devtools::install_github(c("immunogenomics/presto", "jinworks/CellChat"))
devtools::install_github("remisiki/STdrug")
```

Alternatively, use `renv`:
```r
install.packages("renv")
renv::init(bioconductor=T)
# Restart R session
renv::install("remisiki/STdrug")
```

# Tutorial

## Data preparation

Download drug reference data required for STDrug input from [Dropbox](https://www.dropbox.com/scl/fo/sc7tyjuw9k5ci5v0svyw7/AIya9k7PQH786X8bleDW7KY?rlkey=j1fzh131dyl0xffy5pa6j26dl&st=5ti8ct69&dl=0). It is recommended to create a folder named `data` and extract the reference files under `data/reference`. After downloading, the folder should have the following structure:

```
data
└── reference
    ├── drug_validation
    │   ├── liver.csv
    │   └── prostate.csv
    ├── l1000
    │   ├── GSE70138.tar.gz
    │   └── GSE92742.tar.gz
    ├── tahoe
    │   └── drug_ref.rds
    ├── fda.txt
    ├── gdsc.csv
    └── sider.csv
```

Extract tarball files using `tar`:
```bash
tar -xzvf data/reference/l1000/GSE70138.tar.gz -C data/reference/l1000
tar -xzvf data/reference/l1000/GSE92742.tar.gz -C data/reference/l1000
```

After extraction, the folder structure should have the following structure:
```
data
└── reference
    ├── drug_validation
    │   ├── liver.csv
    │   └── prostate.csv
    ├── l1000
    │   ├── GSE70138
    │   │   ├── GSE70138_Broad_LINCS_cell_info_2017-04-28.txt
    │   │   ├── GSE70138_Broad_LINCS_gene_info_2017-03-06.txt
    │   │   ├── GSE70138_Broad_LINCS_inst_info_2017-03-06.txt
    │   │   ├── GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx
    │   │   └── GSE70138_Broad_LINCS_sig_info_2017-03-06.txt
    │   └── GSE92742
    │       ├── GSE92742_Broad_LINCS_cell_info.txt
    │       ├── GSE92742_Broad_LINCS_gene_info.txt
    │       ├── GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
    │       └── GSE92742_Broad_LINCS_sig_info.txt
    ├── tahoe
    │   └── drug_ref.rds
    ├── fda.txt
    ├── gdsc.csv
    └── sider.csv
```

(Optional) Download the sample data analyzed in the manuscript from [Dropbox](https://www.dropbox.com/scl/fo/momgw9d38zx60d8t5sta8/ALdWi7fPJH1MSuVWYPvqFwc?rlkey=5w1fx7ft7fwlg03t5yk013d5w&st=ai5chwd7&dl=0). You can also put them under `data`.
```
data
├── HCC01N.h5ad
├── HCC01N.rds
├── HCC01T.h5ad
├── HCC01T.rds
├── HCC02N.h5ad
├── HCC02N.rds
├── HCC02T.h5ad
├── HCC02T.rds
├── HCC03N.h5ad
├── HCC03N.rds
├── HCC03T.h5ad
├── HCC03T.rds
├── HCC04N.h5ad
├── HCC04N.rds
├── HCC04T.h5ad
└── HCC04T.rds
```

# Quick Start

## Spatial domain identification

The first step of STDrug is to identify spatial domains that match patient tumor tissue and adjacent normal tissue. Run Python script following the tutorial [Identify spatial domains using STDrug for multiple samples](https://remisiki.github.io/STdrug/spatial-domain-identification-example).

This module should produce output files in the following structure:

```
./output
|-- checkpoint
|  |-- stads_cluster.h5ad // AnnData of integrated spatial data with spatial clustering
|-- partition.csv // Spatial domain annotation and meta data
```

## Drug repurposing

Following the spatial domain identification module, STDrug uses a comprehensive drug ranking algorithm to repurpose drugs personalized for each patient. In this step, use R script to run the module following the tutorial [Use STDrug to calculate drug score for multiple samples](https://remisiki.github.io/STdrug/drug-score-calculation-example).

STDrug generates drug outputs structured as follows. The repurposed top drugs can be inspected from `./output/drugs_<patient>.csv`.

```
./output
|-- checkpoint
|  |-- cci_ratio_<patient>.csv // Cell-cell interation results for patient
|  |-- drug_scores_<patient>.csv // Spatial domain specific drug score results for patient
|  |-- stads_cluster.h5ad // AnnData of integrated spatial data with spatial clustering
|  |-- stads_cluster.rds // Seurat object of integrated spatial data with spatial clustering
|-- drugs_<patient>.csv // Drug score and ranking for patient, higher drug score means better treatment potential
|-- partition.csv // Spatial domain annotation and meta data

```

# Contributor

# Citation