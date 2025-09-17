# STADS: A Personalized Computational Method to Use Spatial Transcriptomics to Aid Drug-reposition Recommendation

Drug repurposing is a cost-effective strategy for accelerating therapeutic discovery, yet existing single-cell RNA-seq (scRNA-seq)-based methods often overlook the spatial context critical for capturing tissue-specific drug responses. We introduce STADS (Spatial Transcriptomics to Aid Drug-reposition Strategy), a personalized computational framework that leverages spatial transcriptomics data to improve drug repurposing.

![pipeline](figure/STADS%20pipeline%20cartoon_V8.jpeg)

**Illustration of STADS architecture.** STADS utilizes paired diseased and normal tissues from the same patients as input for its spatial domain identification module. This module first performs batch correction and sample alignment before applying a GCN combined with the coherent point drift (CPD) algorithm to identify corresponding spatial domains across conditions. These paired spatial domains then serve as inputs for the drug repurposing module, which identifies potentially reversible genes by comparing differentially expressed genes (DEGs) between spatial domains and integrating drug perturbation data from the L1000 dataset. To prioritize key reversible genes, STADS leverages weights extracted from an XGBoost model trained on potential drug information retrieved from GPT-4o. Additionally, STADS accounts for spatial domain interactions in its drug score calculation. The final drug score is computed by integrating spatial domain proportions and interactions, the significance and weighted influence of reversible genes from XGBoost, as well as drug side effect profiles and sensitivity data. The potential drugs are then validated using empirical evidence from literature and clinical trials, LLM-validated potential drug information, in-silico validation using EHR, and in-vitro validation using cell line experiments.

# Installation

## Docker

TBD

## Manual install

STADS consists of a spatial domain matching module and a drug score calculation module, written in Python and R correspondingly. For a manual installation and custom package usage, please follow the directions here.

The package is tested under Python `3.11.7` and R `4.3.1`.

Install Python requirements using `pip`:

```bash
# To be added
pip install -r requirements.txt
```

Install R environments using `renv`:

```r
# To be added
install.packages("renv")
renv::init()
renv::activate()
renv::install("./")
```

# Tutorial

## Data preparation

Download input data required for STDrug input from Dropbox (TBA). After downloading, the folder should have the following files:

- List of input files

Download sample data analyzed in the manuscript: HCC liver cancer (TBA) and PCa prostate cancer (TBA).

## Spatial domain identification

The first step of STDrug is to identify spatial domains that match patient tumor tissue and adjacent normal tissue. Run Python script using

```bash
python3 -u domain_matching/main.py \
  --data-name '<name of dataset>' \
  --nclust '<number of spatial domains>' \
  --output '<output directory>'
```

The parameters are defined as:

| Option    | Description                                            | Example      |
|-----------|--------------------------------------------------------|--------------|
| data-name | Name of sample dataset (Custom dataset TBA)            | hcc/prostate |
| nclust    | Number of spatial domains                              | 5            |
| output    | Directory to save spatial domains and checkpoint files | ./output     |

This module should produce output files in the following structure:

```
./output
|-- checkpoint
|  |-- stads_cluster.h5ad // AnnData of integrated spatial data with spatial clustering
|-- partition.csv // Spatial domain annotation and meta data
```

## Drug repurposing

Following the spatial domain identification module, STDrug uses a comprehensive drug ranking algorithm to repurpose drugs personalized for each patient. In this step, use R script to run the module

```bash
# TODO merge input path under one folder
Rscript drug_score/Main.R \
  --data-name '<name of dataset>' \
  --tissue '<tissue kind>' \
  --output '<output directory>' \
  --lincs-phase1-path '<>' \
  --lincs-phase2-path '<>' \
  --tahoe-path '<>' \
  --cluster-output-path '<output directory of spatial domains>' \
  --checkpoint-path '<checkpoint save directory>' \
  --drug-annotation-path '<>' \
  --gdsc-path '<>' \
  --sider-path '<>'
```

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