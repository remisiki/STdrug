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

# Contributor

# Citation