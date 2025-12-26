import os
from pathlib import Path
from typing import *

import SpaGCN as spg
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import torch
from pycpd import DeformableRegistration

from . import Cluster, SpatialGraph


class CoherentPointDrift(DeformableRegistration):
    def __init__(self, verbose: bool = True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.verbose = verbose

    def iterate(self):
        self.expectation()
        self.maximization()
        if self.verbose:
            print(f"Iteration {self.iteration}/{self.max_iterations}, diff={self.diff}")
        self.iteration += 1


def buildSpatialGraph(adata: ad.AnnData):
    # Embedding
    spatial_graphs = {}
    kernel_widths = {}
    batches = adata.obs["batch"].unique()
    for batch in batches:
        print(f"Build graph on batch {batch}")
        adata_batch = adata[adata.obs["batch"] == batch].copy()

        # Construct the spatial graph
        x_pixel = adata_batch.obsm["spatial"][:, 0]
        y_pixel = adata_batch.obsm["spatial"][:, 1]

        # Calculate spatial distances between samples
        spatial_graphs[batch] = spg.calculate_adj_matrix(
            x_pixel, y_pixel, histology=False
        )

        # Choose RBF kernel width
        end = 1000
        while True:
            l = spg.search_l(0.5, spatial_graphs[batch], start=1e-4, end=end)
            if l is not None:
                break
            else:
                end *= 2
        kernel_widths[batch] = l

    print("Spatial graph construction is done")
    return spatial_graphs, kernel_widths


def trainGcn(
    adata: ad.AnnData,
    spatial_graphs,
    kernel_widths,
    nclust: int,
    seed: float = 100,
    use_rep: str = "X_pca",
    key_added: str = "X_spagcn",
):
    # Run SpaGCN for each batch
    adata.obsm[key_added] = np.zeros(adata.obsm[use_rep].shape)
    batches = adata.obs["batch"].unique()
    for batch in batches:
        print(f"Train GCN on batch {batch}")
        adata_batch = adata[adata.obs["batch"] == batch].copy()

        # Search for best resolution
        res = SpatialGraph.searchResolution(
            adata_batch,
            spatial_graphs[batch],
            kernel_widths[batch],
            nclust,
            r_seed=seed,
            t_seed=seed,
            n_seed=seed,
        )

        # Init SpaGCN
        clf = SpatialGraph.SpaGCN()
        clf.set_l(kernel_widths[batch])
        # ----------Start training----------
        adj = spatial_graphs[batch]
        # Use custom node embeddings
        clf.embed = adata_batch.obsm[use_rep]
        clf.adj_exp = np.exp(-1 * (adj**2) / (2 * (clf.l**2)))
        clf.model = spg.models.simple_GC_DEC(clf.embed.shape[1], clf.embed.shape[1])
        clf.model.fit(
            clf.embed,
            clf.adj_exp,
            weight_decay=0,
            opt="admin",
            res=res,
        )
        # ----------End training----------
        with torch.no_grad():
            spagcn_embed = clf.model.gc(
                torch.FloatTensor(clf.embed), torch.FloatTensor(clf.adj_exp)
            ).numpy()
        adata.obsm[key_added][adata.obs["batch"] == batch,] = spagcn_embed.copy()
    return adata


def coherentPointDrift(
    adata: ad.AnnData,
    use_rep: str,
    key_added: Optional[str] = None,
    verbose: bool = True,
    **kwargs,
):
    if key_added is None:
        key_added = use_rep + "_cpd"
    adata.obsm[key_added] = adata.obsm[use_rep].copy()
    patients = adata.obs["patient"].cat.categories
    for patient in patients:
        if verbose:
            print(f"Coherent point drift on patient {patient}")
        # Subset data by patient
        patient_mask = adata.obs["patient"] == patient
        # Subset data to case and control
        case_mask = patient_mask & (adata.obs["condition"] == "T")
        control_mask = patient_mask & (adata.obs["condition"] == "N")
        # Extract original embeddings
        case_embed = adata.obsm[key_added][case_mask].copy()
        control_embed = adata.obsm[key_added][control_mask].copy()
        # Apply deformable transformation on case embeddings
        reg = CoherentPointDrift(
            verbose=verbose, X=control_embed, Y=case_embed, **kwargs
        )
        transformed_embed, _ = reg.register()
        # Assign new embeddings
        adata.obsm[key_added][case_mask] = transformed_embed.copy()
    return adata


def findSpatialDomains(
    tumor_samples: ad.AnnData | List[ad.AnnData],
    normal_samples: ad.AnnData | List[ad.AnnData],
    nclust: int,
    output_dir: Optional[Path] = None,
    tolerance: float = 0.001,
    seed: int = 0,
) -> pd.DataFrame:
    """
    Identifies spatial domains within tumor and normal samples using STDrug clustering.

    Args:
        tumor_samples (ad.AnnData | List[ad.AnnData]): Spatal data for tumor samples in AnnData format.
        normal_samples (ad.AnnData | List[ad.AnnData]): Spatial data for normal samples in AnnData format.
        nclust (int): The number of clusters to form during the domain identification process.
        output_dir (Optional[Path], optional): Directory to save output files. Defaults to None, which means no output will be saved.
        tolerance (float, optional): The tolerance level for coherent point drift distance. Defaults to 0.001.
        seed (int, optional): Random seed for reproducibility. Defaults to 0.
    """

    checkpoint_dir = None
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        checkpoint_dir = Path(output_dir, "checkpoint")
        os.makedirs(checkpoint_dir, exist_ok=True)

    # Check data format
    print("Check data format")
    if len(tumor_samples) != len(normal_samples):
        raise Exception(
            "The number of tumor samples must match the number of normal samples"
        )
    for i in range(len(tumor_samples)):
        if "batch" not in tumor_samples[i].obs.columns:
            tumor_samples[i].obs["batch"] = f"sample{i}t"
        if "patient" not in tumor_samples[i].obs.columns:
            tumor_samples[i].obs["patient"] = "patient"
        if "condition" not in tumor_samples[i].obs.columns:
            tumor_samples[i].obs["condition"] = "tumor"
        if "batch" not in normal_samples[i].obs.columns:
            normal_samples[i].obs["batch"] = f"sample{i}n"
        if "patient" not in normal_samples[i].obs.columns:
            normal_samples[i].obs["patient"] = "patient"
        if "condition" not in normal_samples[i].obs.columns:
            normal_samples[i].obs["condition"] = "normal"

    # Preprocess
    # Concat data
    adata = [*tumor_samples, *normal_samples]
    adata = ad.concat(adata, uns_merge="unique")
    adata.obs_names_make_unique()
    # Normalize
    print("Normalize data")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # PCA
    print("Run PCA")
    sc.pp.pca(adata, random_state=seed)

    # Run harmony
    print("Run harmony")
    sce.pp.harmony_integrate(
        adata,
        key=["batch"],
        basis="X_pca",
        adjusted_basis="X_harmony",
        max_iter_harmony=50,
        random_state=seed,
    )

    # Train GCN
    print("Train GCN")
    spatial_graphs, kernel_widths = buildSpatialGraph(adata)
    adata = trainGcn(
        adata,
        spatial_graphs,
        kernel_widths,
        nclust=nclust,
        seed=(seed + 100),
        use_rep="X_harmony",
    )

    # Clustering
    print("Find clusters")
    # Run umap on gcn embeddings
    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        use_rep="X_spagcn",
        key_added="stads",
        random_state=seed,
    )
    sc.tl.umap(adata, neighbors_key="stads", min_dist=0.1, spread=2, random_state=seed)
    adata.obsm["X_umap_stads"] = adata.obsm["X_umap"].copy()
    # CPD
    adata = coherentPointDrift(
        adata, use_rep="X_umap_stads", alpha=2, beta=5, tolerance=tolerance
    )
    # Leiden
    adata = Cluster.findCluster(
        adata,
        nclust=nclust,
        use_rep="X_umap_stads_cpd",
        flavor="igraph",
        res_max=1,
        seed=seed,
    )

    # Save data
    print(f"Save results to {output_dir}")
    if checkpoint_dir is not None:
        adata.write_h5ad(Path(checkpoint_dir, "stads_cluster.h5ad"))
    result_df = adata.obs[["batch", "cluster_stads"]]
    result_df.columns = ["batch", "cluster"]
    if output_dir is not None:
        result_df.to_csv(Path(output_dir, "partition.csv"))

    return result_df
