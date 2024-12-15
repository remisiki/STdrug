import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from scipy.spatial.distance import pdist, squareform
from typing import *


def countBalancedClusters(adata: ad.AnnData, use_rep: str):
    is_balanced = (
        (
            adata.obs.groupby([use_rep, "patient", "condition"], observed=False).size()
            > 0
        )
        .groupby(level=0, observed=False)
        .all()
    )
    paired_cluster_ids = is_balanced[is_balanced].index.tolist()
    non_paired_cluster_ids = is_balanced[~is_balanced].index.tolist()
    return paired_cluster_ids, non_paired_cluster_ids


def calcCentroid(X: np.ndarray):
    return np.mean(X, axis=0)


def dist(df: pd.DataFrame):
    # Calculate pairwise distances
    distances = pdist(df, metric="euclidean")
    # Convert to square matrix
    return pd.DataFrame(squareform(distances), index=df.index, columns=df.index)


def mergeNonPairedClusters(
    adata: ad.AnnData,
    use_rep: str,
    cluster_key: str,
    from_cluster_ids: List,
    to_cluster_ids: List,
):
    centroids = pd.DataFrame(
        {
            cluster_id: calcCentroid(
                adata.obsm[use_rep][adata.obs[cluster_key] == cluster_id]
            )
            for cluster_id in from_cluster_ids + to_cluster_ids
        }
    ).T
    distances = dist(centroids)
    merge_map = (
        distances.apply(lambda x: x.mask(x.index.isin(from_cluster_ids), np.inf))
        .idxmin()
        .loc[from_cluster_ids]
        .to_dict()
    )
    for from_id, to_id in merge_map.items():
        adata.obs.loc[adata.obs[cluster_key] == from_id, cluster_key] = to_id
    adata.obs[cluster_key] = adata.obs[cluster_key].cat.remove_unused_categories()
    return adata


def balanceClusters(
    data: ad.AnnData, use_rep: str, key_added: str, nclust: int, curr_n: int
):
    is_balanced = False
    paired_cluster_ids, non_paired_cluster_ids = countBalancedClusters(data, key_added)
    n_paired_cluster = len(paired_cluster_ids)
    n_non_paired_cluster = len(non_paired_cluster_ids)
    if n_paired_cluster > nclust:
        cluster_size = data.obs[key_added].value_counts().to_dict()
        paired_cluster_size = {id: cluster_size[id] for id in paired_cluster_ids}
        big_cluster_ids = [
            k
            for k, v in sorted(
                paired_cluster_size.items(), key=lambda item: item[1], reverse=True
            )[:nclust]
        ]
        small_cluster_ids = [k for k in paired_cluster_size if k not in big_cluster_ids]
        paired_cluster_ids = big_cluster_ids
        non_paired_cluster_ids.extend(small_cluster_ids)
        n_paired_cluster = len(paired_cluster_ids)
        n_non_paired_cluster = len(non_paired_cluster_ids)
    if n_paired_cluster == nclust:
        if n_non_paired_cluster == 0:
            print("All clusters are paired")
        else:
            print(
                f"{n_paired_cluster}/{curr_n} clusters are paired, merge non-paired clusters"
            )
            data = mergeNonPairedClusters(
                data,
                use_rep,
                key_added,
                from_cluster_ids=non_paired_cluster_ids,
                to_cluster_ids=paired_cluster_ids,
            )
        is_balanced = True
    else:
        print(f"{n_paired_cluster}/{curr_n} clusters are paired, adjust leiden nclust")
    return is_balanced, data


def findCluster(
    adata: ad.AnnData,
    nclust: int,
    res_min: float = 0,
    res_max: float = 1,
    use_rep: str = "X_umap_stads_cpd",
    key_added: str = "cluster_stads",
    neighbors_key_added: str = "stads_cpd",
    max_iter: int = 100,
    flavor: str = "igraph",
    max_comm_size: int = 1000,
    seed: float = 0,
    balanced: bool = True,
):
    """
    Run leiden until reach the number of clusters
    """
    data = adata.copy()
    # Find neighbours
    sc.pp.neighbors(
        data,
        n_neighbors=15,
        use_rep=use_rep,
        key_added=neighbors_key_added,
        random_state=seed,
    )
    # Record a map of n vs resolution
    resolution_records = {}
    i = 0
    for target_n in range(nclust, nclust + 5):
        recorded_n = min(
            (int(k) for k in resolution_records.keys() if int(k) >= target_n),
            default=None,
        )
        if recorded_n is None:
            cur_res_min = res_min
            cur_res_max = res_max
        else:
            cur_res_min = res_min
            cur_res_max = (resolution_records[str(recorded_n)] - res_min) * 2
        print(f"Finding leiden cluster for nclust = {target_n}")
        while i < max_iter:
            res = (cur_res_min + cur_res_max) / 2
            if flavor == "igraph":
                sc.tl.leiden(
                    data,
                    resolution=res,
                    neighbors_key=neighbors_key_added,
                    flavor="igraph",
                    key_added=key_added,
                    random_state=seed,
                )
            elif flavor == "leidenalg":
                sc.tl.leiden(
                    data,
                    resolution=res,
                    neighbors_key=neighbors_key_added,
                    flavor="leidenalg",
                    key_added=key_added,
                    random_state=seed,
                    max_comm_size=max_comm_size,
                )
            n = data.obs[key_added].cat.categories.shape[0]
            print(f"Iteration {i}, resolution {res}, n = {n}")
            resolution_records[str(n)] = res
            i += 1
            if n == target_n:
                break
            elif np.isclose(cur_res_min, cur_res_max):
                raise Exception(f"Cannot find a partition for nclust={nclust}")
            elif n < target_n:
                cur_res_min = res
            elif n > target_n:
                cur_res_max = res
        data.obs[key_added] = (
            data.obs[key_added]
            .astype("category")
            .cat.rename_categories([f"C{i}" for i in range(target_n)])
        )
        if balanced:
            is_balanced, data = balanceClusters(
                data, use_rep, key_added, nclust, target_n
            )
            if is_balanced:
                break
        else:
            break

    return data
