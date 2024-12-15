import anndata as ad
import scanpy as sc

from pathlib import Path
from typing import *


def preprocessBatchData(data: ad.AnnData):
    data.var_names_make_unique()
    # Quality control
    n_cell = data.shape[0]
    sc.pp.filter_cells(data, min_genes=3)
    n_cell_qc = data.shape[0]
    print(f"{n_cell_qc}/{n_cell} cells remained after quality control")
    return data


def categorizeMeta(data: ad.AnnData, columns: List[str]):
    data.obs[columns] = data.obs[columns].astype("category")
    return data


def addMeta(
    data: ad.AnnData,
    batch: Optional[str] = None,
    patient: Optional[str] = None,
    condition: Optional[str] = None,
):
    # Convert meta column type to int to prevent save error
    numeric_cols = ["in_tissue", "array_row", "array_col", "n_genes"]
    for col in numeric_cols:
        if col in data.obs.columns:
            data.obs[col] = data.obs[col].astype("int")
    if "spatial" in data.obsm:
        data.obsm["spatial"] = data.obsm["spatial"].astype("int")
    # Assign meta info
    if batch is not None:
        data.obs["batch"] = batch
    if patient is not None:
        data.obs["patient"] = patient
    if condition is not None:
        data.obs["condition"] = condition
    # Assign spatial image key the same as batch name, so that image name will
    # keep the same after merge
    if "spatial" in data.uns and batch is not None and len(data.uns["spatial"]) > 0:
        data.uns["spatial"][batch] = data.uns["spatial"].pop(
            list(data.uns["spatial"].keys())[0]
        )
    return data


def loadData(
    data_name: str = "hcc",
    use_batch: Optional[str] = None,
    preprocess: Optional[Callable | str] = "default",
    return_type: str = "concat",
):
    adata = None
    adatas = {}

    if not preprocess:
        preprocess = lambda x: x
    elif preprocess == "default":
        preprocess = preprocessBatchData

    if data_name == "hcc":
        data_dir = "/nfs/dcmb-lgarmire/shared/public/PMC8683021"
        patients = ["HCC01", "HCC02", "HCC03", "HCC04"]
        conditions = ["N", "T"]
        if use_batch is not None:
            patient = use_batch[:5]
            condition = use_batch[-1:]
            data = sc.read_visium(Path(data_dir, patient, condition))
            data = addMeta(data, use_batch, patient, condition)
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                for condition in conditions:
                    batch = patient + condition
                    data = sc.read_visium(Path(data_dir, patient, condition))
                    data = addMeta(data, batch, patient, condition)
                    data = preprocess(data)
                    adatas[batch] = data
    elif data_name == "hcc_sim":
        data_dir = "/nfs/dcmb-lgarmire/shared/public/PMC8683021/simulation"
        patients = ["HCC01", "HCC02", "HCC03", "HCC04"]
        conditions = ["N", "T"]
        if use_batch is not None:
            patient = use_batch[:5]
            condition = use_batch[-1:]
            data = sc.read_h5ad(Path(data_dir, f"{use_batch}.h5ad"))
            data = addMeta(data, batch, patient, condition)
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                for condition in conditions:
                    batch = patient + condition
                    data = sc.read_h5ad(Path(data_dir, f"{batch}.h5ad"))
                    data = addMeta(data, batch, patient, condition)
                    data = preprocess(data)
                    adatas[batch] = data
    elif data_name == "pancreas":
        data_dir = "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/pancreas/GSE272362"
        patients = ["PT_10", "PT_11", "PT_2"]
        conditions = ["N", "T"]
        if use_batch is not None:
            patient = use_batch[:5]
            condition = use_batch[-1:]
            data = sc.read_visium(Path(data_dir, use_batch))
            data = addMeta(data, use_batch, patient, condition)
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                for condition in conditions:
                    batch = patient + "_" + condition
                    data = sc.read_visium(Path(data_dir, batch))
                    data = addMeta(data, batch, patient, condition)
                    data = preprocess(data)
                    adatas[batch] = data
    elif data_name == "nsclc":
        data_dir = "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/nsclc"
        patients = ["P10", "P11", "P15", "P16", "P17", "P19", "P24", "P25"]
        conditions = ["N", "T"]
        batch_map = {
            "P10": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                    "T3",
                    "T4",
                ],
            },
            "P11": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                    "T3",
                    "T4",
                ],
            },
            "P15": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P16": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P17": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P19": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P24": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P25": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
        }
        if use_batch is not None:
            patient = use_batch[:3]
            condition = "N" if use_batch[-2:] in batch_map[patient]["N"] else "T"
            data = sc.read_visium(Path(data_dir, use_batch))
            data = addMeta(data, use_batch, patient, condition)
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                for condition in conditions:
                    for batch in batch_map[patient][condition]:
                        batch = f"{patient}_{batch}"
                        data = sc.read_visium(Path(data_dir, batch))
                        data = addMeta(data, batch, patient, condition)
                        data = preprocess(data)
                        adatas[batch] = data
    elif data_name == "luad":
        data_dir = "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/luad"
        patients = ["P10", "P15", "P16", "P24", "P25"]
        conditions = ["N", "T"]
        batch_map = {
            "P10": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                    "T3",
                    "T4",
                ],
            },
            "P15": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P16": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P24": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P25": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
        }
        if use_batch is not None:
            patient = use_batch[:3]
            condition = "N" if use_batch[-2:] in batch_map[patient]["N"] else "T"
            data = sc.read_visium(Path(data_dir, use_batch))
            data = addMeta(data, use_batch, patient, condition)
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                for condition in conditions:
                    for batch in batch_map[patient][condition]:
                        batch = f"{patient}_{batch}"
                        data = sc.read_visium(Path(data_dir, batch))
                        data = addMeta(data, batch, patient, condition)
                        data = preprocess(data)
                        adatas[batch] = data
    elif data_name == "lusc":
        data_dir = "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/lusc"
        patients = ["P11", "P17", "P19"]
        conditions = ["N", "T"]
        batch_map = {
            "P11": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                    "T3",
                    "T4",
                ],
            },
            "P17": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
            "P19": {
                "N": [
                    "B1",
                    "B2",
                ],
                "T": [
                    "T1",
                    "T2",
                ],
            },
        }
        if use_batch is not None:
            patient = use_batch[:3]
            condition = "N" if use_batch[-2:] in batch_map[patient]["N"] else "T"
            data = sc.read_visium(Path(data_dir, use_batch))
            data = addMeta(data, use_batch, patient, condition)
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                for condition in conditions:
                    for batch in batch_map[patient][condition]:
                        batch = f"{patient}_{batch}"
                        data = sc.read_visium(Path(data_dir, batch))
                        data = addMeta(data, batch, patient, condition)
                        data = preprocess(data)
                        adatas[batch] = data
    elif data_name == "prostate":
        data_dir = "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/prostate/erickson"
        patients = ["patient1", "patient2"]
        conditions = ["N", "T"]
        batch_map = {
            "patient1": {
                "N": [
                    "V1_1",
                    "V1_2",
                ],
                "T": [
                    "H1_1",
                    "H1_2",
                    "H1_4",
                    "H1_5",
                    "H2_1",
                    "H2_2",
                    "H2_5",
                ],
            },
            "patient2": {
                "N": [
                    "H3_2",
                    "H3_4",
                    "H3_5",
                    "V1_1",
                    "V1_2",
                    "V1_3",
                    "V1_4",
                    "V1_5",
                    "V1_6",
                    "V2_1",
                    "V2_2",
                ],
                "T": [
                    "H2_1",
                    "H2_2",
                    "H3_1",
                    "H3_6",
                ],
            },
        }
        if use_batch is not None:
            patient = use_batch[:8]
            condition = "N" if use_batch[-4:] in batch_map[patient]["N"] else "T"
            data = sc.read_visium(Path(data_dir, use_batch))
            data = addMeta(data, use_batch, patient, condition)
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                for condition in conditions:
                    for batch in batch_map[patient][condition]:
                        batch = f"{patient}_{batch}"
                        data = sc.read_visium(Path(data_dir, batch))
                        data = addMeta(data, batch, patient, condition)
                        data = preprocess(data)
                        adatas[batch] = data
    elif data_name == "prostate_hirz":
        data_dir = "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/prostate/hirz"
        patients = ["patient1", "patient2"]
        if use_batch is not None:
            patient = use_batch[:8]
            data = sc.read_h5ad(Path(data_dir, f"{patient}.h5ad"))
            data = preprocess(data)
            adata = data
        else:
            for patient in patients:
                data = sc.read_h5ad(Path(data_dir, f"{patient}.h5ad"))
                data = preprocess(data)
                adatas[patient] = data

    if use_batch is None:
        adata = ad.concat(adatas, uns_merge="unique")
        adata.obs_names_make_unique()
        adata = categorizeMeta(adata, ["patient", "condition", "batch"])
        if return_type == "concat":
            return adata
        elif return_type == "batch":
            return adatas
        elif return_type == "all":
            return adata, adatas
    else:
        return adata
