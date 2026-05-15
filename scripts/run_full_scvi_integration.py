# ============================================================
# CRC MetaThermo full-scale scVI integration
# Machine A / local high-RAM workstation
# ============================================================

import os
import re
import gc
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
import torch
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.sparse import csr_matrix
from tqdm import tqdm

BASE = r"D:\CRC_META"
OUT = r"D:\CRC_META_FULL_SCVI"
os.makedirs(OUT, exist_ok=True)

CORE_DATASETS = ["GSE178318", "GSE231559", "GSE298084"]
SUPPORT_DATASETS = ["GSE299737"]

USE_DATASETS = CORE_DATASETS
MAX_EPOCHS = 80
N_LATENT = 30

print("Versions")
print("scanpy:", sc.__version__)
print("scvi:", scvi.__version__)
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())

def infer_label_site(dataset, sample_id):
    s = sample_id.lower()

    if dataset == "GSE178318":
        if "pbmc" in s:
            return "blood", "blood"
        if "_lm" in s or "lm" in s:
            return "metastasis", "liver"
        if "_crc" in s or "crc" in s:
            return "primary", "colon"

    if dataset == "GSE298084":
        if "pbmc" in s or "blood" in s:
            return "blood", "blood"
        if "liver" in s:
            return "metastasis", "liver"
        if "colon" in s:
            return "primary", "colon"

    if dataset == "GSE231559":
        normal_liver = [
            "GSM7290760", "GSM7290764", "GSM7290765", "GSM7290766",
            "GSM7290770", "GSM7290776", "GSM7290780", "GSM7290784"
        ]
        normal_colon = ["GSM7290762", "GSM7290768", "GSM7290771"]
        primary_colon = [
            "GSM7290763", "GSM7290769", "GSM7290772",
            "GSM7290773", "GSM7290774", "GSM7290777"
        ]
        metastasis_liver = [
            "GSM7290761", "GSM7290767", "GSM7290775",
            "GSM7290778", "GSM7290779", "GSM7290781",
            "GSM7290782", "GSM7290783", "GSM7290785"
        ]
        if any(x.lower() in s for x in normal_liver):
            return "normal", "liver"
        if any(x.lower() in s for x in normal_colon):
            return "normal", "colon"
        if any(x.lower() in s for x in primary_colon):
            return "primary", "colon"
        if any(x.lower() in s for x in metastasis_liver):
            return "metastasis", "liver"

    if dataset == "GSE299737":
        if "_ln" in s or "ln-" in s:
            return "metastasis", "lymph_node"
        if "_td" in s or "td-" in s:
            return "metastasis", "tumor_deposit"
        if "_ti" in s or "ti-" in s:
            return "primary", "mucosal_side"
        if "_to" in s or "to-" in s:
            return "primary", "serosal_side"
        if "_t_" in s or "_t_auto" in s:
            return "primary", "colon"

    return "unknown", "unknown"

def read_sample_folder(dataset, sample_dir):
    mtx = os.path.join(sample_dir, "matrix.mtx")
    features = os.path.join(sample_dir, "features.tsv")
    barcodes = os.path.join(sample_dir, "barcodes.tsv")
    metadata = os.path.join(sample_dir, "metadata.csv")

    if not os.path.exists(mtx):
        raise FileNotFoundError(mtx)

    X = mmread(mtx).T.tocsr()
    feat = pd.read_csv(features, sep="\t", header=None)
    bc = pd.read_csv(barcodes, sep="\t", header=None)

    if feat.shape[1] >= 2:
        gene_ids = feat.iloc[:, 0].astype(str).values
        gene_names = feat.iloc[:, 1].astype(str).values
    else:
        gene_ids = feat.iloc[:, 0].astype(str).values
        gene_names = feat.iloc[:, 0].astype(str).values

    sample_id = os.path.basename(sample_dir)

    obs = pd.DataFrame(index=bc.iloc[:, 0].astype(str).values)

    if os.path.exists(metadata):
        meta = pd.read_csv(metadata)
        if "barcode" in meta.columns:
            meta.index = meta["barcode"].astype(str)
            obs = obs.join(meta, how="left")

    obs["dataset"] = dataset
    obs["sample_id"] = sample_id
    obs["batch"] = dataset
    obs["sample_batch"] = dataset + "__" + sample_id

    label, site = infer_label_site(dataset, sample_id)
    obs["label"] = obs["label"].fillna(label) if "label" in obs.columns else label
    obs["site"] = obs["site"].fillna(site) if "site" in obs.columns else site

    var = pd.DataFrame(index=gene_names)
    var["gene_ids"] = gene_ids

    a = ad.AnnData(X=X, obs=obs, var=var)
    a.obs_names_make_unique()
    a.var_names_make_unique()
    return a

# ============================================================
# Stage 1: Build full core AnnData
# ============================================================

raw_path = os.path.join(OUT, "01_CRC_CORE_raw_counts.h5ad")

if os.path.exists(raw_path):
    print("Loading existing raw AnnData:", raw_path)
    adata = sc.read_h5ad(raw_path)
else:
    adatas = []
    log_rows = []

    for dataset in USE_DATASETS:
        ds_dir = os.path.join(BASE, dataset)
        sample_dirs = [
            os.path.join(ds_dir, x)
            for x in os.listdir(ds_dir)
            if os.path.isdir(os.path.join(ds_dir, x))
        ]

        print(f"Dataset {dataset}: {len(sample_dirs)} samples")

        for sdir in tqdm(sample_dirs, desc=f"Reading {dataset}"):
            try:
                a = read_sample_folder(dataset, sdir)
                adatas.append(a)

                log_rows.append({
                    "dataset": dataset,
                    "sample_id": os.path.basename(sdir),
                    "cells": a.n_obs,
                    "genes": a.n_vars,
                    "label": str(a.obs["label"].iloc[0]),
                    "site": str(a.obs["site"].iloc[0]),
                    "status": "ok",
                    "error": ""
                })
            except Exception as e:
                log_rows.append({
                    "dataset": dataset,
                    "sample_id": os.path.basename(sdir),
                    "cells": None,
                    "genes": None,
                    "label": "",
                    "site": "",
                    "status": "failed",
                    "error": str(e)
                })

    pd.DataFrame(log_rows).to_csv(os.path.join(OUT, "01_read_sample_log.csv"), index=False)

    print("Concatenating samples...")
    adata = ad.concat(
        adatas,
        join="inner",
        index_unique="-",
        label="concat_batch",
        keys=[f"{a.obs['dataset'].iloc[0]}__{a.obs['sample_id'].iloc[0]}" for a in adatas]
    )

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    print(adata)
    adata.write(raw_path)

    del adatas
    gc.collect()

# ============================================================
# Stage 2: QC summaries
# ============================================================

print("Core object:", adata)

pd.crosstab(adata.obs["dataset"], adata.obs["label"]).to_csv(
    os.path.join(OUT, "02_core_dataset_by_label_cells.csv")
)

pd.crosstab(adata.obs["dataset"], adata.obs["site"]).to_csv(
    os.path.join(OUT, "02_core_dataset_by_site_cells.csv")
)

sample_summary = (
    adata.obs
    .groupby(["dataset", "sample_id", "label", "site"], observed=True)
    .size()
    .reset_index(name="cells")
)
sample_summary.to_csv(os.path.join(OUT, "02_core_sample_summary.csv"), index=False)

# ============================================================
# Stage 3: HVG for scVI
# ============================================================

hvg_path = os.path.join(OUT, "03_CRC_CORE_hvg3000_counts.h5ad")

if os.path.exists(hvg_path):
    print("Loading existing HVG object:", hvg_path)
    adata_hvg = sc.read_h5ad(hvg_path)
else:
    print("Selecting HVGs...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=3000,
        flavor="seurat_v3",
        batch_key="dataset",
        subset=False
    )

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    adata_hvg.write(hvg_path)

# ============================================================
# Stage 4: scVI model
# ============================================================

print("Setting up scVI...")
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    batch_key="dataset",
    categorical_covariate_keys=["sample_id", "label", "site"]
)

model_dir = os.path.join(OUT, "04_scvi_model")

if os.path.exists(model_dir):
    print("Loading existing scVI model...")
    model = scvi.model.SCVI.load(model_dir, adata=adata_hvg)
else:
    model = scvi.model.SCVI(
        adata_hvg,
        n_latent=N_LATENT,
        n_layers=2,
        n_hidden=128,
        gene_likelihood="nb"
    )

    print("Training scVI...")
    model.train(max_epochs=MAX_EPOCHS, accelerator="auto")
    model.save(model_dir, overwrite=True)

# ============================================================
# Stage 5: latent space, UMAP, clustering
# ============================================================

print("Getting latent representation...")
adata_hvg.obsm["X_scVI"] = model.get_latent_representation()

sc.pp.neighbors(adata_hvg, use_rep="X_scVI", n_neighbors=20)
sc.tl.umap(adata_hvg, min_dist=0.3)
sc.tl.leiden(adata_hvg, resolution=0.5, key_added="cluster")

integrated_path = os.path.join(OUT, "05_CRC_CORE_scVI_integrated.h5ad")
adata_hvg.write(integrated_path)

# ============================================================
# Stage 6: export tables and figures
# ============================================================

umap = pd.DataFrame(
    adata_hvg.obsm["X_umap"],
    columns=["UMAP_1", "UMAP_2"],
    index=adata_hvg.obs_names
)
umap_meta = pd.concat([umap, adata_hvg.obs.copy()], axis=1)
umap_meta.to_csv(os.path.join(OUT, "06_scVI_umap_coordinates_metadata.csv"))

pd.crosstab(adata_hvg.obs["dataset"], adata_hvg.obs["label"]).to_csv(
    os.path.join(OUT, "06_dataset_by_label.csv")
)
pd.crosstab(adata_hvg.obs["cluster"], adata_hvg.obs["label"]).to_csv(
    os.path.join(OUT, "06_cluster_by_label.csv")
)
pd.crosstab(adata_hvg.obs["cluster"], adata_hvg.obs["dataset"]).to_csv(
    os.path.join(OUT, "06_cluster_by_dataset.csv")
)

sc.settings.figdir = OUT

for color, fname, title in [
    ("dataset", "Figure3A_scVI_UMAP_by_dataset.tiff", "scVI integrated UMAP by dataset"),
    ("label", "Figure3B_scVI_UMAP_by_label.tiff", "scVI integrated UMAP by label"),
    ("site", "Figure3C_scVI_UMAP_by_site.tiff", "scVI integrated UMAP by site"),
    ("cluster", "Figure3D_scVI_UMAP_by_cluster.tiff", "scVI integrated UMAP by cluster"),
]:
    sc.pl.umap(adata_hvg, color=color, frameon=False, show=False, size=2, title=title)
    plt.savefig(os.path.join(OUT, fname), dpi=300, bbox_inches="tight")
    plt.close()

report = {
    "raw_object": str(adata),
    "hvg_object": str(adata_hvg),
    "datasets": USE_DATASETS,
    "max_epochs": MAX_EPOCHS,
    "n_latent": N_LATENT,
    "output_dir": OUT,
    "cuda_available": torch.cuda.is_available()
}

with open(os.path.join(OUT, "00_run_report.json"), "w") as f:
    json.dump(report, f, indent=2)

print("DONE")
print("Output:", OUT)
print("Integrated object:", integrated_path)
print("Dataset by label:")
print(pd.crosstab(adata_hvg.obs["dataset"], adata_hvg.obs["label"]))