# CRC-Metastasis-SingleCell

## Multi-cohort single-cell integration reveals a stress-associated immune trajectory and constrained perturbation landscape in colorectal cancer metastasis

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![scVI](https://img.shields.io/badge/scVI-single--cell-green)]()
[![License](https://img.shields.io/badge/license-Academic-lightgrey)]()

---

## Overview

This repository contains the **minimal reproducible computational workflow** supporting our study:

> **Integrated single-cell analysis of colorectal cancer metastasis reveals a stress-associated trajectory and constrained pharmacological reversal landscape**

The study systematically integrates independent colorectal cancer (CRC) single-cell RNA sequencing datasets to construct a unified atlas of metastasis and identify transcriptional programs associated with immune remodeling, stress adaptation, and therapeutic perturbation.

The repository provides the **minimal code necessary to reproduce principal analyses and manuscript figures**. It intentionally excludes internal automation systems, API orchestration frameworks, server infrastructure, and commercial platform components.

---

## Scientific Objectives

### Objective 1 — Systematic dataset identification and qualification

- systematic GEO search
- dataset screening and metadata extraction
- API-assisted classification
- identification of integration-compatible datasets

Final selected datasets:

| Dataset | Role | Biological Context |
|----------|------|-------------------|
| GSE178318 | Core | primary tumor, liver metastasis, blood |
| GSE231559 | Core | primary tumor, liver metastasis, normal tissue |
| GSE298084 | Core | primary tumor, liver metastasis, blood |
| GSE299737 | Support | metastasis-gradient context |
| GSE110009 | Support | validation |

---

### Objective 2 — Construction of an integrated CRC single-cell atlas

- scVI-based harmonization
- batch correction
- clustering and annotation
- immune and tumor trajectory inference
- pseudotime modeling
- metastasis-associated transcriptional programs

---

### Objective 3 — Drug perturbation landscape analysis

- metastasis-associated gene signatures
- perturbation overlap mapping
- drug–gene interaction analysis
- hub gene detection
- mechanism-flow analysis

---

## Repository Structure

```text
CRC_MetaThermo_v1_code/
│
├── scripts/
│   ├── figure5*        Differential expression & pathway analysis
│   ├── figure6*        Tumor & immune trajectory analysis
│   ├── figure7*        Drug perturbation and network analysis
│   ├── figure9*        Conceptual model generation
│   └── run_full_scvi_integration.py
│
├── environment.yml
├── requirements.txt
└── README.md
```

---

## Figure Reproducibility Map

| Figure | Description | Primary Script |
|--------|-------------|----------------|
| Figure 5 | Differential expression & pathways | `figure5*.py` |
| Figure 6 | Tumor and immune trajectory | `figure6*.py` |
| Figure 7 | Therapeutic perturbation landscape | `figure7*.py` |
| Figure 9 | Conceptual model | `figure9*.R`, `figure9*.py` |

---

## Data Availability

Original datasets are publicly available through the NCBI Gene Expression Omnibus (GEO):

- GSE178318
- GSE231559
- GSE298084
- GSE299737
- GSE110009

Processed data supporting this study, including:

- integrated scVI atlas
- metadata annotations
- differential expression outputs
- trajectory statistics
- perturbation results

are archived separately in the accompanying data release.

---

## Computational Environment

The workflow was developed using:

- Python 3.10
- Scanpy
- scVI-tools
- AnnData
- NumPy
- Pandas
- Matplotlib
- NetworkX

Install dependencies:

### Conda

```bash
conda env create -f environment.yml
conda activate crc_metathermo
```

### Pip

```bash
pip install -r requirements.txt
```

---

## Minimal Reproducibility Workflow

### 1. Integration

```bash
python scripts/run_full_scvi_integration.py
```

### 2. Differential expression and pathway analysis

```bash
python scripts/figure5A_primary_vs_metastasis_DE.py
python scripts/figure5D_pathway_enrichment.py
```

### 3. Trajectory analysis

```bash
python scripts/figure6_tumor_trajectory.py
python scripts/figure6_immune_trajectory.py
```

### 4. Drug perturbation analysis

```bash
python scripts/figure7_drug_reversal_advanced.py
```

---

## Scope of Release

This repository represents a **minimal scientific reproducibility release**.

The following components are intentionally excluded:

- internal automation systems
- private API orchestration
- commercial infrastructure
- server configurations
- authentication systems
- platform-specific deployment code

---

## Citation

If using this repository, please cite:

```text
Javad Alizargar et al.
Integrated single-cell analysis of colorectal cancer metastasis reveals a stress-associated trajectory and constrained pharmacological reversal landscape.
(under review)
```

---

## Contact

For scientific questions regarding the analysis:

**Javad Alizargar**  
National Taiwan University / NTUH  
Taipei, Taiwan
