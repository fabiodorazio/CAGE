# CAGE

Utilities for **CAGE-seq** (Cap Analysis of Gene Expression) analysis, including small helper scripts in R and a (work-in-progress) Nextflow workflow scaffold. The repository is organized around the folders `Nextflow/`, `annotations/`, `bin/`, `data/`, and `src/`. :contentReference[oaicite:0]{index=0}

---

## What’s in here

### R helpers (in `src/`)
A few lightweight R utilities for working with CAGE-derived objects / outputs:

- **Annotate CAGE tag clusters** (promoter-focused) via `ChIPseeker::annotatePeak()`  
  Function: `annotate.cage.peaks()` :contentReference[oaicite:1]{index=1}
- **Subset/merge “early vs late” (maternal/zygotic) CAGE clusters** and compute strand-aware dominant-CTSS shift  
  Function: `subset.cage.object()` :contentReference[oaicite:2]{index=2}
- **Build promoter windows centered on maternal or zygotic dominant CTSS** (for motif/sequence analyses, etc.)  
  Function: `return.pattern.distr(method = "maternal" | "zygotic")` :contentReference[oaicite:3]{index=3}
- **Convert `.ctss` to BED-like output**  
  Function: `ctss.to.bed()` :contentReference[oaicite:4]{index=4}
- **Extract “sharp” vs “broad” promoters** based on interquantile width, returning resized `GRanges`  
  Function: `.extract.sharp.promoters(method = "sharp" | "broad")` :contentReference[oaicite:5]{index=5}

### Nextflow (in `Nextflow/`)
There is a `Nextflow/` directory in the repository root intended for pipeline code/configs. :contentReference[oaicite:6]{index=6}  
*(I couldn’t reliably fetch the file listing from GitHub’s folder view in this environment, so the README below keeps Nextflow instructions general and focuses on the R helpers that are accessible.)*

---

## Requirements

### For the R helpers
You’ll typically need:

- R (>= 4.x recommended)
- Bioconductor packages commonly used in this code:
  - **GenomicRanges**
  - **ChIPseeker** (used by `annotate.cage.peaks()`) :contentReference[oaicite:7]{index=7}
  - A `TxDb` object and an organism annotation database (example in code uses `org.Dr.eg.db`) :contentReference[oaicite:8]{index=8}

> Note: Some objects referenced in functions (e.g. `txdb`, `Drerio`, and `toGRanges()`) are assumed to exist in your analysis environment and/or come from your existing workflow. :contentReference[oaicite:9]{index=9}

### For Nextflow (optional)
If you plan to run the workflow code in `Nextflow/`, you’ll need a working Nextflow installation.

---

## Installation

This repo is not packaged as an R package. Clone and source what you need:

```bash
git clone https://github.com/fabiodorazio/CAGE.git
cd CAGE
