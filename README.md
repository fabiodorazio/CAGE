# CAGE

All **CAGE-seq** (Cap Analysis of Gene Expression) analysis, including small helper scripts a Nextflow workflow scaffold.

---

## What’s in here

### R helpers (in `src/`)
A few lightweight R utilities for working with CAGE-derived objects / outputs:

- **Annotate CAGE tag clusters** (promoter-focused) via `CAGEr`  
  Function: `annotate.cage.peaks()`[](bin/)
- **Subset/merge “early vs late” (maternal/zygotic) CAGE clusters** and compute strand-aware dominant-CTSS shift  
  Function: `subset.cage.object()`
- **Build promoter windows centered on maternal or zygotic dominant CTSS** (for motif/sequence analyses, etc.)  
  Function: `return.pattern.distr(method = "maternal" | "zygotic")`
- **Convert `.ctss` to BED-like output**  
  Function: `ctss.to.bed()` :
- **Extract “sharp” vs “broad” promoters** based on interquantile width, returning resized `GRanges`  
  Function: `.extract.sharp.promoters(method = "sharp" | "broad")`

### Nextflow (in `Nextflow/`)
The CAGE pipeline can be run through Nextflow

---

## Requirements

### For the R helpers

- R (>= 4.x recommended)
- Bioconductor packages commonly used in this code:
  - **GenomicRanges**
  - **ChIPseeker** (used by `annotate.cage.peaks()`)
  - A `TxDb` object and an organism annotation database (example in code uses `org.Dr.eg.db`)


### For Nextflow
If you plan to run the workflow code in `Nextflow/`, you’ll need a working Nextflow installation.
(DSL2 enabled; scripts set `nextflow.enable.dsl=2`).

### Tooling (via conda env provided)
`Nextflow/environment.yml` pins the core tools used by the workflow, including:
- nextflow 20.10.0
- fastqc, multiqc
- star (2.6.1d), bowtie (1.2.3)
- samtools, bedtools, cutadapt, sortmerna

The workflow is a DSL2 scaffold with modular steps. From the module names and scripts currently present, it supports:
1. **Input handling**
   - Read FASTQs either from a glob pattern (`params.input`) or from a `--samplesheet` CSV
   - Optional samplesheet validation using a `check_samplesheet.py` wrapper process

2. **QC**
   - **FastQC** on input reads, written to `${outdir}/fastqc/` (zip files grouped under `zips/`)

3. **Reference indexing**
   - Build a **STAR** genome index from `--fasta` and `--gtf` if no STAR index is provided
   - Build a **bowtie1** index from `--fasta` if no bowtie index is provided

4. **Alignment**
   - Align reads with **STAR** (sorted BAM output) into `${outdir}/STAR/`
   - Or align reads with **bowtie1** (BAM output + logs) into `${outdir}/bowtie/`


---

## Installation

This repo is not packaged as an R package. Clone and source what you need:

```bash
git clone https://github.com/fabiodorazio/CAGE.git
cd CAGE
