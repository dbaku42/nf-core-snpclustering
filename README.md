# SNP clustering pipeline (Nextflow)

This repository contains a Nextflow pipeline for unsupervised clustering of SNP genotype data, with dimensionality reduction and rich visualisation.

The workflow is designed for VCF input (e.g. 1000 Genomes data) and performs:

- VCF quality / missingness filtering (`bcftools`)
- LD pruning (`plink2`)
- Genotype matrix extraction (samples × variants)
- Handling of missing values (`NA` filtering and optional encoding as `-1`)
- Clustering (K-means / DBSCAN)
- Dimensionality reduction (PCA / t-SNE / UMAP)
- Cluster quality metrics (silhouette, Calinski–Harabasz, Dunn)
- Graphical reports (PCA / t-SNE / UMAP coloured by cluster)

---

## 1. Pipeline overview

Main steps:

1. **Preprocessing (BCFTOOLS_FILTER)**  
   - Input: compressed VCF (`.vcf.gz`).
   - Filters variants by minor allele frequency and missingness using `bcftools view`.
   - Output: filtered VCF (`*.filtered.vcf.gz`).

2. **LD pruning (LD_PRUNE)**  
   - Converts filtered VCF to PLINK2 PGEN format.  
   - Performs LD pruning (e.g. `--indep-pairwise 50 5 0.2`).  
   - Exports a pruned VCF containing approximately independent SNPs.

3. **Genotype matrix extraction (EXTRACT_MATRIX)**  
   - Queries genotypes with `bcftools query`.  
   - Converts genotypes (e.g. `0/0`, `0/1`, `1/1`, phased `0|1`, etc.) to **dosage**:
     - 0 → homozygous reference  
     - 1 → heterozygous  
     - 2 → homozygous alternate  
   - Missing calls can be filtered out (`--na_filter`) or encoded as `-1` for clustering.

4. **Clustering (CLUSTERING)**  
   - Reads the sample × variant dosage matrix.
   - Reduce dimensionality with PCA
   - Supports **K-means** (default) and optional **DBSCAN**.  
   - Allows tuning of:
     - Number of clusters (`--n_clusters`)
     - DBSCAN parameters (`--db_eps`, `--db_min_samples`)
   - Writes:
     - Cluster labels per sample (`*_clusters.csv`)
     - Optionally intermediate embeddings (if used inside clustering scripts).

5. **Cluster metrics (CLUSTER_METRICS)**  
   - Computes standard clustering indices for multiple `k` (e.g. 2–10):
     - Silhouette score  
     - Calinski–Harabasz index  
     - Dunn index  
   - Generates metric plots vs. `k` (e.g. elbow-style curves) to help choose the number of clusters.

6. **Report / visualisation (REPORT)**  
   - Uses the final genotype matrix and cluster labels to compute:
     - **PCA** (target explained variance)  
     - **t-SNE** embedding  
     - **UMAP** embedding (with configurable metric and neighbours)  
   - Produces:
     - `*_pca_clusters.png`  
     - `*_tsne_clusters.png`  
     - `*_umap_clusters.png`  
     - `*_combined_embeddings.png` (PCA / t-SNE / UMAP side by side)

---

## 2. Requirements

- **Nextflow** (>= 22.x recommended)
- **Java** (for Nextflow)

### Tools
- `bcftools`
- `plink2`

### Python (for clustering / plots / report)
Python 3 with:
- `numpy`
- `pandas`
- `matplotlib`
- `scikit-learn`
- `umap-learn`

You can manage the Python environment with `conda`, `mamba`, or `venv`.

### Optional: FlashPCA2 (containerized)
If you choose `--pca_engine flashpca`, the pipeline will run PCA using **FlashPCA2**.
FlashPCA2 is provided via a container (recommended), so you do **not** need to install it on the host system.

To run FlashPCA2 you need **one** of:
- Docker, or
- Apptainer/Singularity (HPC-friendly)


---

## 3. FlashPCA2 via container (optional)

The pipeline supports two PCA engines:
- `--pca_engine sklearn` (default, runs with Python/scikit-learn)
- `--pca_engine flashpca` (runs FlashPCA2 via container)

When `--pca_engine flashpca` is selected, the pipeline will automatically use the FlashPCA2 container.

### 3.1 Use a prebuilt container image (recommended)

Run with Docker:
```bash
nextflow run main.nf \
  --vcf path/to/input.vcf.gz \
  --outdir results \
  --pca_engine flashpca \
  --flashpca_container ghcr.io/dbaku42/flashpca2:0.1 \
  --flashpca_bin flashpca \
  --n_pcs \
  --algorithm kmeans|dbscan \
  -with-docker

