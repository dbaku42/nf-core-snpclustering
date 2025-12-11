# SNP clustering pipeline (Nextflow)

Pipeline for SNP clustering with:
- VCF filtering and LD pruning
- Genotype matrix extraction
- Clustering (KMeans / DBSCAN)
- PCA / t-SNE / UMAP visualizations
- Cluster metrics (silhouette, Calinski-Harabasz, Dunn)

## Basic usage

```bash
nextflow run main.nf \
  --vcf path/to/input.vcf.gz \
  --outdir path/to/results \
  --n_clusters 3 \
  --na_filter 0.20 \
  --pca_target_var 0.7 \
  --umap_n_neighbors 25 \
  --umap_metric manhattan
