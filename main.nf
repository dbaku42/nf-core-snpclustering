nextflow.enable.dsl=2

/*
Pipeline: VCF
  01 preprocess        (filter MAF + missing)
  02 pruning           (LD pruning)
  03 extraction_matrix (PLINK .raw genotype matrix)   [solo per PCA sklearn]
  03b plink1_convert   (PLINK2 -> PLINK1 bed/bim/fam) [solo per flashpca2]
  04 pca               (sklearn OR flashpca2)
  05 clustering        (KMeans or DBSCAN)
  06 cluster_metrics   (elbow + silhouette etc.)
  07 cluster_viz       (UMAP or t-SNE)
  08 report            (simple HTML report)
*/

params.vcf            = params.vcf            ?: '/mnt/data/example.vcf.gz'
params.outdir         = params.outdir         ?: './results'
params.threads        = params.threads        ?: 4

// preprocess
params.maf            = params.maf            ?: 0.01
params.missing        = params.missing        ?: 0.10

// pruning (PLINK2 indep-pairwise)
params.ld_window_kb   = params.ld_window_kb   ?: 50
params.ld_step        = params.ld_step        ?: 5
params.ld_r2          = params.ld_r2          ?: 0.2
params.threads = params.threads ?: 4
params.py_container = params.py_container ?: 'python:3.11-slim'
params.plink2_container = params.plink2_container ?: 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
params.bcftools_container = params.bcftools_container ?: 'quay.io/biocontainers/bcftools:1.17--h9ee0642_0'
params.flashpca_container = params.flashpca_container ?: 'flashpca2:latest'   // tua immagine
params.flashpca_bin = params.flashpca_bin ?: '/home/donald/flashpca/flashpca'

// PCA
params.n_pcs          = params.n_pcs          ?: 40
params.ipca_batch_size = params.ipca_batch_size ?: 256
params.imputer        = params.imputer        ?: 'minus1'

// NEW: scegli backend PCA
//   'sklearn'   -> usa extraction_matrix (.raw) + pca.py (come adesso)
//   'flashpca2' -> usa plink1_convert (bed/bim/fam) + flashpca2
params.pca_engine     = params.pca_engine     ?: 'sklearn'   // sklearn|flashpca2

// clustering
params.k_min = params.k_min ?: 2
params.dbscan_eps = params.dbscan_eps ?: 0.5
params.dbscan_min_samples = params.dbscan_min_samples ?: 5
params.viz_perplexity = params.viz_perplexity ?: 30
params.viz_tsne_iter = params.viz_tsne_iter ?: 1000
params.viz_umap_neighbors = params.viz_umap_neighbors ?: 15
params.viz_umap_min_dist = params.viz_umap_min_dist ?: 0.1

params.algorithm      = params.algorithm      ?: 'kmeans'   // kmeans|dbscan
params.n_clusters     = params.n_clusters     ?: 3
params.k_min          = params.k_min          ?: 1
params.k_max          = params.k_max          ?: 12
params.dbscan_eps     = params.dbscan_eps     ?: 0.5
params.dbscan_min_samples = params.dbscan_min_samples ?: 5

// viz
params.viz_perplexity     = params.viz_perplexity     ?: 30     // t-SNE
params.viz_tsne_iter      = params.viz_tsne_iter      ?: 1000
params.viz_umap_neighbors = params.viz_umap_neighbors ?: 15
params.viz_umap_min_dist  = params.viz_umap_min_dist  ?: 0.1

// containers (override as needed)
params.bcftools_container = params.bcftools_container ?: 'quay.io/biocontainers/bcftools:1.17--h9ee0642_0'
params.plink2_container   = params.plink2_container   ?: 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
params.py_container       = params.py_container       ?: 'python:3.11-slim'

// NEW: container con flashpca2 (devi metterci flashpca nel PATH)
params.flashpca_container = params.flashpca_container ?: 'flashpca2:latest'

// Modules
include { preprocess_ch }        from './modules/preprocess.nf'
include { pruning_ch }           from './modules/pruning.nf'
include { extraction_matrix_ch } from './modules/extract_matrix.nf'
include { pca_ch }               from './modules/pca.nf'

include { clustering_ch }        from './modules/clustering.nf'
include { cluster_metrics_ch }   from './modules/cluster_metrics.nf'
include { cluster_viz_ch }       from './modules/cluster_viz.nf'
include { report_ch }            from './modules/report.nf'

// Basic validation
if( !file(params.vcf).exists() ) {
    exit 1, "VCF file not found: ${params.vcf}"
}
if( params.maf < 0 || params.maf > 0.5 ) {
    exit 1, "params.maf must be between 0 and 0.5"
}
if( params.missing < 0 || params.missing > 1 ) {
    exit 1, "params.missing must be between 0 and 1"
}
if( !['sklearn','flashpca2'].contains(params.pca_engine) ) {
    exit 1, "params.pca_engine must be 'sklearn' or 'flashpca2'"
}

workflow {

    Channel
      .fromPath(params.vcf)
      .set { vcf_ch }

    // 01
    filtered_vcf_ch = preprocess_ch(vcf_ch)

    // 02
    pruned_plink_ch = pruning_ch(filtered_vcf_ch)

    // 03/04: PCA backend switch
    pca_out_ch = pca_ch(pruned_plink_ch)


    // 05
    clust_out   = clustering_ch(pca_out_ch)
    clusters_ch = clust_out.clusters

    // 06
    metrics_out = cluster_metrics_ch(pca_out_ch, clusters_ch)

    // 07
    viz_out     = cluster_viz_ch(pca_out_ch, clusters_ch)

    // 08
    report_ch(
        pca_out_ch,
        clusters_ch,
        metrics_out.k_sweep,
        metrics_out.selected,
        viz_out.umap_embedding,
        viz_out.tsne_embedding,
        viz_out.umap_plot,
        viz_out.tsne_plot,
        metrics_out.elbow,
        metrics_out.silhouette,
        metrics_out.davies_bouldin,
        viz_out.pca_plot
    )
}
