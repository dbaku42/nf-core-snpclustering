nextflow.enable.dsl=2

// --- I/O ---
params.vcf     = params.vcf     ?: '/mnt/data/example.vcf.gz'
params.outdir  = params.outdir  ?: './results'
params.threads = params.threads ?: 4

// --- Preprocess ---
params.maf     = params.maf     ?: 0.01
params.missing = params.missing ?: 0.10

// --- LD pruning ---
params.ld_window_kb = params.ld_window_kb ?: 50
params.ld_step      = params.ld_step      ?: 5
params.ld_r2        = params.ld_r2        ?: 0.2

// --- PCA ---
params.n_pcs           = params.n_pcs           ?: 40
params.ipca_batch_size = params.ipca_batch_size ?: 256
params.imputer         = params.imputer         ?: 'minus1'
params.pca_engine      = params.pca_engine      ?: 'sklearn'

// --- Clustering ---
params.algorithm          = params.algorithm          ?: 'kmeans'
params.n_clusters         = params.n_clusters         ?: 3
params.k_min              = params.k_min              ?: 2
params.k_max              = params.k_max              ?: 12
params.dbscan_eps         = params.dbscan_eps         ?: 0.5
params.dbscan_min_samples = params.dbscan_min_samples ?: 5
params.n_init             = params.n_init             ?: 100
params.init_method        = params.init_method        ?: 'random'

// --- Viz ---
params.viz_perplexity     = params.viz_perplexity     ?: 30
params.viz_tsne_iter      = params.viz_tsne_iter      ?: 1000
params.viz_umap_neighbors = params.viz_umap_neighbors ?: 15
params.viz_umap_min_dist  = params.viz_umap_min_dist  ?: 0.1

// --- Containers ---
params.container_py       = params.container_py       ?: false
params.plink2_container   = params.plink2_container   ?: 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'
params.bcftools_container = params.bcftools_container ?: 'quay.io/biocontainers/bcftools:1.17--h9ee0642_0'
params.flashpca_container = params.flashpca_container ?: 'flashpca2:latest'
params.flashpca_bin       = params.flashpca_bin       ?: '/home/donald/flashpca/flashpca'

// --- Modules ---
include { preprocess_ch }        from './modules/preprocess.nf'
include { pruning_ch }           from './modules/pruning.nf'
include { extraction_matrix_ch } from './modules/extract_matrix.nf'
include { pca_ch }               from './modules/pca.nf'
include { clustering_ch }        from './modules/clustering.nf'
include { cluster_metrics_ch }   from './modules/cluster_metrics.nf'
include { cluster_viz_ch }       from './modules/cluster_viz.nf'
include { report_ch }            from './modules/report.nf'

// --- Validation ---
if( !file(params.vcf).exists() ) {
    exit 1, "VCF file not found: ${params.vcf}"
}

workflow {

    Channel.fromPath(params.vcf).set { vcf_ch }

    // 01-04
    filtered_vcf_ch = preprocess_ch(vcf_ch)
    pruned_plink_ch = pruning_ch(filtered_vcf_ch)
    pca_out_ch      = pca_ch(pruned_plink_ch)
    // pca_out_ch: tuple path(features.tsv), path(scaled.tsv), path(pca_scores.tsv), path(pca_info.json)

    // 05 clustering
    clust_out   = clustering_ch(pca_out_ch)
    clusters_ch = clust_out.clusters
    // clusters_ch: tuple val(id), path(features.tsv), path(clusters.csv), path(clustering_info.json)

    // Canali ausiliari taggati con id per i join
    pca_scores_ch = pca_out_ch.map { features, scaled, pca_scores, pca_info ->
        tuple(features.baseName, pca_scores)
    }
    pca_info_ch = pca_out_ch.map { features, scaled, pca_scores, pca_info ->
        tuple(features.baseName, pca_info)
    }

    // Canale arricchito con pca_scores — usato da CLUSTER_METRICS e CLUSTER_VIZ
    // [id, mat, clusters, clustering_info, pca_scores]
    clusters_with_pca_ch = clusters_ch.join(pca_scores_ch)

    // 06 cluster metrics
    metrics_out = cluster_metrics_ch(clusters_with_pca_ch)

    // 07 cluster viz
    viz_out = cluster_viz_ch(clusters_with_pca_ch)

    // 08 report HTML — join di tutti i canali per id
    report_input_ch = clusters_ch
        .join( pca_info_ch )
        .join( metrics_out.k_sweep       .map{ f -> tuple(f.baseName.replaceAll(/_k_sweep$/,         ''), f) } )
        .join( metrics_out.selected      .map{ f -> tuple(f.baseName.replaceAll(/_selected$/,        ''), f) } )
        .join( metrics_out.elbow         .map{ f -> tuple(f.baseName.replaceAll(/_elbow$/,           ''), f) } )
        .join( metrics_out.silhouette    .map{ f -> tuple(f.baseName.replaceAll(/_silhouette$/,      ''), f) } )
        .join( metrics_out.davies_bouldin.map{ f -> tuple(f.baseName.replaceAll(/_davies_bouldin$/,  ''), f) } )
        .join( metrics_out.calinski      .map{ f -> tuple(f.baseName.replaceAll(/_calinski$/,        ''), f) } )
        .join( viz_out.umap_embedding    .map{ f -> tuple(f.baseName.replaceAll(/_umap_embedding$/,  ''), f) } )
        .join( viz_out.tsne_embedding    .map{ f -> tuple(f.baseName.replaceAll(/_tsne_embedding$/,  ''), f) } )
        .join( viz_out.umap_plot         .map{ f -> tuple(f.baseName.replaceAll(/_umap_clusters$/,   ''), f) } )
        .join( viz_out.tsne_plot         .map{ f -> tuple(f.baseName.replaceAll(/_tsne_clusters$/,   ''), f) } )
        .join( viz_out.pca_plot          .map{ f -> tuple(f.baseName.replaceAll(/_pca_clusters$/,    ''), f) } )
        .map { id, mat, clusters_csv, clustering_info,
               pca_info,
               k_sweep, selected,
               elbow, sil, db, cal,
               umap_tsv, tsne_tsv, umap_png, tsne_png, pca_png ->
            tuple(id, mat, clusters_csv, clustering_info,
                  pca_info, k_sweep, selected,
                  elbow, sil, db, cal,
                  umap_tsv, tsne_tsv, umap_png, tsne_png, pca_png)
        }

    report_out = report_ch(report_input_ch)
}
