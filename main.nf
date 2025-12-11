nextflow.enable.dsl = 2

/*
 * PARAMETERS
 */
params.vcf               = params.vcf               ?: null
params.outdir            = params.outdir            ?: './results'

params.maf               = params.maf               ?: 0.01
params.missing           = params.missing           ?: 0.1
params.na_filter        = params.na_filter        ?: 0.10
params.threads           = params.threads           ?: 4
params.container         = params.container         ?: ''
params.container_py      = params.container_py      ?: ''

params.n_clusters        = params.n_clusters        ?: 3
params.cluster_method    = params.cluster_method    ?: 'kmeans'
params.db_eps            = params.db_eps            ?: 0.5
params.db_min_samples    = params.db_min_samples    ?: 5

params.pca_target_var    = params.pca_target_var    ?: 0.8
params.tsne_perplexity   = params.tsne_perplexity   ?: 30.0
params.tsne_n_iter       = params.tsne_n_iter       ?: 1000
params.umap_n_neighbors  = params.umap_n_neighbors  ?: 15
params.umap_min_dist     = params.umap_min_dist     ?: 0.1
params.umap_n_components = params.umap_n_components ?: 2

/*
 * MODULE INCLUDES
 */
include { BCFTOOLS_FILTER } from './modules/preprocess.nf'
include { LD_PRUNE        } from './modules/ld_prune.nf'
include { EXTRACT_MATRIX  } from './modules/extract_matrix.nf'
include { CLUSTERING      } from './modules/clustering.nf'
include { CLUSTER_METRICS } from './modules/cluster_metrics.nf'
include { REPORT          } from './modules/report.nf'


/*
 * MAIN WORKFLOW
 */
workflow {

    if( !params.vcf ) {
        error "Missing --vcf"
    }

    Channel
        .fromPath(params.vcf)
        .set { vcf_ch }

    /*
     * 1) PREPROCESS (bcftools: filter by MAF and missingness)
     *    Output: filtered VCF (.filtered.vcf.gz)
     */
    filtered_vcf_ch = BCFTOOLS_FILTER(vcf_ch)

    /*
     * 2) LD PRUNING (plink2)
     *    Input: filtered VCF
     *    Output: pruned VCF (.plink.pruned.vcf)
     */
    pruned_vcf_ch = LD_PRUNE(filtered_vcf_ch)

    /*
     * 3) GENOTYPE MATRIX EXTRACTION
     *    Input: pruned VCF
     *    Output: TSV genotype matrix (sample + var_*)
     */
    matrix_ch = EXTRACT_MATRIX(pruned_vcf_ch)

    /*
     * 4) CLUSTERING
     *    Input: genotype TSV
     *    Output: tuple (id, matrix_path, clusters_csv)
     */
    cluster_out_ch = CLUSTERING(matrix_ch)

    /*
     * 5) CLUSTERING METRICS
     */
    metrics_out_ch = CLUSTER_METRICS(cluster_out_ch)

    /*
     * 6) REPORT (PCA / t-SNE / UMAP colored by cluster)
     */
    REPORT(cluster_out_ch)

    /*
     * No emit from main workflow: files are published to params.outdir via publishDir.
     */
}
