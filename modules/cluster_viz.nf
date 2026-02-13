process CLUSTER_VIZ {

    publishDir "${params.outdir}/viz", mode: 'copy'

    tag { id }

    if( params.container_py ) container { params.container_py }
    cpus params.threads
    memory '8 GB'

    input:
      tuple val(id), path(mat), path(clusters)

    output:
      tuple val(id),
            path("${id}_umap_embedding.tsv"),
            path("${id}_tsne_embedding.tsv"),
            path("${id}_umap_clusters.png"),
            path("${id}_tsne_clusters.png"),
            path("${id}_pca_clusters.png")

    script:
    """
    python3 ${projectDir}/scripts/cluster_viz.py \
      --pca ${mat} \
      --clusters ${clusters} \
      --umap_neighbors ${params.viz_umap_neighbors} \
      --umap_min_dist ${params.viz_umap_min_dist} \
      --tsne_perplexity ${params.viz_perplexity} \
      --tsne_iter ${params.viz_tsne_iter} \
      --out_prefix ${id}
    """
}

workflow cluster_viz_ch {

  take:
    pca_out_ch
    clusters_ch

  main:
    out = CLUSTER_VIZ(clusters_ch)

  emit:
    umap_embedding = out.map{ id, umap_tsv, tsne_tsv, umap_png, tsne_png, pca_png -> umap_tsv }
    tsne_embedding = out.map{ id, umap_tsv, tsne_tsv, umap_png, tsne_png, pca_png -> tsne_tsv }
    umap_plot      = out.map{ id, umap_tsv, tsne_tsv, umap_png, tsne_png, pca_png -> umap_png }
    tsne_plot      = out.map{ id, umap_tsv, tsne_tsv, umap_png, tsne_png, pca_png -> tsne_png }
    pca_plot       = out.map{ id, umap_tsv, tsne_tsv, umap_png, tsne_png, pca_png -> pca_png }
}
