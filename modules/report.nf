process REPORT {

    publishDir "${params.outdir}/plots", mode: 'copy'

    tag { id }

    if( params.container_py ) container { params.container_py }
    cpus params.threads
    memory '8 GB'

    input:
      tuple val(id), path(mat), path(clusters)

    output:
      path("${id}_pca_clusters.png")
      path("${id}_tsne_clusters.png")
      path("${id}_umap_clusters.png")
      path("${id}_combined_embeddings.png")
      path("${id}_report.txt")

    script:
    """
    python3 ${projectDir}/scripts/generate_plots.py \
      --pca ${mat} \
      --clusters ${clusters} \
      --out_prefix ${id}

    printf "Analysis ID: %s\n" "${id}" > ${id}_report.txt
    printf "Input matrix: %s\n" "${mat}" >> ${id}_report.txt
    printf "Clusters: %s\n" "${clusters}" >> ${id}_report.txt
    """
}
