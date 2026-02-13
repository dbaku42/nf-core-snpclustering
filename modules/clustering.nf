process CLUSTERING {

    publishDir "${params.outdir}/clusters", mode: 'copy'

    tag { mat.baseName }

    if( params.container_py ) container { params.container_py }
    cpus params.threads
    memory '8 GB'

    input:
      path mat

    output:
      tuple val("${mat.baseName}"), path(mat), path("${mat.baseName}_clusters.csv")

    script:
    """
    python3 ${projectDir}/scripts/clustering.py \
      --input ${mat} \
      --n_clusters ${params.n_clusters} \
      --method ${params.cluster_method} \
      --db_eps ${params.db_eps} \
      --db_min_samples ${params.db_min_samples} \
      --out_prefix ${mat.baseName}
    """
}
workflow clustering_ch {

  take:
    pca_out_ch

  main:
    out = CLUSTERING(pca_out_ch)

  emit:
    clusters = out
}
