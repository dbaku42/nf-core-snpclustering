process CLUSTERING {

    publishDir "${params.outdir}/clusters", mode: 'copy'

    // FIX 1: tag usa features.baseName, non mat.baseName (mat era la tupla intera)
    tag { features.baseName }

    if( params.container_py ) container { params.container_py }
    cpus   params.threads
    memory '8 GB'

    input:
        // FIX 2: destruttura la tupla 4-file emessa da pca_ch invece di path mat singolo
        tuple path(features), path(scaled), path(pca_scores), path(pca_info)

    output:
        // FIX 3: tupla a 4 elementi — aggiunto _clustering_info.json
        // FIX 4: prefix = features.baseName (stringa pulita, no parentesi)
        tuple val(features.baseName),
              path(features),
              path("${features.baseName}_clusters.csv"),
              path("${features.baseName}_clustering_info.json")

    script:
    // FIX 5: --method usa params.algorithm (definito in main.nf), non params.cluster_method (undefined)
    // FIX 6: aggiunti --n_init e --init-method (trattino, non underscore)
    // FIX 7: file passati singolarmente via --input, prefix pulito via --out_prefix
    """
    python3 ${projectDir}/scripts/clustering.py \
      --input ${features} ${scaled} ${pca_scores} ${pca_info} \
      --n_clusters ${params.n_clusters} \
      --method ${params.algorithm} \
      --n_init ${params.n_init ?: 100} \
      --init-method ${params.init_method ?: 'random'} \
      --db_eps ${params.dbscan_eps} \
      --db_min_samples ${params.dbscan_min_samples} \
      --out_prefix ${features.baseName}
    """
}
<<<<<<< Updated upstream
workflow clustering_ch {

  take:
    pca_out_ch

  main:
    out = CLUSTERING(pca_out_ch)

  emit:
    clusters = out
=======

workflow clustering_ch {

    take:
        pca_out_ch

    main:
        out = CLUSTERING(pca_out_ch)

    emit:
        clusters = out
>>>>>>> Stashed changes
}
