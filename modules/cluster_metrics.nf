process CLUSTER_METRICS {

    publishDir "${params.outdir}/metrics", mode: 'copy'

    tag { id }

    if( params.container_py ) container { params.container_py }
    cpus   params.threads
    memory '8 GB'

    input:
        // FIX: tupla a 5 elementi — aggiunto pca_scores (joinato in main.nf)
        // mat (features.tsv) è la lista SNP di FlashPCA, non ha sample_id
        // pca_scores.tsv ha sample_id + colonne PC -> quello che serve per le metriche
        tuple val(id),
              path(mat),
              path(clusters),
              path(clustering_info),
              path(pca_scores)

    output:
<<<<<<< Updated upstream
      tuple val(id),
        path("${id}_metrics.tsv"),
        path("${id}_k_scan.tsv"),
        path("${id}_elbow.png"),
        path("${id}_silhouette_vs_k.png"),
        path("${id}_calinski_vs_k.png"),
        path("${id}_dunn_vs_k.png"),
        path("${id}_metrics.png")

=======
        tuple val(id),
              path("${id}_k_sweep.csv"),
              path("${id}_selected.json"),
              path("${id}_elbow.png"),
              path("${id}_silhouette.png"),
              path("${id}_davies_bouldin.png"),
              path("${id}_calinski.png")
>>>>>>> Stashed changes

    script:
    // FIX: --features usa pca_scores (ha sample_id), non mat (lista SNP)
    """
    python3 ${projectDir}/scripts/cluster_metrics.py \
      --features    ${pca_scores} \
      --clusters    ${clusters} \
      --k-min       ${params.k_min} \
      --k-max       ${params.k_max} \
      --out-k-sweep ${id}_k_sweep.csv \
      --out-selected ${id}_selected.json \
      --out-prefix  ${id}
    """
}
<<<<<<< Updated upstream
workflow cluster_metrics_ch {

  take:
    pca_out_ch
    clusters_ch

  main:
    out = CLUSTER_METRICS(clusters_ch)

  emit:
    selected       = out.map{ id, metrics_tsv, k_scan_tsv, elbow_png, sil_png, cal_png, dunn_png, metrics_png -> metrics_tsv }
    k_sweep        = out.map{ id, metrics_tsv, k_scan_tsv, elbow_png, sil_png, cal_png, dunn_png, metrics_png -> k_scan_tsv }
    elbow          = out.map{ id, metrics_tsv, k_scan_tsv, elbow_png, sil_png, cal_png, dunn_png, metrics_png -> elbow_png }
    silhouette     = out.map{ id, metrics_tsv, k_scan_tsv, elbow_png, sil_png, cal_png, dunn_png, metrics_png -> sil_png }
    davies_bouldin = out.map{ id, metrics_tsv, k_scan_tsv, elbow_png, sil_png, cal_png, dunn_png, metrics_png -> metrics_png }
}

=======

workflow cluster_metrics_ch {

    take:
        clusters_ch   // tuple: [id, mat, clusters, clustering_info, pca_scores]

    main:
        out = CLUSTER_METRICS(clusters_ch)

    emit:
        k_sweep        = out.map{ id, ksweep, sel, elbow, sil, db, cal -> ksweep }
        selected       = out.map{ id, ksweep, sel, elbow, sil, db, cal -> sel }
        elbow          = out.map{ id, ksweep, sel, elbow, sil, db, cal -> elbow }
        silhouette     = out.map{ id, ksweep, sel, elbow, sil, db, cal -> sil }
        davies_bouldin = out.map{ id, ksweep, sel, elbow, sil, db, cal -> db }
        calinski       = out.map{ id, ksweep, sel, elbow, sil, db, cal -> cal }
}
>>>>>>> Stashed changes
