process CLUSTER_METRICS {

    publishDir "${params.outdir}/metrics", mode: 'copy'

    tag { id }

    if( params.container_py ) container { params.container_py }
    cpus params.threads
    memory '8 GB'

    input:
      tuple val(id), path(mat), path(clusters)

    output:
      tuple val(id),
        path("${id}_metrics.tsv"),
        path("${id}_k_scan.tsv"),
        path("${id}_elbow.png"),
        path("${id}_silhouette_vs_k.png"),
        path("${id}_calinski_vs_k.png"),
        path("${id}_dunn_vs_k.png"),
        path("${id}_metrics.png")


    script:
    """
    python3 ${projectDir}/scripts/cluster_metrics.py \
      --features ${mat} \
      --clusters ${clusters} \
      --k_min 2 \
      --k_max 8 \
      --out_prefix ${id}
    """
}
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

