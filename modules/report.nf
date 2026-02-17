process REPORT {

    publishDir "${params.outdir}/report", mode: 'copy'

    tag { id }

    if( params.container_py ) container { params.container_py }
    cpus   params.threads
    memory '8 GB'

    input:
        // Tutti i file in un'unica tupla — il join è fatto in main.nf
        tuple val(id),
              path(mat),
              path(clusters),
              path(clustering_info),
              path(pca_info_json),
              path(k_sweep_csv),
              path(selected_json),
              path(elbow_png),
              path(silhouette_png),
              path(davies_bouldin_png),
              path(calinski_png),
              path(umap_embedding),
              path(tsne_embedding),
              path(umap_png),
              path(tsne_png),
              path(pca_cluster_png)

    output:
        tuple val(id), path("${id}_report.html")

    script:
    """
    python3 ${projectDir}/scripts/generate_report.py \
      --pca-info           ${pca_info_json} \
      --clusters           ${clusters} \
      --k-sweep            ${k_sweep_csv} \
      --selected-metrics   ${selected_json} \
      --umap-embedding     ${umap_embedding} \
      --tsne-embedding     ${tsne_embedding} \
      --umap-png           ${umap_png} \
      --tsne-png           ${tsne_png} \
      --elbow-png          ${elbow_png} \
      --silhouette-png     ${silhouette_png} \
      --davies-bouldin-png ${davies_bouldin_png} \
      --pca-cluster-png    ${pca_cluster_png} \
      --out                ${id}_report.html
    """
}
<<<<<<< Updated upstream
workflow report_ch {

  take:
    clusters_ch

  main:
    out = REPORT(clusters_ch)

  emit:
    pca_plot  = out.map{ id, pca_png, tsne_png, umap_png, comb_png, txt -> pca_png }
    tsne_plot = out.map{ id, pca_png, tsne_png, umap_png, comb_png, txt -> tsne_png }
    umap_plot = out.map{ id, pca_png, tsne_png, umap_png, comb_png, txt -> umap_png }
    combined  = out.map{ id, pca_png, tsne_png, umap_png, comb_png, txt -> comb_png }
    report    = out.map{ id, pca_png, tsne_png, umap_png, comb_png, txt -> txt }
=======

workflow report_ch {

    take:
        combined_ch   // unico canale con tupla completa — join fatto in main.nf

    main:
        out = REPORT(combined_ch)

    emit:
        report = out.map{ id, html -> html }
>>>>>>> Stashed changes
}
