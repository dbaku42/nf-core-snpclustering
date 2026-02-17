
process PCA_FLASHPCA {
  tag 'pca'
  cpus { params.threads as int }
  container params.flashpca_container
  publishDir "${params.outdir}/04_pca", mode: 'copy'

  input:
  tuple path(bed), path(bim), path(fam)

  output:
  tuple path('features.tsv'), path('scaled.tsv'), path('pca_scores.tsv'), path('pca_info.json')

  script:
  def flashpca_bin = params.flashpca_bin ?: 'flashpca'
  """
  set -euo pipefail

  command -v ${flashpca_bin} >/dev/null 2>&1 || { echo "ERROR: '${flashpca_bin}' not found in PATH" >&2; exit 127; }
  command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' not found" >&2; exit 127; }

  prefix="${bed.baseName}"

  ${flashpca_bin} \\
    --bfile "\$prefix" \\
    --ndim ${params.n_pcs ?: 40} \\
    --numthreads ${task.cpus} \\
    --outpc outpc.raw \\
    --outmeansd scaled.tsv

  awk '{print \$2}' "${bim}" > features.tsv

  python3 ${projectDir}/bin/flashpca_outpc_to_tsv.py \\
    --outpc outpc.raw \\
    --n-pcs ${params.n_pcs ?: 40} \\
    --out-pca pca_scores.tsv \\
    --out-info pca_info.json \\
    --id-mode fid_iid
  """
}

workflow pca_ch {
  take:
    pruned_plink_ch
  main:
    out = PCA_FLASHPCA(pruned_plink_ch)
  emit:
    out
}
