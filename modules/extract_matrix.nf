process EXTRACTION_MATRIX {
    tag 'genotype_matrix'
    cpus { params.threads as int }
    container params.plink2_container
    publishDir "${params.outdir}/03_extraction_matrix", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path 'genotypes.raw'

    script:
    """
    set -euo pipefail

    prefix="${bed.baseName}"

    plink2 --threads ${task.cpus} \
        --bfile "$prefix" \
        --export A \
        --out genotypes
    """

}

workflow extraction_matrix_ch {
    take:
    plink_ch

    main:
    mat = EXTRACTION_MATRIX(plink_ch)

    emit:
    mat
}
