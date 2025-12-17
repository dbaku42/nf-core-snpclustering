process PREPROCESS_VCF {
    tag { vcf.baseName }
    cpus { params.threads as int }
    container params.bcftools_container
    publishDir "${params.outdir}/01_preprocess", mode: 'copy'

    input:
    path vcf

    output:
    path 'filtered.vcf.gz',     emit: vcf
    path 'filtered.vcf.gz.tbi', emit: tbi

    script:
    def maf = params.maf
    def miss = params.missing
    """
    set -euo pipefail

    # Add INFO tags if missing and filter.
    bcftools +fill-tags ${vcf} -Oz -o filled.vcf.gz -- -t MAF,F_MISSING
    bcftools view -i 'MAF>=${maf} && F_MISSING<=${miss}' -Oz -o filtered.vcf.gz filled.vcf.gz
    tabix -f -p vcf filtered.vcf.gz
    """
}

workflow preprocess_ch {
    take:
    vcf_ch

    main:
    out = PREPROCESS_VCF(vcf_ch)

    emit:
    out.vcf
}
