process BCFTOOLS_FILTER {

    publishDir "${params.outdir}/vcf_filtered", mode: 'copy'

    tag { vcf.baseName }

    if( params.container ) container { params.container }
    cpus params.threads
    memory '8 GB'

    input:
      path vcf

    output:
      path "*.filtered.vcf.gz"

    script:
    """
    bcftools view \
      -Oz -o ${vcf.baseName}.filtered.vcf.gz \
      -q ${params.maf}:minor \
      -e 'F_MISSING > ${params.missing}' \
      ${vcf}

    bcftools index -t ${vcf.baseName}.filtered.vcf.gz
    """
}
