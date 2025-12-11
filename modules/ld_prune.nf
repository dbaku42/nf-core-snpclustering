process LD_PRUNE {

    publishDir "${params.outdir}/vcf_pruned", mode: 'copy'

    tag { vcf.baseName }

    if( params.container ) container { params.container }
    cpus params.threads
    memory '16 GB'

    input:
      path vcf

    output:
      path "*.plink.pruned.vcf"

    script:
    """
    # 1) VCF -> PGEN with unique IDs CHR:POS and duplicate removal
    plink2 \
      --vcf ${vcf} \
      --double-id \
      --allow-extra-chr \
      --set-all-var-ids @:# \
      --rm-dup force-first \
      --make-pgen \
      --out ${vcf.baseName}.plink

    # 2) LD pruning
    plink2 \
      --pfile ${vcf.baseName}.plink \
      --allow-extra-chr \
      --indep-pairwise 50 5 0.2 \
      --out ${vcf.baseName}.plink

    # 3) Extract pruned variants and export as VCF
    plink2 \
      --pfile ${vcf.baseName}.plink \
      --allow-extra-chr \
      --extract ${vcf.baseName}.plink.prune.in \
      --export vcf \
      --out ${vcf.baseName}.plink.pruned
    """
}
