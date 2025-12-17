process LD_PRUNE {
    tag { vcf.baseName }
    cpus { params.threads as int }
    container params.plink2_container
    publishDir "${params.outdir}/02_pruning", mode: 'copy'

    input:
    path vcf

    output:
    tuple path('pruned.bed'), path('pruned.bim'), path('pruned.fam'), emit: plink
    path 'pruned.prune.in', emit: prune_in

    script:
    def window = params.ld_window_kb
    def step   = params.ld_step
    def r2     = params.ld_r2
    """
    set -euo pipefail

       plink2 --threads ${task.cpus} \
       --vcf ${vcf} \
       --allow-extra-chr \
       --max-alleles 2 \
       --snps-only just-acgt \
       --set-all-var-ids '@:#:\$r:\$a' \
       --rm-dup force-first \
       --make-bed \
       --out tmp


    # LD pruning
    plink2 --threads ${task.cpus} \
           --bfile tmp \
           --indep-pairwise ${window} ${step} ${r2} \
           --out pruned

    # Keep only pruned set
    plink2 --threads ${task.cpus} \
           --bfile tmp \
           --extract pruned.prune.in \
           --make-bed \
           --out pruned
    """
}

workflow pruning_ch {
    take:
    vcf_ch

    main:
    out = LD_PRUNE(vcf_ch)

    emit:
    out.plink
}
