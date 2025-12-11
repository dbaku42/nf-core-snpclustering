process EXTRACT_MATRIX {

    publishDir "${params.outdir}/matrix", mode: 'copy'

    tag { vcf.baseName }

    if( params.container ) container { params.container }
    cpus params.threads
    memory '8 GB'

    input:
      path vcf

    output:
      path "${vcf.baseName}.genotypes.tsv"

    script:
    """
    # 1) Extract genotypes with bcftools (rows = variants, columns = samples)
    bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' ${vcf} > ${vcf.baseName}.raw.tsv

    # 2) Convert to sample x variant dosage matrix
    #    - Encode missing as -1
    #    - Drop variants with missing rate >= params.na_filter
    python3 - << 'PY'
import csv

raw = "${vcf.baseName}.raw.tsv"
out = "${vcf.baseName}.genotypes.tsv"
na_filter = float("${params.na_filter}")

with open(raw) as fh:
    reader = csv.reader(fh, delimiter='\t')
    header = next(reader)
    samples = header[2:]
    # rows_variants will hold one list of dosages per variant (length = n_samples)
    rows_variants = []
    for row in reader:
        gts = row[2:]
        dosages = []
        missing_count = 0
        for gt in gts:
            if gt in ("./.", ".", ""):
                dosages.append(-1)
                missing_count += 1
            else:
                gt_str = str(gt).replace("|", "/")
                alleles = gt_str.split("/")
                try:
                    dos = sum(1 for a in alleles if a == "1")
                except Exception:
                    dos = -1
                    missing_count += 1
                dosages.append(dos)
        if len(gts) == 0:
            continue
        missing_rate = missing_count / len(gts)
        if missing_rate < na_filter:
            rows_variants.append(dosages)

# If no variants remain, write a file with only the sample column
if not rows_variants:
    with open(out, "w", newline="") as outfh:
        w = csv.writer(outfh, delimiter="\t")
        w.writerow(["sample"])
        for s in samples:
            w.writerow([s])
else:
    # Transpose: variants -> columns (samples x variants)
    transposed = list(zip(*rows_variants))
    with open(out, "w", newline="") as outfh:
        w = csv.writer(outfh, delimiter="\t")
        w.writerow(["sample"] + [f"var_{i}" for i in range(len(rows_variants))])
        for s, vals in zip(samples, transposed):
            w.writerow([s] + list(vals))
PY
    """
}
