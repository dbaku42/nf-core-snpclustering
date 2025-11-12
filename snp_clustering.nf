#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ===============================
// Parameters (overridable via CLI)
// ===============================
params.vcf_file             = "input.vcf.gz"     // --vcf_file <path>
params.out_dir              = "results"          // --out_dir <dir>

params.variant_type         = "SNPs"             // "SNPs" | "INDELs" | "ALL"
params.maf_filter           = 0.01               // --maf_filter 0.01
params.na_filter            = 0.10               // --na_filter 0.10

// LD pruning (PLINK2)
params.ld_prune             = true               // --ld_prune true|false
params.plink_window_kb      = 250                // --plink_window_kb 250
params.plink_polyploid_mode = 'missing'          // --plink_polyploid_mode 'missing'|'error'|'allow' (default 'missing')
// plink step è 1 in modalità 'kb' con PLINK2 (indicato in script)
params.plink_r2             = 0.2                // --plink_r2 0.2

// PCA + t-SNE + clustering
params.pca_target_variance  = 0.80               // --pca_target_variance 0.8
params.tsne_perplexity      = 30                 // --tsne_perplexity 50
params.n_clusters           = 4                  // --n_clusters 3
params.max_k_elbow          = 10                 // --max_k_elbow 10

// Resources
params.mem_pca              = '16 GB'            // --mem_pca '24 GB'
params.mem_cluster          = '8 GB'             // --mem_cluster '12 GB'

// ===============================
// Workflow
// ===============================
workflow {
    // 0) Input
    vcf_ch = Channel.fromPath(params.vcf_file)

    // 1) Prefiltro VCF (MAF + tipo variante)
	vcf_id_ch = vcf_ch.map { file -> tuple(file.baseName.replaceAll(/\.(vcf|vcf\.gz)$/, ''), file) }
	prefiltered_ch = prefilter_vcf(vcf_id_ch)

    // 2) LD pruning opzionale con PLINK2 -> VCF pruned
    if (params.ld_prune) {
        pgen_ch        = vcf_to_plink2(prefiltered_ch)               // (id, .pgen, .pvar, .psam)
        pruned_pgen_ch = plink2_ld_prune(pgen_ch)                     // (id, .pruned.pgen/.pvar/.psam, .prune.in/.out)
        pruned_vcf_ch  = plink2_export_pruned_vcf(
                           pruned_pgen_ch.map { id, pgen, pvar, psam, pin, pout -> tuple(id, pgen, pvar, psam) }
                         )                                            // (id, .pruned.vcf.gz)
        vcf_for_matrix_ch = pruned_vcf_ch                             // (id, vcf.gz)
    } else {
        vcf_for_matrix_ch = prefiltered_ch                            // (id, vcf.gz)
    }

    // 3) Estrazione matrice genotipi (samples x variants), imputazione mean
    matrix_ch = extract_genotype_matrix(vcf_for_matrix_ch)            // (id, genotype_matrix_imputed.npy, samples.txt)

    // 4) PCA su varianza cumulativa target (streaming)
    pca_ch = perform_pca(matrix_ch)                                   // (id, pca_coords.csv, pca_plot.png)

    // 5) t‑SNE sui punteggi PCA
    tsne_ch = perform_tsne_visualization(pca_ch)                      // (id, tsne_coords.csv, tsne_plot.png)

    // 6) Join PCA + t‑SNE per clustering e plot colorati
    pca_kv  = pca_ch .map { id, pca_csv,  pca_png  -> tuple(id, [pca_csv,  pca_png]) }
    tsne_kv = tsne_ch.map { id, tsne_csv, tsne_png -> tuple(id, [tsne_csv, tsne_png]) }
    clust_in = pca_kv.join(tsne_kv)                                   // join su id
                    .map { id, p, t -> tuple(id, p[0], p[1], t[0], t[1]) } // (id, pca_csv, pca_png, tsne_csv, tsne_png)
    perform_clustering(clust_in)                                      // cluster_labels.csv + plot
}

// ===============================
// Process 1: Prefilter VCF
// ===============================
process prefilter_vcf {
    publishDir "${params.out_dir}/prefiltered_vcf", mode: 'copy'
    conda 'bioconda::bcftools=1.15.1'
    tag "${vcf_file.baseName}"

    input:
    tuple val(id), path(vcf_file)

    output:
    tuple val(id), path("filtered.vcf.gz")

    when:
    true

    script:
    def vtype = params.variant_type.toUpperCase()
    def type_filter = vtype == 'ALL' ? '' : "--types ${params.variant_type.toLowerCase()}"
    // prefisso 'id' senza estensioni per PLINK2
    def id = vcf_file.getName().replaceAll(/\.vcf\.gz$/, '').replaceAll(/\.vcf$/, '')
    """
    set -euo pipefail
    bcftools view -i 'MAF > ${params.maf_filter}' ${type_filter} -O z -o filtered.vcf.gz ${vcf_file}
    """
}

// ===============================
// Process 2a: VCF -> PLINK2 pgen
// ===============================
process vcf_to_plink2 {
    publishDir "${params.out_dir}/plink", mode: 'copy'
    conda 'bioconda::plink2=2.00a3'
    tag "${name}"

    input:
    tuple val(name), path(vcf_gz)

    output:
    tuple val(name), path("${name}.pgen"), path("${name}.pvar"), path("${name}.psam")

    script:
    """
    set -euo pipefail
    # assicurati che il VCF sia indicizzato con tabix (utile per chunking/controlli)
    if [[ ! -f "${vcf_gz}.tbi" ]]; then
      if command -v tabix >/dev/null 2>&1; then
        tabix -p vcf ${vcf_gz} || true
      fi
    fi
    # Prefisso --out senza estensioni
    plink2 --vcf ${vcf_gz} --polyploid-mode ${params.plink_polyploid_mode} --make-pgen --out ${name}
    """
}

// ===============================
// Process 2b: PLINK2 LD prune (kb, step=1)
// ===============================
process plink2_ld_prune {
    publishDir "${params.out_dir}/plink_pruned", mode: 'copy'
    conda 'bioconda::plink2=2.00a3'
    tag "${name}"

    input:
    tuple val(name), path(pgen), path(pvar), path(psam)

    output:
    tuple val(name),
          path("${name}.pruned.pgen"),
          path("${name}.pruned.pvar"),
          path("${name}.pruned.psam"),
          path("${name}.prune.in"),
          path("${name}.prune.out")

    script:
    """
    set -euo pipefail
    # LD pruning in kb, step fisso a 1
    plink2 --pfile ${name} --indep-pairwise ${params.plink_window_kb}kb 1 ${params.plink_r2} --out ${name}
    # Applica il pruning
    plink2 --pfile ${name} --extract ${name}.prune.in --make-pgen --out ${name}.pruned
    """
}

// ===============================
// Process 2c: Esporta VCF pruned
// ===============================
process plink2_export_pruned_vcf {
    publishDir "${params.out_dir}/plink_pruned", mode: 'copy'
    conda 'bioconda::plink2=2.00a3'
    tag "${name}"

    input:
    tuple val(name), path(pruned_pgen), path(pruned_pvar), path(pruned_psam)

    output:
    tuple val(name), path("${name}.pruned.vcf.gz")

    script:
    """
    set -euo pipefail
    plink2 --pfile ${name}.pruned --recode vcf bgz --out ${name}.pruned
    """
}

// ===============================
// Process 3: Estrazione matrice genotipi (samples x variants)
// ===============================
process extract_genotype_matrix {
    publishDir "${params.out_dir}/genotype_data", mode: 'copy'
    conda 'bioconda::scikit-allel=1.3.2 conda-forge::pandas=1.3.4 conda-forge::scikit-learn=1.0.2'
    tag "${name}"

    input:
    tuple val(name), path(vcf_gz)

    output:
    tuple val(name), path("genotype_matrix_imputed.npy"), path("samples.txt")

    script:
    """
    #!/usr/bin/env python
    import allel, numpy as np, pandas as pd
    from sklearn.impute import SimpleImputer

    v = allel.read_vcf('${vcf_gz}', fields=['samples','calldata/GT'])
    if v is None or len(v['samples']) == 0:
        np.save('genotype_matrix_imputed.npy', np.array([])); open('samples.txt','w').close(); raise SystemExit(0)

    gt = allel.GenotypeArray(v['calldata/GT'])     # (variants, samples, ploidy)
    gm = gt.to_n_alt()                              # (variants, samples), missing -> -1

    missing_rate = (gm == -1).sum(axis=1) / gm.shape[1]
    keep = missing_rate < ${params.na_filter}
    gm_f = gm[keep, :]

    if gm_f.shape[0] == 0:
        np.save('genotype_matrix_imputed.npy', np.array([])); open('samples.txt','w').close(); raise SystemExit(0)

    # Imputazione mean per variante e formato finale (samples x variants)
    X = gm_f.T
    X_imp = SimpleImputer(missing_values=-1, strategy='mean').fit_transform(X)

    np.save('genotype_matrix_imputed.npy', X_imp)
    pd.Series(v['samples']).to_csv('samples.txt', index=False, header=False)
    """
}

// ===============================
// Process 4: PCA (IncrementalPCA, varianza target)
// ===============================
process perform_pca {
    publishDir "${params.out_dir}/pca_results", mode: 'copy'
    conda 'conda-forge::scikit-learn=1.0.2 conda-forge::pandas=1.3.4 conda-forge::matplotlib-base=3.5.0 conda-forge::seaborn-base=0.11.2'
    cpus 2
    memory { params.mem_pca }
    tag "${name}"

    input:
    tuple val(name), path(imputed_matrix), path(samples_list)

    output:
    tuple val(name), path("pca_coords.csv"), path("pca_plot.png")

    script:
    """
    #!/usr/bin/env python
    import numpy as np, pandas as pd
    from sklearn.decomposition import IncrementalPCA
    import matplotlib.pyplot as plt, seaborn as sns

    X_mm = np.load('${imputed_matrix}', mmap_mode='r')   # (samples, variants)
    n_samples, n_features = X_mm.shape if X_mm.size else (0, 0)
    if X_mm.size == 0 or n_samples < 2 or n_features < 2:
        pd.DataFrame(columns=['Sample']).to_csv('pca_coords.csv', index=False)
        plt.figure(); plt.text(0.5,0.5,'Not enough data for PCA', ha='center', va='center'); plt.savefig('pca_plot.png'); raise SystemExit(0)

    mean = np.asarray(X_mm.mean(axis=0), dtype=np.float32)
    max_comp = int(min(n_samples-1, n_features, 200))
    batch = 1024

    ipca = IncrementalPCA(n_components=max_comp, batch_size=batch)
    for s in range(0, n_samples, batch):
        e = min(s+batch, n_samples)
        chunk = X_mm[s:e, :].astype(np.float32, copy=False); chunk -= mean
        ipca.partial_fit(chunk)

    Z_parts = []
    for s in range(0, n_samples, batch):
        e = min(s+batch, n_samples)
        chunk = X_mm[s:e, :].astype(np.float32, copy=False); chunk -= mean
        Z_parts.append(ipca.transform(chunk))
    Z = np.vstack(Z_parts)

    ev = ipca.explained_variance_ratio_; cum = ev.cumsum()
    target = float(${params.pca_target_variance})
    k = int(np.searchsorted(cum, target) + 1); k = max(1, min(k, Z.shape[1]))

    Zk = Z[:, :k]
    samples = pd.read_csv('${samples_list}', header=None)[0].tolist()
    cols = [f'PC{i+1}' for i in range(k)]
    df = pd.DataFrame(Zk, columns=cols); df['Sample'] = samples
    df = df[['Sample'] + cols]; df.to_csv('pca_coords.csv', index=False)

    plt.figure(figsize=(12,8))
    if k >= 2:
        ax = sns.scatterplot(x='PC1', y='PC2', data=df, s=100)
        plt.xlabel(f'PC1 ({ev[0]:.2%})'); plt.ylabel(f'PC2 ({ev[1]:.2%})')
        for i, r in df.iterrows(): ax.text(r['PC1'], r['PC2'], r['Sample'], fontsize=8)
    else:
        plt.scatter(df['PC1'], [0]*len(df), s=100)
        for i, r in df.iterrows(): plt.text(r['PC1'], 0.02, r['Sample'], fontsize=8)
        plt.xlabel(f'PC1 ({ev[0]:.2%})'); plt.ylabel('PC2 (NA)')
    plt.title(f'PCA (components kept: {k}, total var: {cum[k-1]:.2%})'); plt.grid(True); plt.savefig('pca_plot.png')
    """
}

// ===============================
// Process 5: t‑SNE  on PCA_data
// ===============================
process perform_tsne_visualization {
    publishDir "${params.out_dir}/tsne_results", mode: 'copy'
    conda 'conda-forge::scikit-learn=1.0.2 conda-forge::pandas=1.3.4 conda-forge::matplotlib-base=3.5.0 conda-forge::seaborn-base=0.11.2'
    tag "${name}"

    input:
    tuple val(name), path(pca_coords_csv), path(pca_plot)

    output:
    tuple val(name), path("tsne_coords.csv"), path("tsne_plot.png")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib.pyplot as plt, seaborn as sns
    from sklearn.manifold import TSNE

    df = pd.read_csv('${pca_coords_csv}')
    pcs = [c for c in df.columns if c.startswith('PC')]
    if len(df) <= 1 or not pcs:
        pd.DataFrame(columns=['Sample','TSNE1','TSNE2']).to_csv('tsne_coords.csv', index=False)
        plt.figure(); plt.text(0.5,0.5,'Not enough data for t-SNE', ha='center', va='center'); plt.savefig('tsne_plot.png'); raise SystemExit(0)

    X = df[pcs]
    perpl = min(float(${params.tsne_perplexity}), len(df)-1.0)
    try:
        tsne = TSNE(n_components=2, perplexity=perpl, random_state=42, max_iter=1000)
    except TypeError:
        tsne = TSNE(n_components=2, perplexity=perpl, random_state=42, n_iter=1000)

    T = tsne.fit_transform(X)
    out = pd.DataFrame(T, columns=['TSNE1','TSNE2']); out['Sample'] = df['Sample']
    out[['Sample','TSNE1','TSNE2']].to_csv('tsne_coords.csv', index=False)

    plt.figure(figsize=(12,8))
    ax = sns.scatterplot(x='TSNE1', y='TSNE2', data=out, s=100)
    plt.title(f't-SNE (using {len(pcs)} PCs, perplexity: {perpl})')
    for i, r in out.iterrows(): ax.text(r['TSNE1'], r['TSNE2'], r['Sample'], fontsize=8)
    plt.grid(True); plt.savefig('tsne_plot.png')
    """
}

// ===============================
// Process 6: Clustering on PCA data + colored plots and metrics
// ===============================
process perform_clustering {
    publishDir "${params.out_dir}/clustering_results", mode: 'copy'
    conda 'conda-forge::scikit-learn=1.0.2 conda-forge::pandas=1.3.4 conda-forge::matplotlib-base=3.5.0'
    cpus 2
    memory { params.mem_cluster }
    tag "${name}"

    input:
    tuple val(name), path(pca_coords_csv), path(pca_plot), path(tsne_coords_csv), path(tsne_plot)

    output:
    tuple val(name),
          path("cluster_labels.csv"),
          path("elbow_plot.png"),
          path("silhouette_plot.png"),
          path("pca_cluster_plot.png"),
          path("tsne_cluster_plot.png")

    script:
    """
    #!/usr/bin/env python
    import numpy as np, pandas as pd, matplotlib.pyplot as plt
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.metrics import silhouette_score

    pca_df  = pd.read_csv('${pca_coords_csv}')
    tsne_df = pd.read_csv('${tsne_coords_csv}')
    pc_cols = [c for c in pca_df.columns if c.startswith('PC')]

    if len(pca_df) < 2 or not pc_cols:
        pd.DataFrame(columns=['Sample','Cluster']).to_csv('cluster_labels.csv', index=False)
        for fn in ['elbow_plot.png','silhouette_plot.png','pca_cluster_plot.png','tsne_cluster_plot.png']:
            plt.figure(); plt.text(0.5,0.5,'Not enough data', ha='center', va='center'); plt.savefig(fn); plt.clf()
        raise SystemExit(0)

    X = pca_df[pc_cols].to_numpy(dtype=np.float32, copy=False)
    samples = pca_df['Sample'].tolist()

    ks = list(range(2, min(${params.max_k_elbow}, len(pca_df)-1) + 1))
    inertia, sils = [], []
    for k in ks:
        km = MiniBatchKMeans(n_clusters=k, random_state=42, batch_size=1024, n_init='auto', max_iter=100)
        labels = km.fit_predict(X)
        inertia.append(km.inertia_)
        try:
            sils.append(silhouette_score(X, labels) if 2 <= len(set(labels)) <= (len(X)-1) else np.nan)
        except Exception:
            sils.append(np.nan)

    plt.figure(figsize=(10,6)); plt.plot(ks, inertia, marker='o')
    plt.title('Elbow Method (PCA space)'); plt.xlabel('k'); plt.ylabel('Inertia'); plt.grid(True)
    plt.savefig('elbow_plot.png'); plt.clf()

    plt.figure(figsize=(10,6)); plt.plot(ks, sils, marker='o', color='orange')
    plt.title('Silhouette score vs k (PCA space)'); plt.xlabel('k'); plt.ylabel('Silhouette score'); plt.grid(True)
    plt.savefig('silhouette_plot.png'); plt.clf()

    final_k = min(${params.n_clusters}, max(2, len(pca_df)-1))
    km = MiniBatchKMeans(n_clusters=final_k, random_state=42, batch_size=1024, n_init='auto', max_iter=200)
    labels = km.fit_predict(X)
    pd.DataFrame({'Sample': samples, 'Cluster': labels}).to_csv('cluster_labels.csv', index=False)

    cmap = plt.cm.get_cmap('tab10', final_k)

    # PCA colorato per cluster
    plt.figure(figsize=(12,8))
    if len(pc_cols) >= 2:
        x, y = pca_df['PC1'].to_numpy(), pca_df['PC2'].to_numpy()
        plt.scatter(x, y, c=labels, cmap=cmap, s=100, edgecolor='k')
        for xi, yi, name in zip(x, y, samples): plt.text(xi, yi, name, fontsize=8, ha='left', va='bottom')
        plt.xlabel('PC1'); plt.ylabel('PC2'); plt.title(f'PCA colored by clusters (k={final_k})'); plt.grid(True)
    else:
        x = pca_df['PC1'].to_numpy()
        plt.scatter(x, np.zeros_like(x), c=labels, cmap=cmap, s=100, edgecolor='k')
        for xi, name in zip(x, samples): plt.text(xi, 0.02, name, fontsize=8, ha='left', va='bottom')
        plt.xlabel('PC1'); plt.ylabel('PC2'); plt.title(f'PCA colored by clusters (k={final_k})'); plt.grid(True)
    plt.savefig('pca_cluster_plot.png'); plt.clf()

    # t‑SNE colorato (allineamento su Sample)
    if set(['Sample','TSNE1','TSNE2']).issubset(tsne_df.columns):
        m = tsne_df.merge(pd.DataFrame({'Sample': samples, 'Cluster': labels}), on='Sample', how='inner')
        plt.figure(figsize=(12,8))
        plt.scatter(m['TSNE1'], m['TSNE2'], c=m['Cluster'], cmap=cmap, s=100, edgecolor='k')
        for xi, yi, name in zip(m['TSNE1'], m['TSNE2'], m['Sample']): plt.text(xi, yi, name, fontsize=8, ha='left', va='bottom')
        plt.title(f't-SNE colored by clusters (k={final_k})'); plt.xlabel('TSNE1'); plt.ylabel('TSNE2'); plt.grid(True)
        plt.savefig('tsne_cluster_plot.png'); plt.clf()
    else:
        plt.figure(); plt.text(0.5,0.5,'No t-SNE data', ha='center', va='center'); plt.savefig('tsne_cluster_plot.png'); plt.clf()
    """
}
