#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

try:
    import umap
except ImportError:
    umap = None

def load_matrix_and_clusters(mat_path, clusters_path):
    if not os.path.exists(mat_path):
        raise SystemExit(f"ERROR: features file not found: {mat_path}")

    if not os.path.exists(clusters_path):
        raise SystemExit(f"ERROR: clusters file not found: {clusters_path}")

    # Matrix: TSV, first column = sample, rest = features
    df = pd.read_csv(mat_path, sep="\t")
    if df.shape[1] < 2:
        raise SystemExit("ERROR: matrix must have at least one feature column")
    X_raw = df.iloc[:, 1:]
    X = X_raw.apply(pd.to_numeric, errors="coerce").fillna(-1.0).to_numpy()

    # Clusters: CSV with at least column 'cluster'
    cl = pd.read_csv(clusters_path)
    if "cluster" not in cl.columns:
        raise SystemExit("ERROR: clusters CSV must have a 'cluster' column")
    labels = cl["cluster"].to_numpy()

    # Simple dimension alignment
    n = min(X.shape[0], labels.shape[0])
    if X.shape[0] != labels.shape[0]:
        print(
            f"WARNING: X has {X.shape[0]} rows, labels has {labels.shape[0]} rows; "
            f"truncating to {n}."
        )
        X = X[:n, :]
        labels = labels[:n]

    return X, labels

def scatter_embed(
    X_emb,
    labels,
    title,
    out_png,
    xlabel="dim 1",
    ylabel="dim 2",
    figsize=(5, 4),
    s=10,
    alpha=0.8,
    legend_fontsize=8,
    legend_markerscale=2,
):
    cmap = plt.get_cmap("prism")
    unique = np.unique(labels)

    plt.figure(figsize=figsize)

    # Controllo di sicurezza: se PCA ha restituito meno di 2 dimensioni (es. varianza target bassa)
    if X_emb.shape[1] < 2:
        print(f"WARNING: {title} resulted in fewer than 2 dimensions. Skipping plot.")
        plt.close()
        return

    for i, lab in enumerate(unique):
        m = labels == lab
        plt.scatter(
            X_emb[m, 0],
            X_emb[m, 1],
            s=s,
            alpha=alpha,
            color=cmap(i / max(len(unique) - 1, 1)),
            label=f"cluster {lab}",
        )

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(markerscale=legend_markerscale, fontsize=legend_fontsize)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

def main():
    ap = argparse.ArgumentParser(
        description="Generate PCA / t-SNE / UMAP plots colored by cluster"
    )

    # Input/Output
    ap.add_argument("--pca", required=True, help="TSV matrix (samples x features)")
    ap.add_argument("--clusters", required=True, help="CSV with column 'cluster'")
    ap.add_argument("--out_prefix", required=True, help="Output prefix")
    
    # Parametri Estetici
    ap.add_argument("--figsize", default="5,4", help="Figure size (e.g., '5,4')")
    ap.add_argument("--s", type=int, default=10, help="Marker size")
    ap.add_argument("--alpha", type=float, default=0.8, help="Marker alpha")
    ap.add_argument("--xlabel", default="dim 1", help="X-axis label")
    ap.add_argument("--ylabel", default="dim 2", help="Y-axis label")
    ap.add_argument("--legend_fontsize", type=int, default=8, help="Legend font size")
    ap.add_argument("--legend_markerscale", type=int, default=2, help="Legend marker scale")

    # Parametri Algoritmici (Nuovi)
    # PCA: se specifichi target_var (float), usa quello, altrimenti usa 2 componenti fisse
    ap.add_argument("--pca_target_var", type=float, default=None, 
                    help="PCA explained variance ratio (e.g., 0.95). If set, overrides n_components=2.")
    
    # t-SNE
    ap.add_argument("--tsne_perplexity", type=float, default=30.0, help="t-SNE perplexity")
    
    # UMAP
    ap.add_argument("--umap_n_neighbors", type=int, default=15, help="UMAP n_neighbors")
    ap.add_argument("--umap_min_dist", type=float, default=0.1, help="UMAP min_dist")
    ap.add_argument("--umap_metric", type=str, default="euclidean", help="UMAP metric (e.g., euclidean, manhattan, cosine)")

    args = ap.parse_args()

    # Parsing figsize
    try:
        figsize = tuple(map(int, args.figsize.split(',')))
    except ValueError:
        print("WARNING: invalid figsize format, using default (5,4)")
        figsize = (5, 4)

    X, labels = load_matrix_and_clusters(args.pca, args.clusters)

    # --- PCA ---
    # Logica: Se pca_target_var è impostato, lo passiamo a n_components (float = % varianza)
    # Altrimenti n_components = 2 (int = numero componenti)
    pca_n_comp = args.pca_target_var if args.pca_target_var is not None else 2
    
    print(f"Running PCA with n_components={pca_n_comp}...")
    pca = PCA(n_components=pca_n_comp, random_state=0)
    X_pca = pca.fit_transform(X)
    
    scatter_embed(
        X_pca,
        labels,
        "PCA (clusters)",
        f"{args.out_prefix}_pca_clusters.png",
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        figsize=figsize,
        s=args.s,
        alpha=args.alpha,
        legend_fontsize=args.legend_fontsize,
        legend_markerscale=args.legend_markerscale,
    )

    # --- t-SNE ---
    print(f"Running t-SNE with perplexity={args.tsne_perplexity}...")
    tsne = TSNE(n_components=2, perplexity=args.tsne_perplexity, random_state=0)
    X_tsne = tsne.fit_transform(X)
    scatter_embed(
        X_tsne,
        labels,
        "t-SNE (clusters)",
        f"{args.out_prefix}_tsne_clusters.png",
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        figsize=figsize,
        s=args.s,
        alpha=args.alpha,
        legend_fontsize=args.legend_fontsize,
        legend_markerscale=args.legend_markerscale,
    )

    # --- UMAP ---
    if umap is None:
        print("WARNING: umap-learn not installed. Skipping UMAP.")
        X_umap = None
    else:
        print(f"Running UMAP with n_neighbors={args.umap_n_neighbors}, min_dist={args.umap_min_dist}, metric={args.umap_metric}...")
        um = umap.UMAP(
            n_neighbors=args.umap_n_neighbors,
            min_dist=args.umap_min_dist,
            n_components=2,
            metric=args.umap_metric,
            random_state=0,
        )
        X_umap = um.fit_transform(X)
        scatter_embed(
            X_umap,
            labels,
            "UMAP (clusters)",
            f"{args.out_prefix}_umap_clusters.png",
            xlabel=args.xlabel,
            ylabel=args.ylabel,
            figsize=figsize,
            s=args.s,
            alpha=args.alpha,
            legend_fontsize=args.legend_fontsize,
            legend_markerscale=args.legend_markerscale,
        )

    # --- Combined Figure ---
    embeddings = [X_pca, X_tsne]
    names = ["PCA", "t-SNE"]
    if X_umap is not None:
        embeddings.append(X_umap)
        names.append("UMAP")
    
    # Filtra embedding validi (che hanno almeno 2 dimensioni)
    valid_embeddings = []
    valid_names = []
    for emb, nm in zip(embeddings, names):
        if emb.shape[1] >= 2:
            valid_embeddings.append(emb)
            valid_names.append(nm)
    
    if not valid_embeddings:
        print("No valid embeddings with >=2 dimensions for combined plot.")
        return

    fig, axes = plt.subplots(1, len(valid_embeddings), figsize=(4 * len(valid_embeddings), 4))
    if len(valid_embeddings) == 1:
        axes = [axes]

    cmap = plt.get_cmap("prism")
    unique = np.unique(labels)

    for ax, emb, name in zip(axes, valid_embeddings, valid_names):
        for i, lab in enumerate(unique):
            m = labels == lab
            ax.scatter(
                emb[m, 0],
                emb[m, 1],
                s=args.s,
                alpha=args.alpha,
                color=cmap(i / max(len(unique) - 1, 1)),
                label=f"{lab}" if ax is axes[0] else None,
            )
        ax.set_title(name)
        ax.set_xticks([])
        ax.set_yticks([])

    axes[0].legend(title="cluster", fontsize=args.legend_fontsize, markerscale=args.legend_markerscale)
    plt.tight_layout()
    fig.savefig(f"{args.out_prefix}_combined_embeddings.png", dpi=150)
    plt.close()

if __name__ == "__main__":
    main()
