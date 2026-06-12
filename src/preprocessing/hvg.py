# Highly Variable Genes - Keep the most informative genes. Single-cell datasets
# may contain: 20,000+ genes BUT: most genes are not useful many contain noise many barely vary
# Using all genes: increases computation increases noise hurts learning

import scanpy as sc

def select_hvg(
    adata,
    n_top_genes=2000
):

    # ============================================================
    # Highly Variable Gene (HVG) Selection
    # ============================================================

    # Single-cell datasets contain thousands of genes,
    # many of which contribute little useful biological information.
    #
    # Highly variable genes (HVGs) are genes whose expression
    # varies significantly across cells.
    #
    # These genes are more informative for identifying:
    # - cell states
    # - cell populations
    # - biological differences
    #
    # Selecting HVGs reduces noise and computational complexity,
    # while preserving the most important biological signals.

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000
    )

    # Keep only the highly variable genes.

    adata = adata[:, adata.var.highly_variable]

    # Print updated dataset dimensions.

    return adata
