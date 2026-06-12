# quality control - Remove low-quality or noisy cells to ensure graph isnt noisy
# remove cells with too few genes, remove cells with very high mitochondrial expression
# (generally represent dying, stressed, or low-quality cells), remove genes
# expressed in very few cells


import scanpy as sc

def run_qc(adata):

    # ============================================================
    # Calculate Quality Control (QC) Metrics
    # ============================================================

    # QC metrics help evaluate the quality of each cell.
    #
    # Examples include:
    # - number of genes detected per cell
    # - total sequencing counts per cell
    # - expression distribution statistics
    #
    # These metrics are stored in:
    # adata.obs

    sc.pp.calculate_qc_metrics(
        adata,
        inplace=True
    )

    # View the first few rows of cell metadata
    # containing QC statistics.

    adata.obs.head()

    # ============================================================
    # Visualise Quality Control Distributions
    # ============================================================

    # Violin plots are used to visualise:
    #
    # - n_genes_by_counts:
    #   number of genes detected in each cell
    #
    # - total_counts:
    #   total sequencing reads/counts per cell
    #
    # This helps identify:
    # - low-quality cells
    # - sequencing outliers
    # - abnormal cells

    # sc.pl.violin(
    #     adata,
    #     ['n_genes_by_counts', 'total_counts'],
    #     jitter=0.4
    # )


    # ============================================================
    # Filter Low-Quality Cells and Genes
    # ============================================================

    # Remove cells expressing fewer than 200 genes.
    # These are often:
    # - dead cells
    # - empty droplets
    # - poor-quality captures

    sc.pp.filter_cells(
        adata,
        min_genes=200
    )

    # Remove genes expressed in fewer than 3 cells.
    # These genes are usually uninformative or noisy.

    sc.pp.filter_genes(
        adata,
        min_cells=3
    )

    # Print updated dataset dimensions after filtering.
    return adata