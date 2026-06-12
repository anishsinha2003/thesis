import scanpy as sc

def run_pca(adata):

    # ============================================================
    # Principal Component Analysis (PCA)
    # ============================================================

    # Even after highly variable gene selection,
    # the dataset still contains thousands of features.
    #
    # PCA reduces the dimensionality of the dataset by
    # transforming the original gene expression space into
    # a smaller number of principal components.
    #
    # These components capture the major sources of
    # biological variation across cells.
    #
    # Dimensionality reduction improves:
    # - computational efficiency
    # - clustering
    # - graph construction
    # - downstream machine learning performance

    sc.tl.pca(adata)

    # Visualise how much variance is explained
    # by each principal component.

    # sc.pl.pca_variance_ratio(
    #     adata,
    #     log=True
    # )

    return adata