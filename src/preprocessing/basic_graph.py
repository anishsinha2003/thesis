import scanpy as sc


def build_knn_graph(
    adata,
    n_neighbors=10,
    n_pcs=40
):
    """
    Construct a k-nearest neighbour graph using PCA embeddings.

    Parameters
    ----------
    adata : AnnData
        Preprocessed single-cell dataset.

    n_neighbors : int
        Number of neighbours connected to each cell.

    n_pcs : int
        Number of principal components used for similarity calculations.
    """

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs
    )

    return adata


def run_umap(adata):
    """
    Generate a UMAP embedding from the neighbour graph.
    """

    sc.tl.umap(adata)

    return adata


def run_leiden(
    adata,
    resolution=1.0
):
    """
    Perform Leiden clustering on the neighbour graph.
    """

    sc.tl.leiden(
        adata,
        resolution=resolution
    )

    return adata