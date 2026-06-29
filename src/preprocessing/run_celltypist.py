import scanpy as sc
import celltypist.models
import celltypist.annotate

def run_celltypist(adata, model_name="Immune_All_Low.pkl"):

    # ============================================================
    # CellTypist annotation
    # ============================================================

    # CellTypist is a machine learning model that annotates cells based on their gene expression profile.
    # It uses a pre-trained model to predict the cell type of each cell.
    # The model is trained on a large dataset of single-cell RNA-seq data.

    model = celltypist.models.Model.load(model_name)
    predictions = celltypist.annotate(adata, model=model)

    # adata.obs['celltypist_predictions'] = predictions.predictions
    adata = predictions.to_adata()

    return adata