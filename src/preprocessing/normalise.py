# Make cells comparable. When sequencing happens, the machine does NOT directly measure:
# “this gene has value 10”. Instead, it reads tiny fragments of RNA millions of times
# and each fragment detected is called a "read". More read means more sequencing coverage
# So before normalisation the value at the expression value of each gene in the cells (in the matrix)
# is how many reads mapped to each gene (basically read counts).
# however, different cells often get different sequencing depth (read counts)
# so if cell A has total 50,000 reads and cell B has 5000 reads This DOES NOT necessarily mean: Cell A is biologically more active
# It might simply mean: the sequencer captured more RNA fragments from Cell A
# Without normalisation: Cell A appears more “expressed” even if biologically: Cell A and Cell B are actually similar
# What normalisation does: It rescales expression values so: all cells are on similar scale

import scanpy as sc

def normalise(adata):

    # ============================================================
    # Normalisation of Gene Expression Data
    # ============================================================

    # Different cells can contain different total sequencing counts
    # due to variations in sequencing depth and RNA capture efficiency.
    #
    # For example:
    # - one cell may contain 50,000 total reads
    # - another cell may contain 5,000 total reads
    #
    # This does not necessarily reflect true biological differences.
    #
    # Normalisation rescales each cell so that all cells have
    # approximately the same total expression level, making them
    # directly comparable for downstream analysis.

    sc.pp.normalize_total(
        adata,
        target_sum=1e4
    )

    # target_sum = 1e4 means each cell is scaled so that
    # the total counts across all genes sum to 10,000.


    # ============================================================
    # Log Transformation
    # ============================================================

    # Gene expression values are highly skewed, where a small number
    # of genes may have extremely large expression values.
    #
    # Log transformation compresses large values and stabilises
    # variance across the dataset, improving the performance of:
    #
    # - dimensionality reduction
    # - graph construction
    # - clustering
    # - machine learning models
    #
    # log1p(x) computes:
    #
    # log(x + 1)
    #
    # The +1 prevents issues with zero values.

    sc.pp.log1p(adata)

    return adata
