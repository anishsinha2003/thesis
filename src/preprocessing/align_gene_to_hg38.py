# read in adata, overlap with gencode gff file, and update a new obs column gene_name_hg38

import anndata as ad
import scanpy as sc
import pandas as pd

def preprocess_gtf_file_hg38(gff3_path):
    # ============================================================
    # Preprocess GFF3 File
    # Using gencode v36 annotation gtf file, we extract the ensembl id,
    # ============================================================

    gene_id_to_symbol = {}

    with open(gff3_path, 'r') as gff3_file:
        for line in gff3_file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'gene':
                continue
            attributes_field = fields[8]
            attributes = {}
            for attr in attributes_field.strip().split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key] = value
            # Try GFF3 conventions for gene_id and gene symbol
            gene_id = attributes.get('gene_id', attributes.get('ID', None))
            gene_symbol = attributes.get('gene_name', attributes.get('gene', None))
            if not gene_symbol:
                # Sometimes in GFF3 gene symbol is in Name or gene_name
                gene_symbol = attributes.get('Name', attributes.get('gene_name', None))
            if gene_id and gene_symbol:
                gene_id_to_symbol[gene_id] = gene_symbol

    gene_id_symbol_df = pd.DataFrame(list(gene_id_to_symbol.items()), columns=['gene_id', 'gene_symbol'])
    gene_id_symbol_df['gene_id_cleaned'] = gene_id_symbol_df['gene_id'].str.split('.', expand=True)[0]
    return gene_id_symbol_df


def preprocess_gtf_file_hg37(gtf_file):
    # ============================================================
    # Preprocess GTF File
    # Using gencode v37 annotation gtf file, we extract the ensembl id,
    # ============================================================

    # Parse gencode.v19.chr_patch_hapl_scaff.annotation.gtf to extract symbol, Ensembl gene ID, and gene type

    def parse_gtf_attributes(attr_string):
        """Parse the attribute column of a GTF file into a dict."""
        attrs = {}
        for attr in attr_string.strip().split(';'):
            if attr.strip() == '':
                continue
            key_value = attr.strip().split(' ', 1)
            if len(key_value) != 2:
                continue
            key, value = key_value
            attrs[key] = value.strip('"')
        return attrs

    gtf_gene_info = []
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type != "gene":
                continue
            attr_dict = parse_gtf_attributes(fields[8])
            # gene_id (Ensembl), gene_name (symbol), gene_type
            gene_id = attr_dict.get('gene_id')
            gene_name = attr_dict.get('gene_name')
            # GTF files may use "gene_type" or "gene_biotype"
            gene_type = attr_dict.get('gene_type', attr_dict.get('gene_biotype'))
            # Only save if gene_id and gene_name present. gene_type can be None sometimes.
            if gene_id and gene_name:
                gtf_gene_info.append({'gene_id': gene_id, 'gene_name': gene_name, 'gene_id_cleaned': gene_id.split('.')[0]})

    gencode_v19_anno_df = pd.DataFrame(gtf_gene_info)
    return gencode_v19_anno_df


def align_gene_to_hg38(adata, gff3_path, gtf_file):
    # read in gencode gff file
    hg38_gene_id_symbol_df = preprocess_gtf_file_hg38(gff3_path)
    hg37_gene_id_symbol_df = preprocess_gtf_file_hg37(gtf_file)

    # lookup dictionaries (ADDED)
    hg38_ensid_to_symbol = dict(
        zip(
            hg38_gene_id_symbol_df["gene_id_cleaned"],
            hg38_gene_id_symbol_df["gene_symbol"]
        )
    )

    hg37_symbol_to_ensid = dict(
        zip(
            hg37_gene_id_symbol_df["gene_name"],
            hg37_gene_id_symbol_df["gene_id_cleaned"]
        )
    )

    # we firstly need to determine if the adata is hg38 or hg37, also gene symbol or gene id.
    # so need to count the overlaps of adata.var_names with hg38_gene_id_symbol_df and hg37_gene_id_symbol_df
    # 1. map with gene id
    hg38_overlap_ensid = adata.var_names.isin(hg38_gene_id_symbol_df["gene_id_cleaned"])
    hg37_overlap_ensid = adata.var_names.isin(hg37_gene_id_symbol_df["gene_id_cleaned"])

    # CHANGED
    if max(hg38_overlap_ensid.sum(), hg37_overlap_ensid.sum()) > 0:

        print("The var name of the adata map with ensembl id")

        # now need to check the number of overlaps
        print(f"The number of overlaps with hg38_gene_id_symbol_df: {hg38_overlap_ensid.sum()}")
        print(f"The number of overlaps with hg37_gene_id_symbol_df: {hg37_overlap_ensid.sum()}")

        if hg38_overlap_ensid.sum() > hg37_overlap_ensid.sum():

            print("The adata is hg38")

            # CHANGED
            adata.var["gene_name_hg38"] = (
                adata.var_names.to_series(index=adata.var.index)
                .map(hg38_ensid_to_symbol)
            )

        else:

            print("The adata is hg37; mapping to hg38 by using ensembl id")

            # we need to print out the number of genes versus number of overlaps
            print(f"The number of genes: {adata.var_names.shape[0]}")
            print(f"The number of overlaps with hg38 gene ensembl ids: {hg38_overlap_ensid.sum()}")

            # calculate the percentage of overlaps, add a warning flag if it's less than 70%
            if hg38_overlap_ensid.sum() / adata.var_names.shape[0] < 0.01:
                print("Warning: The number of overlaps with hg38 gene ensembl ids is less than 1%")
                print("Please check the gencode gff file")

            # CHANGED
            adata.var["gene_name_hg38"] = (
                adata.var_names.to_series(index=adata.var.index)
                .map(hg38_ensid_to_symbol)
            )

            # CHANGED
            adata = adata[:, adata.var["gene_name_hg38"].notna()].copy()

            print(f"Dimension after filtering: {adata.shape}")

        adata.var_names = adata.var["gene_name_hg38"]
        adata.var.index.name = None

        return adata

    # 2. map with gene symbol
    else:

        print("The var name of the adata map with gene symbols")

        hg38_overlap_symbol = adata.var_names.isin(
            hg38_gene_id_symbol_df["gene_symbol"]
        )

        # CHANGED
        hg37_overlap_symbol = adata.var_names.isin(
            hg37_gene_id_symbol_df["gene_name"]
        )

        print(f"The number of overlaps with hg38_gene_id_symbol_df: {hg38_overlap_symbol.sum()}")
        print(f"The number of overlaps with hg37_gene_id_symbol_df: {hg37_overlap_symbol.sum()}")

        if hg38_overlap_symbol.sum() > hg37_overlap_symbol.sum():

            print("The adata is hg38")

            # we do nothing here
            return adata

        else:

            print("The adata is hg37")

            # we need to convert the symbol to ensembl id first, and then map to hg38

            # CHANGED
            adata.var["ensembl_id"] = (
                adata.var_names.to_series(index=adata.var.index)
                .map(hg37_symbol_to_ensid)
            )

            print(
                f"The number of ensembl ids: {adata.var['ensembl_id'].shape[0]}, "
                f"number of genes: {adata.var_names.shape[0]}"
            )

            hg38_overlap_ensid = adata.var["ensembl_id"].isin(
                hg38_gene_id_symbol_df["gene_id_cleaned"]
            )

            print(f"The number of overlaps with hg38 gene ensembl ids: {hg38_overlap_ensid.sum()}")

            # calculate the percentage of overlaps, add a warning flag if it's less than 70%
            if hg38_overlap_ensid.sum() / adata.var_names.shape[0] < 0.01:
                print("Warning: The number of overlaps with hg38 gene ensembl ids is less than 1%")
                print("Please check the gencode gff file")

            # CHANGED
            adata.var["gene_name_hg38"] = (
                adata.var["ensembl_id"]
                .map(hg38_ensid_to_symbol)
            )

            # CHANGED
            adata = adata[:, adata.var["gene_name_hg38"].notna()].copy()

            print(f"Dimension after filtering: {adata.shape}")

        adata.var_names = adata.var["gene_name_hg38"]
        adata.var.index.name = None

        return adata
