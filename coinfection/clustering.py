import scanpy as sc
import pandas as pd

# adata: AnnData with raw counts or normalized expression; adata.obs has metadata
co_mask = adata.obs['infection_status'] == 'co-infection'
adata_co = adata[co_mask].copy()
