import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import infercnvpy as cnv
import numpy as np
import warnings
from anndata import AnnData
from matplotlib.pyplot import rc_context,plot,savefig
import matplotlib.pyplot as pl
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
print(os.getcwd())
# load adata
os.chdir('/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage')
DATA_PATH = "/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage"
H5AD_FILE = os.path.join(DATA_PATH, "cancer_scrna_stage2.h5ad")
adata = sc.read_h5ad(H5AD_FILE)
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

cnv.io.genomic_position_from_gtf(gtf_file="/share/pub/dengcy/refGenome/gencode.v34.annotation.gtf", adata=adata)
adata.var.loc[:, ["chromosome", "start", "end"]].head()
cnv.tl.infercnv(
    adata,
    reference_key="cnv_annotation",
    reference_cat=["Macrophage","Fibroblasts"],
    window_size=250,
)
cnv.pl.chromosome_heatmap(adata, groupby="cnv_annotation",save="stage2_infercnv_chromsomeHeatmap.pdf")
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True,save="stage2_infercnv_chromsomeHeatmap2.pdf")
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)
#fig, ((ax1, ax2), (ax3, ax4)) = pl.subplots(2, 2, figsize=(11, 11))
#ax4.axis("off")
cnv.pl.umap(
    adata,
    color="cnv_leiden",
    show=False,
    save="stage2_cnv_leiden_umap.pdf"
)
cnv.pl.umap(adata, color="cnv_score", show=False,save="stage2_cnv_score_umap.pdf")
cnv.pl.umap(adata, color="cnv_annotation", show=False,save="stage2_cnv_annotation_umap.pdf")
adata.obs["cnv_status"] = "tumor"
adata.obs.loc[
    adata.obs["cnv_leiden"].isin(["3", "5", "7", "1"]), "cnv_status"
] = "normal"
#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw=dict(wspace=0.5))
cnv.pl.umap(adata, color="cnv_status",show=False,save="stage2_cnv_NT_umap.pdf")
#sc.pl.umap(adata, color="cnv_status",show=False,save="stage2_NT_umap.pdf")
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "tumor", :],save="stage2_infercnv_chromsomeHeatmap_tumor.pdf")
cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "normal", :],save="stage2_infercnv_chromsomeHeatmap_normal.pdf")
adata.obs.to_csv('/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_infercnv_stage2.csv')
