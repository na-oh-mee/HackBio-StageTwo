Single Cell RNA Sequencing Analysis

Single-cell RNA sequencing (scRNA-seq) technologies allow the dissection of gene expression at single-cell resolution, which greatly revolutionizes transcriptomic studies. Here, it is the details of single cell sequencing analysis of Peripheral Blood Mononuclear Cell (PBMC) dataset using Scanpy.
Annotated cells found in datasets are:
Neutrophil,
Plasma cells,
Platelets,
Dendritic cells,
Gamma delta T cells,
T cells, and
B cells naive.
Core python data library used for single cell analysis:
Scanpy is a powerful python library used to analyze single cell RNA seq, involves preprocessing, dimensionality reduction, clustering, and visualization.  
Anndata is a specialized data structure for storing single cell data.
igraph is used for graph construction and analysis, used internally by Scanpy to build k-nearest neighbor graphs (kNN graph).
Celltypist is a machine-learning tool for automatic cell-type annotation, uses a reference database to label clusters
Decoupler is functional analysis on gene expression data (especially single-cell).

Steps for single cell analysis include:
1. Preprocessing
2. Quality Control
3. Filtering and Normalization
4. Dimensinality reduction
5. Cell Type Annotation

# import core single cell tools
import scanpy as sc
import anndata as ad
import igraph as ig
Preprocessing
# Step 1: Loading of Dataset
bone_marrow_adata = sc.read('bone_marrow.h5ad')
Cells and Genes present in dataset
bone_marrow_adata.var.head()
bone_marrow_adata.obs.head()
# Step 2: QUALITY CONTROL
This is important to keep high quality cells and informative genes. Removes low quality cells, doublets and empty droplets
# Quality control of both cells and genes
bone_marrow_adata.var_names_make_unique()
bone_marrow_adata.obs_names_make_unique()
n_counts (total UMI / read counts per cell)

n_genes (number of genes detected per cell; usually >0)

% mitochondrial counts — high % suggests dying/stressed cells (mitochondrial gene names often MT- or mt-)

% ribosomal / hemoglobin — useful for blood or tissue contamination detection

Doublet scores — computationally inferred doublets (two cells captured as one)


# Mitochondrial genes
sc.pl.violin( bone_marrow_adata, ["pct_counts_MT"], jitter=0.4, multi_panel=False)
# number of genes in cells
sc.pl.violin(
    bone_marrow_adata,
    ["n_genes_by_counts"],
    jitter=0.4,
    multi_panel=False,
)
# Filtration of cells and genes
sc.pp.filter_cells(adata, min_genes=1000)
sc.pp.filter_genes(adata, min_cells=1000)
sc.pp.scrublet(bone_marrow_adata) #doublet detection
# Step 3: Noramlize
Normalization: is to correct for differences in sequencing depth / library size so expression values are comparable across cells.

# Normalizing to median total counts
sc.pp.normalize_total()
# Logarithmize the data
sc.pp.log1p()
# Step 4: Find Highly variable genes
Highly Variable Gene Selection: are genes that show biologically meaningful variation across cells. HVGs focuses on true biological differences, removes noise, speeds up computation, as well as improve clustering accuracy.

sc.pl.highly_variable_genes(bone_marrow_adata )
sc.pp.highly_variable_genes(bone_marrow_adata)

# Step 5: Dimensionality reduction: PCA and UMAP
Dimensionality reduction: This includes principal component analysis (PCAs), which sepearted naive states from inactive states and healthy from infected states
sc.tl.pca(bone_marrow_adata)
sc.pl.pca_variance_ratio(bone_marrow_adata, n_pcs=10, log=False)
Nonlinear Dimensionality Reduction (UMAP): creates a 2D embedding for visualization. Cluster become visible, rare population stands out
# UMAP of Hemoglobin %cell counts
sc.pl.umap(
    bone_marrow_adata,
    color=["pct_counts_HB"],
    size=8,
)
# Step 6: Cluster and Annonate
Clustering (Leiden/ Louvain): Using the neighborhood graph, algorithms like Leiden: find densely connected groups of cells, maximize modularity (network theory) and discover distinct cell populations.
Each cluster represent cell types, states or sub population
In leiden, the greater the resolution, the number the clusters
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_02", resolution=0.02)
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res2", resolution=2)
sc.pl.umap(
    bone_marrow_adata,
    color=["leiden_res0_02", "leiden_res0_5", "leiden_res2"],
)


Cell Type Annotation: Cell annotation is the process of assigning biological meaning (cell type or state) to each computationally identified cluster in your dataset.

It is the most critical step because all downstream conclusions depend on whether annotation is correct. Cell annotation combines computational methods, biological knowledge, and marker genes.

Cell annotation answers:

Are these naive or activated states?

Are these monocytes or dendritic cells?

Is this a rare stem-cell-like population?

Without annotation, the clusters are meaningless.
import pandas as pd
import decoupler as dc
# Query Omnipath and get PanglaoDB
markers = dc.op.resource(name="PanglaoDB", organism="human")

# Keep canonical cell type markers alone
markers = markers[markers["canonical_marker"]]

# Remove duplicated entries
markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

#Format because dc only accepts cell_type and genesymbol
markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
markers = markers[["source", "target"]]
markers.head()
# Cell Type Annotation
cluster_annotations = bone_marrow_adata_rank[bone_marrow_adata_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()

cluster_annotations

sc.pl.umap(bone_marrow_adata, color='cell_type', legend_loc='on data', title='Cell Type Annotation')


Other visualization include:
# Generate the dot plot
sc.pl.dotplot(
    bone_marrow_adata,
    var_names=marker_genes_dict, # Use the dictionary for grouped genes
    groupby='leiden', # Or whatever column holds your cell clusters (e.g., 'clusters', 'cell_type')
    dendrogram=True, # Set to True if you want to cluster the genes/groups
    color_map='Reds',
    use_raw=False
)


