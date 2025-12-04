\# HackBio-StageTwo  
Single Cell RNA Sequencing Analysis

Single-cell RNA sequencing (scRNA-seq) technologies allow the dissection of gene expression at single-cell resolution, which greatly revolutionizes transcriptomic studies. Here, it is the details of single cell sequencing analysis using Scanpy.  
Core python data library used for single cell analysis:  
Scanpy is a powerful python library used to analyze single cell RNA seq, involves preprocessing, dimensionality reduction, clustering, and visualization.    
Anndata is a specialized data structure for storing single cell data.  
igraph is used for graph construction and analysis, used internally by Scanpy to build k-nearest neighbor graphs (kNN graph).  
Celltypist is a machine-learning tool for automatic cell-type annotation, uses a reference database to label clusters  
Decoupler is functional analysis on gene expression data (especially single-cell).

Steps for single cell analysis include:  
1\. Preprocessing  
2\. Quality Control  
3\. Filtering and Normalization  
4\. Dimensinality reduction  
5\. Cell Type Annotation

\#\#\# import core single cell tools  
'''  
     import scanpy as sc  
                         '''

'''  
     import anndata as ad  
                        '''

'''  
     import igraph as ig  
                       '''  
                       

Preprocessing  
\#\#\# Step 1: Loading of Dataset  
'''  
        bone\_marrow\_adata \= sc.read('bone\_marrow.h5ad')  
                                              '''

\#\#\#Cells and Genes present in dataset

'''  
    bone\_marrow\_adata.var.head()  
                           '''

'''  
     bone\_marrow\_adata.obs.head()  
                           '''

This is important to keep high quality cells and informative genes. Removes low quality cells, doublets and empty droplets

\# Quality control of both cells and genes

'''  
     bone\_marrow\_adata.var\_names\_make\_unique()  
                           '''

'''  
     bone\_marrow\_adata.obs\_names\_make\_unique()  
                           '''

n\_counts (total UMI / read counts per cell)

n\_genes (number of genes detected per cell; usually \>0)

% mitochondrial counts — high % suggests dying/stressed cells (mitochondrial gene names often MT- or mt-)

% ribosomal / hemoglobin — useful for blood or tissue contamination detection

Doublet scores — computationally inferred doublets (two cells captured as one)

\#\#\# Mitochondrial genes  
'''  
       sc.pl.violin( bone\_marrow\_adata, \["pct\_counts\_MT"\], jitter=0.4, multi\_panel=False)  
                                                  '''  
\#\#\# number of genes in cells  
'''  
      sc.pl.violin(  
    bone\_marrow\_adata,  
    \["n\_genes\_by\_counts"\],  
    jitter=0.4,  
    multi\_panel=False)  
                                '''  
\#\#\# Filtration of cells and genes

'''  
       sc.pp.filter\_cells(adata, min\_genes=1000)  
                              '''

'''  
      sc.pp.filter\_genes(adata, min\_cells=1000)  
                                   '''

'''  
       sc.pp.scrublet(bone\_marrow\_adata) \# doublet detection  
                                '''

\#\#\# Step 3: Noramlize  
Normalization: is to correct for differences in sequencing depth / library size so expression values are comparable across cells.

\# Normalizing to median total counts  
'''  
      sc.pp.normalize\_total()  
                           '''

\# Logarithmize the data  
'''  
   sc.pp.log1p()  
          '''

\#\#\# Step 4: Find Highly variable genes  
Highly Variable Gene Selection: are genes that show biologically meaningful variation across cells. HVGs focuses on true biological differences, removes noise, speeds up computation, as well as improve clustering accuracy.

'''  
       sc.pl.highly\_variable\_genes(bone\_marrow\_adata)  
                        '''

'''  
      sc.pp.highly\_variable\_genes(bone\_marrow\_adata)  
                     '''

\#\#\# Step 5: Dimensionality reduction: PCA and UMAP  
Dimensionality reduction: This includes principal component analysis (PCAs), which sepearted naive states from inactive states and healthy from infected states

'''  
     sc.tl.pca(bone\_marrow\_adata)  
                  '''

'''  
    sc.pl.pca\_variance\_ratio(bone\_marrow\_adata, n\_pcs=10, log=False)  
                        '''

Nonlinear Dimensionality Reduction (UMAP): creates a 2D embedding for visualization. Cluster become visible, rare population stands out

\# UMAP of Hemoglobin %cell counts

'''  
    sc.pl.umap(  
    bone\_marrow\_adata,  
    color=\["pct\_counts\_HB"\],  
    size=8)  
           '''

\#\#\# Step 6: Cluster and Annonate  
Clustering (Leiden/ Louvain): Using the neighborhood graph, algorithms like Leiden: find densely connected groups of cells, maximize modularity (network theory) and discover distinct cell populations.  
Each cluster represent cell types, states or sub population  
In leiden, the greater the resolution, the number the clusters

'''  
     sc.tl.leiden(bone\_marrow\_adata, flavor="igraph", n\_iterations=2, key\_added="leiden\_res0\_02", resolution=0.02)   
                   '''

'''  
    sc.tl.leiden(bone\_marrow\_adata, flavor="igraph", n\_iterations=2, key\_added="leiden\_res0\_5", resolution=0.5)  
                 '''

'''  
    sc.tl.leiden(bone\_marrow\_adata, flavor="igraph", n\_iterations=2, key\_added="leiden\_res2", resolution=2)  
              '''

'''  
    sc.pl.umap(  
    bone\_marrow\_adata,  
    color=\["leiden\_res0\_02", "leiden\_res0\_5", "leiden\_res2"\])  
            '''

Cell Type Annotation: Cell annotation is the process of assigning biological meaning (cell type or state) to each computationally identified cluster in your dataset.

It is the most critical step because all downstream conclusions depend on whether annotation is correct. Cell annotation combines computational methods, biological knowledge, and marker genes.

Cell annotation answers:

Are these naive or activated states?

Are these monocytes or dendritic cells?

Is this a rare stem-cell-like population?

Without annotation, the clusters are meaningless.

'''  
  import pandas as pd  
  '''

'''  
  import decoupler as dc  
  '''

\#\#\# Query Omnipath and get PanglaoDB  
'''  
   markers \= dc.op.resource(name="PanglaoDB", organism="human")  
 '''

\#\#\# Keep canonical cell type markers alone  
'''  
    markers \= markers\[markers\["canonical\_marker"\]\]  
    '''

\#\#\# Remove duplicated entries  
'''  
   markers \= markers\[\~markers.duplicated(\["cell\_type", "genesymbol"\])\]  
     '''

\#\#\# Format because dc only accepts cell\_type and genesymbol  
'''  
    markers \= markers.rename(columns={"cell\_type": "source", "genesymbol": "target"})  
     '''  
'''  
    markers \= markers\[\["source", "target"\]\]  
    '''  
'''  
    markers.head()  
    '''  
\#\#\# Cell Type Annotation  
'''  
   cluster\_annotations \= bone\_marrow\_adata\_rank\[bone\_marrow\_adata\_rank\["stat"\] \> 0\].groupby("group").head(1).set\_index("group")\["name"\].to\_dict()  
     '''

'''  
   cluster\_annotations  
   '''

'''  
   sc.pl.umap(bone\_marrow\_adata, color='cell\_type', legend\_loc='on data', title='Cell Type Annotation')  
    '''

Other visualization include:  
\#\#\# Generate the dot plot

'''  
    sc.pl.dotplot(  
    bone\_marrow\_adata,  
    var\_names=marker\_genes\_dict,   
    groupby='leiden', dendrogram=True, color\_map='Reds',  
    use\_raw=False)  
    '''  
