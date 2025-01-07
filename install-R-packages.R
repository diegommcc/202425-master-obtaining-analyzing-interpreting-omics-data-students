## packages to be installed for the subject
# this is just to save time, since some installations may take long. 
# beyond this packages, it is likely to find other dependencies throughout the 
# course.

# packages from CRAN
install.packages(
  c(
    "dplyr", "ggplot2", "BiocManager", "datasets", "stringr", 
    "remotes", "SpatialDDLS", "ggpubr", "data.table", "readr",
    "Seurat", "nycflights13", "ggvenn", "RCurl", "clustree",
    "contrast", "ggupset", "here"
  )
)


# packages from BioConductor
BiocManager::install(
  c(
    "SummarizedExperiment", "SingleCellExperiment", "SpatialExperiment",
    "limma", "edgeR", "ComplexHeatmap", "scran", "fgsea", "clusterProfiler",
    "decoupleR", "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi", "biomaRt",
    "msigdb", "GO.db", "clusterProfiler", "GenomicRanges",
    "TxDb.Hsapiens.UCSC.hg38.knownGene", "ChIPseeker", "rGREAT"
  )
)

# packages from github repositories
remotes::install_github("diegommcc/SpatialDDLS")
remotes::install_github("diegommcc/SpatialDDLSdata")
remotes::install_github("dmcable/spacexr")
remotes::install_github("saezlab/MetaProViz")

