# Materials for 'Obtaining, analyzing, and interpreting omics data' subject (private version) - Course 2024/25 - Nebrija University

This subject belongs to the [Master in Bioinformatics](https://www.nebrija.com/programas-postgrado/master/bioinformatica/) of the Nebrija University. The repository contains practical exercises for each of the Modules of the subject. 

**Outline of the subject**

1. R basics for Data Science and Bioinformatics
2. Introduction to Omics data
3. Analysis of Transcriptomics Data
4. Analysis of Epigenomics Data
5. Analysis of Proteomics and Metabolomics Data

# Modules and exercises

## Module 1: R basics for Data Science and Bioinformatics

It contains a set of Rmarkdowns explaining the basics of R as a programming language. Due to the practical nature of this part, these Rmarkdowns are the resources for the theory of this part as well as the practical exercises to be done in class.

The objective is to learn the R basics, as it will be the main programming language used for the rest of the subject (among others). Of course, it is assumed that students already have some knowledge in programming using other languages such as Python. 

**Disclaimer:** These RMarkdown files are not intended to serve as standalone explanations (unlike some of the recommended books below). They are designed to be used as guides during class sessions.

### Exercise 1.1: Introduction to R

Topics to be discussed:

* R and RStudio installation. 
* Philosophy behind R.
* Installing packages: CRAN and BioConductor. 
* The R console for interactive calculations. 
* Scripting and non-interactive runs.
* Documentation in R.
* Bonus: working with projects in RStudio. 


### Exercise 1.2: Types of variables, loops, and vectorization

Topics to be discussed:

* Assigning names to objects
* Main data structures in R: vectors, matrices, data frames, lists, etc.
* Accessing to objects: indexing and subsetting.
* Programming basics: conditional expressions, defining functions, scoping...
* Loops and vectorization in R: `apply`-family functions, `aggregate()`...


### Exercise 1.3: Data visualization 

Topics to be discussed: 

* Basic plots using R-base: `plot()`, `abline()`, `hist()`, `par()`, `boxplot()`, changing colors, saving plots to disk. 
* Basic plots using `ggplot2`: philosophy behind `ggplot2`, main types of plots, setting colors, etc.


### Exercise 1.4: Brief introduction to `tidyverse`

Topics to be discussed: 

* Introduction to tidy data.
* Pipes in R.
* The `dplyr` R package: main functions and operations. 


### Exercise 1.5: Reading and saving data

Topics to be discussed: 

* Using R-base functions.
* The `readr` R package.
* The `data.table` R package.
* Saving complex objects to disk.


## Module 2: Introduction to Omics data

### Exercise 2.1: Getting familiar with GEO and SRA

This exercise consists of Writing a small Python script to use `sra-tools` to obtain raw data from the SRA database. 

## Module 3: Analysis of Transcriptomics Data

This part can be subdivided into three main sections: 

* Transcriptomics at the tissue level
* Transcriptomics at the single-cell level
* Spatial Transcriptomics

For each of them, exercises will consist of taking public datasets and analyzing them to get familiar with the computational pipeline required in each case, and to learn how to interpret the results. 

### _Transcriptomics at the tissue level_

### Exercise 3.1: Preprocessing data: from reads to expression matrix

Students will download a small public dataset, and use tools using the Linux command line to preprocess raw reads (`fastq` files) from a experiment. The idea is to obtain the gene expression matrix for further analyses. 

### Exercise 3.2: Filtering, normalization, and multivariate analyses

Using R and BioConductor packages, students will implement the necessary steps to analyze the data from the expression matrix generated before.

### Exercise 3.3: Differential expression and enrichment analyses

Students will learn how to perform differential expression and enrichment analyses to make biological conclusions about the data. 

### _Transcriptomics at the single-cell level_

### Exercise 3.4: Analyzing single-cell RNA-seq data using Seurat

[Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) is the most used software package in R for the analysis of single-cell RNA-seq data. In this exercise, students will download a public dataset of human Peripheral Blood Mononuclear Cells (PBMCs), and perform the analyses from normalization to visualization and clustering. 

### _Spatial Transcriptomics_

### Exercise 3.5: Deconvolution of Spatial Transcriptomics Data

First generation of sequencing-based spatial transcriptomics technologies were unable to reach single-cell resolution, meaning that each profiled spatial location usually consists of several cells (an thus, several cell types). In this exercise, students will use two deconvolution tools for spatial transcriptomics data to uncover the cell proportions underlying each spot using a dataset of mouse draining lymph nodes. 

## Module 4: Analysis of Epigenomics Data

### Exercise 4.1: ATAC-seq: peak calling and quantification

Students will process a BAM file coming from an ATAC-seq experiment to understand the main steps carried out for the analysis of this kind of data: read filtering, quality control, peak calling, and visualization of data. 

### Exercise 4:2 Peak annotation and motif enrichment analysis

After defining peaks, students will perform common analyses for ATAC-seq data: multivariate analyses, differential accessibility analysis, and motif enrichment analysis.

## Module 5: Analysis of Metabolomics data

Students will download a public metabolomics dataset already preprocessed and will analyze the abundance matrix using R. Students wil learn the particular steps carried out to analyze metabolomics data. 

# Materials and resources

This is just a list of interesting books, web pages and papers for the subject. It is not required to read/study all of them, but they may be useful to you now or in the future. 

## R, Data Science, and Statistics

* [R Programming for Data Science](https://bookdown.org/rdpeng/rprogdatascience/) by Roger D. Peng. Great resource to learn R.
* [Introduction to Data Science: Data Wrangling and Visualization with R](https://rafalab.dfci.harvard.edu/dsbook-part-1/) by Rafael Irizarry. This is a great book to start out learning R and Statistics. I strongly recommend reading it and going through its examples and exercises. 
* [Introduction to Data Science: Statistics and Prediction Algorithms Through Case Studies](https://rafalab.dfci.harvard.edu/dsbook-part-2/) by Rafael Irizarry. This is the second part of the book. It dives more into basic concepts in Statistics. 
* [Advanced R](https://adv-r.hadley.nz/) by Hadley Wickham. Only if you are interested to learn more about R as a programming language. 
* [Resources to understand linear models and basic statistics](https://rafalab.dfci.harvard.edu/pages/harvardx.html).
* [Course on linear regression](https://www.coursera.org/learn/regression-models).
* [Basic concepts of statistics clearly explained](https://www.nature.com/collections/qghhqm/pointsofsignificance)
* PCA vs NMF: check [this lecture](https://www.youtube.com/watch?v=9f4Rwt0yqr4&t=1727s) from the MIT where they explain the main differences between NMF and PCA.

## Transcriptomics at the tissue level

* [limma](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf): this is the user guide. It is very complete, shows a lot of examples and user-cases. 
* [edgeR](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf): same as above. 


## Single-cell transcriptomics

* [Best practices for single-cell analysis across modalities](https://www.nature.com/articles/s41576-023-00586-w). By Fabian J. Theis' group, this review covers the main steps and concerns about single-cell data analysis.
* [Current best practices in single-cell RNA-seq analysis: a tutorial](https://www.embopress.org/doi/full/10.15252/msb.20188746). Older version of the previous paper. Still a great resource. 
* [Seurat tutorials](https://satijalab.org/seurat/). They are a great starting point to get familiar with the analysis of single-cell transcriptomics data. 
* [Bioconductor tutorials](https://bioconductor.org/books/3.12/OSCA/index.html). In addition to Seurat, Bioconductor offers several packages to perform analyze single-cell data. In this online book, basides explaining these R packages, they cover several important aspects about each of the steps. It is a great resource if you are interested in digging more into the details of each part. 
* [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials/index.html). Scanpy is the alternative to Seurat/Bioconductor packages in Python. It is a great tool and implements the same algorithms. Useful for those of you more comfortable with Python. 
* [The impact of package selection and versioning on single-cell RNA-seq analysis](https://www.biorxiv.org/content/10.1101/2024.04.04.588111v2). Even though both software packages implement the same algorithms, there are some differences in the default hyperparameters, their implementation, and their final results. This paper points out these differences.
* [GitHub repository listing tools for single-cell RNA-seq analysis](https://github.com/seandavi/awesome-single-cell). You will find tools for almost anything!

About differential expression analysis in single-cell RNA-seq data: 

* [Confronting false discoveries in single-cell differential expression](https://www.nature.com/articles/s41467-021-25960-2)
* [Benchmarking integration of single-cell differential expression](https://www.nature.com/articles/s41467-023-37126-3)
* [A comparison of marker gene selection methods for single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03183-0)
* [Valid Post-clustering Differential Analysis for Single-Cell RNA-Seq](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30269-8)


## Spatial Transcriptomics

* [Seurat tutorial for ST data](https://satijalab.org/seurat/articles/spatial_vignette). ST version of Seurat workflow. 
* [Bioconductor tutorial for ST data](https://bookdown.org/sjcockell/ismb-tutorial-2023/). Bioconductor offers [SpatialExperiemnt](https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html), a great S4 class for analyzing ST data.
* [Scanpy tutorial for ST data](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html).
* [Museum of spatial transcriptomics](https://www.nature.com/articles/s41592-022-01409-2). Summary of the main ST technologies.

List of papers on spatial transcriptomics deconvolution tools:

* [Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution](https://www.nature.com/articles/s41592-022-01480-9)
* [A comprehensive benchmarking with practical guidelines for cellular deconvolution of spatial transcriptomics](https://www.nature.com/articles/s41467-023-37168-7)
* [Benchmarking and integration of methods for deconvoluting spatial transcriptomic data ](https://academic.oup.com/bioinformatics/article/39/1/btac805/6900924)
* [Cell2location maps fine-grained cell types in spatial transcriptomics](https://www.nature.com/articles/s41587-021-01139-4)
* [Robust decomposition of cell type mixtures in spatial transcriptomics](https://www.nature.com/articles/s41587-021-00830-w)
* [SpatialDDLS: an R package to deconvolute spatial transcriptomics data using neural networks](https://academic.oup.com/bioinformatics/article/40/2/btae072/7609037)


## Epigenomics

* [MACS3](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137): this is the most used peak caller for CHiP- and ATAC-seq. 
* In this [GitHub issue](https://github.com/macs3-project/MACS/discussions/435), there is a discussion about the parameters that are recommended when using MACS2/3 for ATAC-seq. 


## Proteomics and Metabolomics

* [Summary of main tools for proteomics data analysis](https://github.com/crazyhottommy/awesome-proteomics).
* [R package to analyze metabolomics data](https://github.com/saezlab/MetaProViz). We will use it in class. 
* [R package to preprocess MS-based metabolomics data](https://www.bioconductor.org/packages/release/bioc/html/xcms.html).
* [Basic tool to convert files into different formats](https://proteowizard.sourceforge.io/).
* [Review about tools for metabolomics data](https://www.mdpi.com/2218-1989/12/4/357).
