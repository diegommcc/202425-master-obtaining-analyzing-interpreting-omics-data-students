# Installing software needed for the subject

## R and RStudio 

The details to install R and RStudio can be found at 01-R-basics/01-Intro-R.Rmd

In addition, the R packages needed during the course can be installed by running the code in install-R-packages.R. 

## The Integrative Genomics Viewer (IGV) 

There is a web version, but in my experience it does not work pretty well. I recommend installing it, it is available for all OS: <https://igv.org/download/html/oldtempfixForDownload.html>. 

## miniconda environments

Below, you can find a list with all the software that you will need to have installed in you computer. Ideally, we should be using Docker containers, but instead we will use miniconda, a lightweight distribution of the popular Anaconda package and environment management system. It provides a minimal installation of Conda, which allows us to manage environments and packages. 

Therefore, we need fist to install `miniconda`. Follow [these instructions](https://docs.anaconda.com/miniconda/install/) according to your OP.

Once `miniconda` is installed, we are going to create different environments in which the software needed for each practical exercise will be installed and available. 

### `rna-seq-env` environment

Run the following commands: 

```bash
conda activate
conda create -n 'rna-seq-env'
conda activate rna-seq-env

conda install bioconda::samtools
conda install bioconda::fastqc
conda install bioconda::trim-galore
conda install bioconda::multiqc
# not available on M1 macos computers
conda install bioconda::rsem 
conda install bioconda::star
conda install bioconda::qualimap

# for exiting the environment
conda deactivate
```


### `sra-tools-env` environment

Make sure that the `sra-tools` version that you install is the latests (v3.11 at least).

```bash
conda activate
conda create -n 'sra-tools-env'
conda activate sra-tools-env

conda install bioconda::sra-tools
conda install anaconda::pandas
## for exiting the environment
conda deactivate
```

### `atac-seq-env`

Install `macs3` before `multiqc` to avoid a problem with Python versions. 

```bash 
conda activate
conda create -n 'atac-seq-env'
conda activate atac-seq-env

conda install bioconda::samtools
conda install bioconda::macs3
conda install bioconda::multiqc
conda install bioconda::picard # not available on M1 macos computers
conda install bioconda::deeptools
conda install bioconda::subread
conda install bioconda::homer
```

It is likely that when `homer` is used, you`ll need to install the genome version used in the analysis. Try to run it without doing the following, but if it raises an error, run: 

```bash
## you will have a different path
perl /home/dmananesc/miniconda3/envs/atac-seq-env/share/homer/.//configureHomer.pl -install hg19
```

