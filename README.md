# P-SCREEN: Plant Single-Cell Resource for Evaluation of Expression Niche
by Brenda Yu and Ben Cole

This guide will walkthrough the steps of processing scRNAseq data and merging multiple scRNAseq datasets into a unified cluster map. This pipeline was successful in processing samples that were obtained via 10x Genomics and Dropseq. The *Arabidopsis thaliana* single-cell root samples used in this project were from [Denyer, et al. (2019) *Developmental Cell*](https://www.sciencedirect.com/science/article/abs/pii/S1534580719301455), [Jean-Baptiste, et al. (2019) *The Plant Cell*](http://www.plantcell.org/content/31/5/993.abstract), [Ryu, et al. (2019) *Plant physiology*](http://www.plantphysiol.org/content/179/4/1444.abstract), [Shulse, et al. (2019) *Cell reports*](https://www.sciencedirect.com/science/article/pii/S2211124719305273), and [Zhang, et al. (2019) *Molecular plant*](https://www.sciencedirect.com/science/article/pii/S1674205219301339). 

## Initial set up

### Packages and programs
Install the necessary programs you will need to run this package. Install the binaries for NCBI's [SRATools](https://github.com/ncbi/sra-tools), DOE Joint Genome Institute's [BBtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/), [Samtools/BCFTools/HTSlib](http://www.htslib.org/download/), [DROP-seq Tools](https://github.com/broadinstitute/Drop-seq), and [Gffread](https://github.com/gpertea/gffread). The two required Python packages are [STAR](https://anaconda.org/bioconda/star) and [UMAP](https://umap-learn.readthedocs.io/en/latest/). You will also need to have R installed for the final analysis.

Clone the PSCREEN repository for the code utilized in this project.
```git clone https://github.com/byu24/at.sc.db.git```

### Directory structure
Below are the descriptions for the directory structure pre-set from the PSCREEN repository.

* `data`: Store the metadata and scRNAseq whitelists here. This is also where the final analyzed objects can be stored.
* `log`: Store all the log files from every step here. Log files should be produced at all steps to ensure that processing was done correctly.
* `reports`: Final analysis graphs should be stored here. These will be graphs produced in R.
* `scratch`: Includes all raw data and intermediate files that are worked on in each step. Subdirectories should be labeled accordingly to each sample.
* `src`: Includes all the working lines of code. Original and example scripts can be found here.


## Download
1. Prior to downloading the data, a metadata CSV file must be set up. Follow the format of `sample_metadata.csv` in the `data` directory. 
	* Set a naming structure for the samples in the first column of the metadata CSV. 
2. Follow the sample `mk_download.sh` in the `src` directory to make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: SRATools is required to be installed beforehand. 
3. All scripts can be submitted using `launch_download.sh`.

## Filter
1. Specify the scRNAseq technique used for the sample in the metadata.CSV file. Dropseq and 10x Genomics utilizes different filtering methods.
2. Follow the sample `mk_filter.sh` in the `src` directory to make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: BBTools is required to be installed beforehand. All scripts can be submitted using `launch_filter.sh`.

## Map
### Setting up the reference genome
Reference genomes must be created prior to mapping samples. Sample code be found under `src` in `create_at10_genome.sh` and `create_refs.sh`. Reference genomes should be stored under `scratch` under `genomes`. 

### Mapping samples
1. Specify the scRNAseq technique used for the sample in the metadata.CSV file. Dropseq and 10x Genomics utilizes different filtering methods.
2. A whitelist should be placed in the `data` directory for the 10x Genomics samples. DropSeq generates whitelists *de novo*. Indicating the sample is Dropseq will trigger `map_dropseq.sh`. The sample code can be found within `mk_map.sh` 
3. Follow the sample `mk_map.sh` in the `src` directory to make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: BBTools, STAR, and Drop-seq Tools are required to be installed beforehand. All scripts can be submitted using `launch_map.sh`.

## Analyze
The last step is divided into three subsections: Process, Merge, and Analyze. Before proceeding, be sure to set up a Python environment with [UMAP](https://umap-learn.readthedocs.io/en/latest/). It is highly recommended to run these steps on a HPC as the computing steps require a lot of memory. 

Install the necessary R packages in the Python environment. 
```
install.packages('future', repos="https://cran.cnr.berkeley.edu/")
install.packages('furrr', repos="https://cran.cnr.berkeley.edu/")
install.packages("Seurat",repos="https://cran.cnr.berkeley.edu/")
install.packages("ggplot2",repos="https://cran.cnr.berkeley.edu/")
```

### Process
The final mapping outputs will be inputted into the Seurat Package in R.









