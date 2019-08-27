# P-SCREEN: Plant Single-Cell Resource for Evaluation of Expression Niche
<p align="center">
	by Brenda Yu<sup>1</sup>, Ben Cole<sup>2</sup>, and Axel Visel<sup>2</sup><br>
	Last updated: <b>August 27, 2019</b>
	<sup>1</sup>University of California, Merced, <sup>2</sup>Department of Energy Joint Genome Institute, Walnut Creek, CA<br>
	</p>

This guide will walkthrough the steps of processing scRNAseq data and merging multiple scRNAseq datasets into a unified cluster map. This pipeline was successful in processing samples that were obtained via 10x Genomics and Dropseq. Four main sections are detailed below: Download, Filter, Map, and Process. The *Arabidopsis thaliana* single-cell root samples used in this project were from [Denyer, et al. (2019) *Developmental Cell*](https://www.sciencedirect.com/science/article/abs/pii/S1534580719301455), [Jean-Baptiste, et al. (2019) *The Plant Cell*](http://www.plantcell.org/content/31/5/993.abstract), [Ryu, et al. (2019) *Plant physiology*](http://www.plantphysiol.org/content/179/4/1444.abstract), [Shulse, et al. (2019) *Cell reports*](https://www.sciencedirect.com/science/article/pii/S2211124719305273), and [Zhang, et al. (2019) *Molecular plant*](https://www.sciencedirect.com/science/article/pii/S1674205219301339). 

## Initial set up

### Packages and programs
Install the necessary programs you will need to run this package. Install the binaries for NCBI's [SRATools](https://github.com/ncbi/sra-tools), DOE Joint Genome Institute's [BBtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/), [Samtools/BCFTools/HTSlib](http://www.htslib.org/download/), [DROP-seq Tools](https://github.com/broadinstitute/Drop-seq), and [Gffread](https://github.com/gpertea/gffread). The two required Python packages are [STAR](https://anaconda.org/bioconda/star) and [UMAP](https://umap-learn.readthedocs.io/en/latest/). You will also need to have R installed for the final analysis.

Clone the PSCREEN repository for the code utilized in this project.
```
git clone https://github.com/byu24/at.sc.db.git
```
### Directory structure
Below are the descriptions for the directory structure pre-set from the PSCREEN repository.

* `data`: Store the metadata and scRNAseq whitelists here. This is also where the final analyzed objects can be stored.
* `log`: Store all the log files from every step here. Log files should be produced at all steps to ensure that processing was done correctly.
* `reports`: Final analysis graphs should be stored here. These will be graphs produced in R.
* `scratch`: Includes all raw data and intermediate files that are worked on in each step. Subdirectories should be labeled accordingly to each sample.
* `src`: Includes all the working lines of code. Original and example scripts can be found here. 


## Download
1. Prior to downloading the data, a metadata CSV file must be set up. Follow the format of `sample_metadata.csv` in the `data` directory. 
	* Set a naming structure for the samples in the first column of the metadata CSV. Use an easy structure that will help identify the originating sample in the final database.
	* Include the SRA IDs and the corresponding SRR IDs in the appropriate column. If a sample has two or more SRR IDs under the same SRA, samples must be concatenated after filtering.
	* The scRNAseq Platform and Version is important to specify. The two methods used in this project were DropSeq and 10xGenomics.
	* If the samples have been sequenced with additional genomes, specify the genome in the "Genome" column. This project had "at10" (*A. thaliana* only), "at10_hg38" (with human genome), and "at10_mm10" (with mouse genome). Specifying the genomes is necessary for accurate mapping.
	* Specify the number of expected cells in the "ExpCell" column. This information is used in the Analyze section.
2. Running `mk_download.sh` in the `src/download` directory will make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: SRATools is required to be installed beforehand. 
3. All scripts can be submitted using `launch_download.sh`.

## Filter
1. Specify the scRNAseq technique used for the sample in the metadata.CSV file. Dropseq and 10x Genomics utilizes different filtering methods.
2. Running `mk_filter.sh` in the `src/filter` directory will make individual file scripts that filters each sample. This script can be programmed to accept 10xGenomics or DropSeq samples. Change all directory files and paths to match the location on your local machine. Note: BBTools is required to be installed beforehand. All scripts can be submitted using `launch_filter.sh`.

## Map
### Setting up the reference genome
Reference genomes must be created prior to mapping samples. Sample code be found under `src` in `create_at10_genome.sh` and `create_refs.sh`. Reference genomes should be stored under `scratch` under `genomes`. 
### Mapping samples
1. Specify the scRNAseq technique used for the sample in the metadata.CSV file. Dropseq and 10x Genomics utilizes different filtering methods.
2. A whitelist should be placed in the `data` directory for the 10x Genomics samples. DropSeq generates whitelists *de novo*. Indicating the sample is Dropseq will trigger `map_dropseq.sh`. The sample code can be found within `mk_map.sh` 
3. Running `mk_map.sh` in the `src/map` directory will make individual file scripts that maps each sample. Change all directory files and paths to match the location on your local machine. Note: BBTools, STAR, and Drop-seq Tools are required to be installed beforehand. All scripts can be submitted using `launch_map.sh`.


## Analyze
The last step is divided into three subsections: Process, Merge, and Analyze. Before proceeding, be sure to set up a Python environment with [UMAP](https://umap-learn.readthedocs.io/en/latest/). It is highly recommended to run these steps on a HPC as the computing steps require a lot of memory. 

Install the necessary R packages in the Python environment. Additional features and examples can be explored for the [Seurat Package](https://satijalab.org/seurat/vignettes.html) on its website.
```
install.packages('future', repos="https://cran.cnr.berkeley.edu/")
install.packages('furrr', repos="https://cran.cnr.berkeley.edu/")
install.packages("ggplot2",repos="https://cran.cnr.berkeley.edu/")
install.packages("Seurat",repos="https://cran.cnr.berkeley.edu/")
```
### Process
This subsection processes the mapping outputs for each sample into a Seurat object. 
1. `process_cluster.r` in the `src/analyze` directory includes all the functions needed to create the initial Seurat object. Edit the functions to reflect the path to the sample mapping matrix directory. Source 'process_cluster.r` at the beginning of your R script for each sample.
```
source('"/path/to/at.sc.db/src/analyze/process_cluster.r"')
```
2. `process_ici.R` in the `src/analyze` directory includes all the functions needed to calculate ICI scores [(Elfroni, et al. (2016))](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4354993/) from the Seurat object. Source the functions at the beginning of your R script for each sample.
```
source('"/path/to/at.sc.db/src/analyze/process_ici.R"')
```
3. Run `mk_process.sh` in the `src/analyze` directory to make sample-specific R scripts and the corresponding batch script. This script will create the Seurat object from the mapping matrices, save an Elbow Plot and UMAP dimension plot, and calculate the ICI scores. The output will be a sample Seurat Robject and the corresponding ICI Score Robject. Be sure that the output directory has enough space and memory to accommodate the Robjects.
### Merge
This subsection merges all the individual samples and ICI scores together into a master Seurat object. The scripts on this step could be integrated into the Process subsection but it was separated due to the high memory and space allocations required.
1. Edit `merge.R` in the `src/analyze` directory to load all the Seurat objects into the R environment. The script will merge all the individual samples into a master Seurat object.
2. Edit `merge.R` in the `src/analyze` directory to load all the ICI score Robjects into the R environment. A function exists in the script to calculate the sum of the sample ICI scores and convert them into a dataframe. The ICI scores will then be merged with the master Seurat object.
3. It is recommended to run `mk_merge.sh` and submit it as a batch job. Be sure to submit onto a HPC cluster with a large memory and space allocation. 
4. The output will be the final Seurat Object with the appropriate statistical calculations completed. Additional statistical calculations can be added to the script by reviewing example Seurat [vignettes](https://satijalab.org/seurat/vignettes.html). The log file produced will also output the metadata information to review and confirm the accuracy of the merging.
### Analyze
This subsection features data manipulation to produce different statistical plots in `analyze.R`. Plotting can be done on a local machine if desired or submitted to a HPC using `mk_analyze.sh`.
1. A `DimPlot` with a `PCA` reduction will be outputted to visualize all the clusters in the master Seurat object.
2. A section to isolate only QC type cells is included and can be customized to isolate other cell types.
3. Individual cluster composition of cell types can be visualized as barplots using `ggplot()`. Clusters can be reordered and grouped by cell type composition. The steps shown in `analyze.R` were completed manually. Additional scripts can be used to reclassify the clusters.
4. If clusters are reclassified, a new corresponding cluster map should be generated using `DimPlot`. 









