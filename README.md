# P-SCREEN: Plant Single-Cell Resource for Evaluation of Expression Niche
by Brenda Yu and Ben Cole

This guide will walkthrough the steps of processing scRNAseq data and merging multiple scRNAseq datasets into a unified cluster map. This pipeline was successful in processing samples that were obtained via 10x Genomics and Dropseq.

## Initial set up

### Packages and programs
Install the necessary programs you will need to run this package. Install the binaries for SRATools, BBtools, HKlibs, Bcftools, Samtools, DROP-SEQ-Tools, and Gffread. The two required Python packages are STAR and UMAP. You will also need to have R installed for the final analysis.

Clone the PSCREEN repository for the code utilized in this project.
```git clone https://github.com/byu24/at.sc.db.git```

### Directory structure
Below are the descriptions for the directory structure pre-set from the PSCREEN repository.

* Data: Store the metadata and scRNAseq whitelists here. This is also where the final analyzed objects can be stored.
* Log: Store all the log files from every step here. Log files should be produced at all steps to ensure that processing was done correctly.
* Reports: Final analysis graphs should be stored here. These will be graphs produced in R.
* Scratch: Includes all raw data and intermediate files that are worked on in each step. Subdirectories should be labeled accordingly to each sample.
* Src: Includes all the working lines of code. Original and example scripts can be found here.


## Download
1. Prior to downloading the data, a metadata CSV file must be set up. Follow the format of `sample_metadata.csv` in the `Data` directory. 
	* Set a naming structure for the samples in the first column of the metadata CSV. 
2. Follow the sample `mk_download.sh` in the `Src` directory to make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: SRATools is required to be installed beforehand. 
3. All scripts can be submitted using `launch_download.sh`.

## Filter
1. Specify the scRNAseq technique used for the sample in the metadata.CSV file. Dropseq and 10x Genomics utilizes different filtering methods.
2. Follow the sample `mk_filter.sh` in the `Src` directory to make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: BBTools is required to be installed beforehand. All scripts can be submitted using `launch_filter.sh`.

## Map
### Setting up the reference genome
Reference genomes must be created prior to mapping samples. Sample code be found under `Src` in `create_at10_genome.sh` and `create_refs.sh`. 
1. Specify the scRNAseq technique used for the sample in the metadata.CSV file. Dropseq and 10x Genomics utilizes different filtering methods.
2. A whitelist should be placed in the `Data` directory for the 10x Genomics samples.
3.
4. Follow the sample `mk_map.sh` in the `Src` directory to make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: BBTools is required to be installed beforehand. All scripts can be submitted using `launch_map.sh`.

