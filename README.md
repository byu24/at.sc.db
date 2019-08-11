# P-SCREEN: Plant Single-Cell Resource for Evaluation of Expression Niche
by Brenda Yu and Ben Cole

This guide will walkthrough the steps of processing scRNAseq data and merging multiple scRNAseq datasets into a unified cluster map. This pipeline was successful in processing samples that were obtained via 10x Genomics and Dropseq.

## Initial set up

### Packages and programs
Install the necessary programs you will need to run this package. Install the binaries for SRATools, HKlibs, Bcftools, Samtools, DROP-SEQ-Tools, and Gffread. The two required Python packages are STAR and UMAP. You will also need to have R installed for the final analysis.

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

2. Follow the sample `mk_download_files.sh` to make individual file scripts that downloads each of the samples from NCBI. This script saves time from having to manually download each raw data file. Change all directory files and paths to match the location on your local machine. Note: SRATools is required to be installed beforehand. 

3. All scripts can be submitted using the `launch_download.sh` script generated from running `mk_download_files.sh`

4. If the raw data file was not originally from NCBI and was not in FASTA/FASTQ format, it is recommended to run `reshuffle.sh` from the from `BBTools` package to ensure that the raw data is in random order for the filtering step.

## Filter
1. 


