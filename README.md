# P-SCREEN: Plant Single-Cell Resource for Evaluation of Expression Niche
by Brenda Yu and Ben Cole

This guide will walkthrough the steps of processing scRNAseq data and merging multiple scRNAseq datasets into a unified cluster map. This pipeline was successful in processing samples that were obtained via 10x Genomics and Dropseq.

## Initial set up

### Packages and programs
Install the necessary programs you will need to run this package. Install the binaries for SRATools, HKlibs, Bcftools, Samtools, DROP-SEQ-Tools, and Gffread. The two required Python packages are STAR and UMAP. 

Clone the PSCREEN repository for the code utilized in this project.
```git clone https://github.com/byu24/at.sc.db.git```

### Directory structure
Below are the descriptions for the directory structure pre-set from the PSCREEN repository.

Data: Store the metadata and scRNAseq whitelists here. This is also where the final analyzed objects can be stored.



## Download

1. Set up the metadata 


