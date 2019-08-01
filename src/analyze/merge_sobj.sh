#!/bin/bash
#SBATCH -C skylake
#SBATCH -A gtrnd
#SBATCH -q jgi_exvivo
#SBATCH -J merge_sobj
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=/global/projectb/scratch/byu24/at.sc.db/log/merge_sobj.out

module load python/3.7-anaconda-2019.07
source activate /global/projectb/scratch/bjcole/env/scRNAseq2
cd $BSCRATCH/at.sc.db/

Rscript --verbose $BSCRATCH/at.sc.db/src/analyze/merge_sobj.R >> $BSCRATCH/at.sc.db/log/merge_sobj.Rout
