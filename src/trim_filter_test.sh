#!/bin/bash
#SBATCH --mail-user=brendayu@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -A gtrnd
#SBATCH -q genepool_shared
#SBATCH -J filter_dc_019
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=filter_dc_019.out

prefix=$1 ## BC 6/18/19 This argument is your friend; makes the below code reusable
genome=$2
version=$3
genomedir=$BSCRATCH/at.sc.db/scratch/${prefix}/${genome}


## Trim and filter fastq files

## BC 6/18/19: Avoid creating another folder in the scratch directory; makes it more difficult to handle later on.
##             Instead, create a "prefix" directory where you store all intermediate files: $BSCRATCH/at.sc.db/scratch/${prefix}/${prefix}_

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/dc_019_1.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/dc_019_2.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_filtered.fastq \
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_filtered.fastq \
  k=21 edist=2 ref=$BSCRATCH/at.sc.db/scratch/at10/rrna_seqs.fa \ #<-- BC20190618: Remove backslash at the end
  
$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_filtered.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_filtered.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R1_polyA.fastq \ #BC20190618: Should call output for second read, R2
# out2=$BSCRATCH/at.sc.db/scratch/filtered; #BC20190619
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R1_polyA.fastq \
  k=21 mink=11 minlen=30 \
  ref=$BSCRATCH/bin/bbmap/resources/adapters.fa \
  edist=2 ktrim=r qtrim=30 \ # <-- BC20190618: Remove backslash at the end

# BC20190618: This section is duplicating the previous step
$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_filtered.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_filtered.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R2.fastq \ # BC20190618: Should call output for second read, R2
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R2.fastq \
  k=21 mink=11 minlen=30 \
  ref=$BSCRATCH/bin/bbmap/resources/adapters.fa \
  edist=2 ktrim=r qtrim=30 \ #<-- BC20190618: Remove backslash at the end

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R1_polyA.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R1_polyA.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R1.fastq \ 
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R1.fastq \ # BC20190618: Should call ouptut for second read, R2
  ftr=$((cb_len + umi_len - 1)) \ #<-- BC20190618: remove backslash at the end

## BC 6/18/2019: Suggested code (remove sbatch arguments at the top, move all downloaded fastq files to the scratch folder.
##               Call the below script (say, "filter_reads.sh" using a separate .bs file for each dataset:
##               #!/bin/bash
##               #SBATCH HEADER
##                
##               sh src/filter_reads.sh dc_019 10x.V2

#!/bin/bash
prefix=$1
technology=$2

if [[ $technology eq "ds" ]] then
  cb_len=12
  umi_len=8
elif [[ $technology eq "10x.V2" ]] then
  cb_len=16
  umi_len=10
elif [[$technology eq "10x.V3" ]] then
  cb_len=16
  umi_len=12
fi

file_prefix=$BSCRATCH/at.sc.db/scratch/${prefix}/${file_prefix}

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=${file_prefix}_1.fastq \
  in2=${file_prefix}_2.fastq \
  out1=${file_prefix}_rRNA_filtered_1.fastq \
  out2=${file_prefix}_rRNA_filtered_2.fastq \
  k=21 edist=2 ref=$BSCRATCH/at.sc.db/scratch/at10/rrna_seqs.fa

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=${file_prefix}_rRNA_filtered_1.fastq \
  in2=${file_prefix}_rRNA_filtered_2.fastq \
  out1=${file_prefix}_rRNA_adapter_filtered_1.fastq \
  out2=${file_prefix}_R2.fastq \
  ref=$BSCRATCH/bin/bbmap/resources/adapters.fa \
  k=21 mink=11 minlen=$((cb_len + umi_len - 1)) edist=2 ktrim=r qtrim=30
