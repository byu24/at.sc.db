prefix=$1
genome=$2
version=$3
genomedir=$BSCRATCH/at.sc.db/scratch/${prefix}/${genome}

## Initialize variables for barcode extraction, extract index information for 10x library
if [ $version = "v2" ]; then
    cb_len=16
    umi_len=10
fi

if [ $version = "v3" ]; then
    cb_len=16
    umi_len=12
fi

## Trim and filter fastq files

$BSCRATCH/bin/bbmap/bbduk.sh \
  in=scratch/${prefix}.fastq \
  out=scratch/filtered/${prefix}_filtered.fastq \
  k=21 edist=2 ref=$BSCRATCH/at.sc.db/scratch/at10/rrna_seqs.fa \
  -Xmx15g
  
$BSCRATCH/bin/bbmap/bbduk.sh \
  in=scratch/${prefix}_filtered.fastq \
  out1=scratch/filtered/${prefix}_R1_polyA.fastq \
  out2=scratch/filtered/${prefix}_R2.fastq \
  k=21 mink=11 minlen=30 \
  ref=$BSCRATCH/bin/bbmap/resources/adapters.fa \
  edist=2 ktrim=r qtrim=30 \
  -Xmx15g

$BSCRATCH/bin/bbmap/bbduk.sh \
  in=scratch/${prefix}_R1_polyA.fastq \
  out=scratch/filtered/${prefix}_R1.fastq \
  ftr=$((cb_len + umi_len - 1)) \
  -Xmx15g
