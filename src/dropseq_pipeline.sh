prefix=$1

module load python3
source activate $BSCRATCH/bin/

#Set STAR
export STAR=$BSCRATCH/bin/env_STAR
export PATH=$STAR/bin:$PATH 

#Set SamTools
export SAMTOOLS=$BSCRATCH/bin/samtools
export PATH=$SAMTOOLS/bin:$PATH 

#Set htslib
export HTSLIB=$BSCRATCH/bin/htslib
export PATH=$HTSLIB/bin:$PATH 

#Set Java path to the folder
export JAVA_HOME=$BSCRATCH/bin/jdk-12.0.1
export PATH=$JAVA_HOME/bin:$PATH 
java -version

#Set DropSeq
export DROPSEQ=$BSCRATCH/bin/dropseq
export PATH=$DROPSEQ/bin:$PATH 

#Create an environmental variable
PICARD='$BSCRATCH/bin/picard.jar'
alias picard="java -jar $PICARD"

#Create a sequence dictionary
picard CreateSequenceDictionary \
R=$BSCRATCH/at.sc.db/scratch/at10/at10.fa \
O=$BSCRATCH/at.sc.db/scratch/at10/at10.dict \
SP=at

#Create a refFlat File
java $BSCRATCH/bin/Drop-seq-2.3.0/src/java/org/broadinstitute/dropseqrna/annotation/ConvertToRefFlat.java \
ANNOTATIONS_FILE=$BSCRATCH/at.sc.db/scratch/at10/at10.gtf \
SEQUENCE_DICTIONARY=$BSCRATCH/at.sc.db/scratch/at10/at10.dict \
OUTPUT=$BSCRATCH/at.sc.db/scratch/at10/at10.refFlat

picard ReduceGTF \
SEQUENCE_DICTIONARY=$BSCRATCH/at.sc.db/scratch/at10/at10.dict \
GTF=$BSCRATCH/at.sc.db/scratch/at10/at10.gtf \
OUTPUT=$BSCRATCH/at.sc.db/scratch/at10/at10.reduced.gtf


#Run metaData Generation Pipeline
bash $BSCRATCH/bin/Drop-seq-2.3.0/src/scripts/create_Drop-seq_reference_metadata.sh \
-n at10 \        #Name for reference metadata set
-r $BSCRATCH/at.sc.db/scratch/at10/at10.fa \    #Reference dataset
-s at \     #Species name
-g $BSCRATCH/at.sc.db/scratch/at10/at10.gtf \    #gene annotation file
-o $BSCRATCH/at.sc.db/scratch/at10 \   #Output directory
-t $BSCRATCH/at.sc.db/scratch/at10 \    #Temp directory

bash $BSCRATCH/bin/dropseq/share/dropseq_tools-2.0.0-0/create_Drop-seq_reference_metadata.sh -n at10 -r $BSCRATCH/at.sc.db/scratch/at10/at10.fa -s at -g $BSCRATCH/at.sc.db/scratch/at10/at10.gtf -d $BSCRATCH/bin/dropseq -o $BSCRATCH/at.sc.db/scratch/at10 -t $BSCRATCH/at.sc.db/scratch/at10











