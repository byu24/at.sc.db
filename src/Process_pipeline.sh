prefix=$1


#Set Java path to the folder
export JAVA_HOME=$BSCRATCH/bin/jdk-12.0.1
export PATH=$JAVA_HOME/bin:$PATH
java -version

#Create an environmental variable
PICARD='$BSCRATCH/bin/picard/picard.jar'
alias picard="java -jar $PICARD"

#Create a sequence dictionary
picard CreateSequenceDictionary \
R=$BSCRATCH/at.sc.db/scratch/at10/at10.fa \
O=$BSCRATCH/at.sc.db/scratch/at10/at10.dict \
SP=at

#Create a refFlat File
picard ConvertToRefFlat
ANNOTATIONS_FILE=$BSCRATCH/at.sc.db/scratch/at10/at10.gtf \
SEQUENCE_DICTIONARY=$BSCRATCH/at.sc.db/scratch/at10/at10.dict \
OUTPUT=$BSCRATCH/at.sc.db/scratch/at10/at10.refFlat

picard ReduceGTF
SEQUENCE_DICTIONARY=$BSCRATCH/at.sc.db/scratch/at10/at10.dict \
GTF=$BSCRATCH/at.sc.db/scratch/at10/at10.gtf \
OUTPUT=$BSCRATCH/at.sc.db/scratch/at10/at10.reduced.gtf