#!/bin/bash
#SBATCH -A gtrnd
#SBATCH -q genepool_shared
#SBATCH -J create_at10_genome
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=$BSCRATCH/at.sc.db/log/create_at10_genome.out

module load python3
source activate $BSCRATCH/bin/env_STAR

genomedir=$BSCRATCH/at.sc.db/scratch/at10
fa_source=ftp://ftp.ensemblgenomes.org/pub/release-43/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gff_source=ftp://ftp.ensemblgenomes.org/pub/release-43/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.43.gff3.gz

export JAVA_HOME=$BSCRATCH/bin/jdk-12.0.1
export PATH=$JAVA_HOME/bin:$PATH
PICARD='$BSCRATCH/bin/picard/picard.jar'
alias picard="java -jar $PICARD"
export PATH=$PATH:$BSCRATCH/bin/gffread/gffread
export PATH=$PATH:$BSCRATCH/bin/salmon/bin

mkdir $genomedir
wget https://www.ncbi.nlm.nih.gov/search/api/sequence/NR_141643.1/?report=fasta -O ${genomedir}/at10_rrna1.fa
wget https://www.ncbi.nlm.nih.gov/search/api/sequence/NR_141642.1/?report=fasta -O ${genomedir}/at10_rrna2.fa
wget https://www.ncbi.nlm.nih.gov/search/api/sequence/X52320.1/?report=fasta -O ${genomedir}/at10_rrna3.fa

# Create genome files
wget $fa_source -O ${genomedir}/arabidopsis.fa.gz
wget $gff_source -O ${genomedir}/at10.gff3.gz

cat ${genomedir}/at10_rrna*.fa > ${genomedir}/rrna_seqs.fa

gunzip ${genomedir}/arabidopsis.fa.gz
gunzip ${genomedir}/at10.gff3.gz

$BSCRATCH/bin/gffread -g ${genomedir}/arabidopsis.fa -T -o ${genomedir}/arabidopsis_raw.gtf ${genomedir}/at10.gff3

sed -E 's/transcript_id \"transcript:([[:alnum:]]+)([[:graph:]]+).*/transcript_id "\1\2 transcript_name "\1\2 gene_id "\1"; gene_name "\1";/' \
  ${genomedir}/arabidopsis_raw.gtf > ${genomedir}/arabidopsis.gtf

picard NormalizeFasta \
  INPUT=${genomedir}/arabidopsis.fa \
  OUTPUT=${genomedir}/at10.fa

picard CreateSequenceDictionary \
  REFERENCE=${genomedir}/at10.fa \
  OUTPUT=${genomedir}/at10.dict \
  SPECIES=Arabidopsis

$BSCRATCH/bin/Drop-seq_tools-2.3.0/FilterGtf \
  GTF=${genomedir}/arabidopsis.gtf \
  SEQUENCE_DICTIONARY=${genomedir}/at10.dict \
  OUTPUT=${genomedir}/at10.gtf

$BSCRATCH/bin/Drop-seq_tools-2.3.0/ConvertToRefFlat \
  ANNOTATIONS_FILE=${genomedir}/at10.gtf \
  SEQUENCE_DICTIONARY=${genomedir}/at10.dict \
  OUTPUT=${genomedir}/at10.refFlat

$BSCRATCH/bin/Drop-seq_tools-2.3.0/ReduceGtf \
  GTF=${genomedir}/at10.gtf \
  SEQUENCE_DICTIONARY=${genomedir}/at10.dict \
  OUTPUT=${genomedir}/at10.reduced.gtf

$BSCRATCH/bin/Drop-seq_tools-2.3.0/CreateIntervalsFiles \
  SEQUENCE_DICTIONARY=${genomedir}/at10.dict \
  REDUCED_GTF=${genomedir}/at10.reduced.gtf \
  PREFIX=at10 \
  OUTPUT=${genomedir}

STAR \
  --runMode genomeGenerate \
  --genomeDir $genomedir \
  --genomeFastaFiles ${genomedir}/at10.fa \
  --sjdbGTFfile ${genomedir}/at10.gtf \
  --outTmpDir ${genomedir}/_STARtmp \
  --sjdbOverhang 100 \
  --runThreadN 4 \
  -Xmx15G

# Create transcript files
$BSCRATCH/bin/gffread -g ${genomedir}/at10.fa -w ${genomedir}/at10_transcripts.fa ${genomedir}/at10.gtf

# Use salmon to create a transcriptome index
salmon index -t ${genomedir}/at10_transcripts.fa -i $genomedir -k 31
grep -E '>' ${genomedir}/at10_transcripts.fa | sed -E 's/>([[:alnum:]]+)([[:graph:]]+?).*/\1\2\t\1/' > ${genomedir}/at10_tgMap.tsv
grep -E AT[MC]G ${genomedir}/at10_tgMap.tsv | awk '{print $2}' > ${genomedir}/at10_ptGenes.tsv
echo "AT2G01010" >> ${genomedir}/at10_riboGenes.tsv
echo "AT2G01020" >> ${genomedir}/at10_riboGenes.tsv
echo "AT3G41768" >> ${genomedir}/at10_riboGenes.tsv
echo "AT3G41979" >> ${genomedir}/at10_riboGenes.tsv
