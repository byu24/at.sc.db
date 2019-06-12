prefix=$1
module load python3

mappingdir=scratch/${prefix}
genomedir=${mappingdir}/at10
version=v3
fa_source=ftp://ftp.ensemblgenomes.org/pub/release-43/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gtf_source=ftp://ftp.ensemblgenomes.org/pub/release-43/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.43.gtf.gz

mkdir $genomedir
wget https://www.ncbi.nlm.nih.gov/search/api/sequence/NR_141643.1/?report=fasta -O ${genomedir}/at10_rrna1.fa
wget https://www.ncbi.nlm.nih.gov/search/api/sequence/NR_141642.1/?report=fasta -O ${genomedir}/at10_rrna2.fa
wget https://www.ncbi.nlm.nih.gov/search/api/sequence/X52320.1/?report=fasta -O ${genomedir}/at10_rrna3.fa

## Create genome files
wget $fa_source -O ${genomedir}/at10.fa.gz
wget $gtf_source -O ${genomedir}/at10.gtf.gz

cat ${genomedir}/at10_rrna*.fa > ${genomedir}/rrna_seqs.fa

gunzip ${genomedir}/at10.fa.gz
gunzip ${genomedir}/at10.gtf.gz

STAR \
  --runMode genomeGenerate \
  --genomeDir $genomedir \
  --genomeFastaFiles ${genomedir}/at10.fa \
  --sjdbGTFfile ${genomedir}/at10.gtf \
  --outTmpDir ${genomedir}/_STARtmp \
  --sjdbOverhang 100 \
  --runThreadN 4 \
  -Xmx15G
