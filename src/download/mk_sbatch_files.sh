cd $BSCRATCH/at.sc.db/

srr_names=($(awk -F, '{if(NR>1) print $9}' data/sample_metadata.csv))
set_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))

len=${#srr_names[@]}
i=0
while [ "$i" -lt "$len" ]; do
  set_name=${set_names[$i]}
  srr_name=${srr_names[$i]}

  echo "#!/bin/bash" > src/download/download_${set_name}.bs
  echo "#SBATCH --mail-user=bjcole@lbl.gov" >> src/download/download_${set_name}.bs
  echo "#SBATCH --mail-type=ALL" >> src/download/download_${set_name}.bs
  echo "#SBATCH -A gtrnd" >> src/download/download_${set_name}.bs
  echo "#SBATCH -q genepool_shared" >> src/download/download_${set_name}.bs
  echo "#SBATCH -J download_SRA_${set_name}" >> src/download/download_${set_name}.bs
  echo "#SBATCH -t 12:00:00" >> src/download/download_${set_name}.bs
  echo "#SBATCH --mem-per-cpu=4000" >> src/download/download_${set_name}.bs
  echo "#SBATCH --ntasks=1" >> src/download/download_${set_name}.bs
  echo "#SBATCH --cpus-per-task=4" >> src/download/download_${set_name}.bs
  echo "#SBATCH --output=download_SRA_${set_name}.out" >> src/download/download_${set_name}.bs
  echo "" >> src/download/download_${set_name}.bs
  echo "module load python3" >> src/download/download_${set_name}.bs
  echo "source activate $BSCRATCH/env_STARsolo" >> src/download/download_${set_name}.bs

  echo "cd $BSCRATCH/at.sc.db/" >> src/download/download_${set_name}.bs

  if ! [[ $set_name =~ js|zc ]]; then 
    echo "$BSCRATCH/sratoolkit.2.9.6-1-ubuntu64/bin/prefetch.2.9.3 ${srr_name} -o scratch/${set_name} -t http" >> src/download/download_${set_name}.bs
    echo "$BSCRATCH/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump scratch/${set_name} \\" >> src/download/download_${set_name}.bs
    echo "  -o ${set_name}.fastq \\" >> src/download/download_${set_name}.bs
    echo "  -O $BSCRATCH/at.sc.db/scratch \\" >> src/download/download_${set_name}.bs
    echo "  -t $BSCRATCH/at.sc.db/scratch \\" >> src/download/download_${set_name}.bs
    echo "  -sS --include-technical" >> src/download/download_${set_name}.bs
  fi

  if [[ $set_name =~ "jsh_016" ]]; then
    echo "wget ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR808/SRR8086586/whole_root_Heatshock_possorted_genome_bam.bam \\" >> src/download/download_${set_name}.bs
  fi

  if [[ $set_name =~ "js_017" ]]; then
    echo "wget ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR808/SRR8086585/whole_root_Control_2_possorted_genome_bam.bam \\" >> src/download/download_${set_name}.bs
  fi

  if [[ $set_name =~ "js_018" ]]; then
    echo "wget ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR808/SRR8086584/whole_root_Control_1_possorted_genome_bam.bam \\" >> src/download/download_${set_name}.bs
  fi

  if [[ $set_name =~ "js" ]]; then
    echo "  -q -O- | samtools view | \\" >> src/download/download_${set_name}.bs
    echo "    awk '{ \\" >> src/download/download_${set_name}.bs
    echo -e "    match(\$0, \"CR:Z:([ATCGN]+)\", arr1)" >> src/download/download_${set_name}.bs
    echo -e "    match(\$0, \"UR:Z:([ATCGN]+)\", arr2)" >> src/download/download_${set_name}.bs
    echo -e "    match(\$0, \"CY:Z:([[:graph:]]+)\", arr3)" >> src/download/download_${set_name}.bs
    echo -e "    match(\$0, \"UY:Z:([[:graph:]]+)\", arr4)" >> src/download/download_${set_name}.bs
    echo -E "    print \$1\"\\n\" arr1[1] arr2[1] \"\\n+\\n\" arr3[1] arr4[1] > \"scratch/${set_name}_1.fastq\""  >> src/download/download_${set_name}.bs
    echo -E "    print \$1\"\\n\"\$10\"\\n+\\n\"\$11 > \"scratch/${set_name}_2.fastq\"}'"  >> src/download/download_${set_name}.bs
  fi

  if [[ $set_name =~ "zc" ]]; then
    echo "wget ftp://download.big.ac.cn/gsa/CRA001559/CRR054647/CRR054647_f1.tar.gz -O scratch/${set_name}_1.tar.gz" >> src/download/download_${set_name}.bs
    echo "tar -xzvf scratch/${set_name}_1.tar.gz -C scratch/" >> src/download/download_${set_name}.bs
    echo "zcat scratch/Root_R1.fastq.gz > scratch/${set_name}_1.fastq" >> src/download/download_${set_name}.bs
    echo "wget ftp://download.big.ac.cn/gsa/CRA001559/CRR054647/CRR054647_r2.fastq.gz -O scratch/${set_name}_2.fastq.gz" >> src/download/download_${set_name}.bs
    echo "zcat scratch/${set_name}_2.fastq.gz > scratch/${set_name}_2.fastq" >> src/download/download_${set_name}.bs
  fi

  i=$(($i + 1))
done

awk -F, '{if(NR > 1) print "sbatch src/download/download_"$1".bs &"}' data/sample_metadata.csv > src/launch_download.sh