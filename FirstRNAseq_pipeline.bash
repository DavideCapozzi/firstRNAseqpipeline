#!/bin/bash

#Progetto in bash per analisi trascrittomica, utilizzando ambiente linux tramite WSL (ubuntu) e miniconda 
SECONDS=0 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 0 scarico i files necessari 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#link al progetto: https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=16063335
#link agli esperimenti:
#https://www.ncbi.nlm.nih.gov/sra/SRR15852399
#https://www.ncbi.nlm.nih.gov/sra/SRR15852429


#User input per scaricare i file fastq
read -p "Type the run of the experiment: " run 
sratoolkit.3.0.7-ubuntu64/bin/prefetch $run
sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump $run


#utilizzo prefetch per velocizzare il successivo download con fasterq-dump:
#sratoolkit.3.0.7-ubuntu64/bin/prefetch SRR15852399
#sratoolkit.3.0.7-ubuntu64/bin/prefetch SRR15852429

#avvio il downaload:
#sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR15852399
#sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump SRR15852429

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 1 Fastqc: analisi qualità delle read

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Ordino le directory
mkdir 0input_files 
for file in './*.fastq'; do
	mv $file 0input_files/
done 

#Lancio fastqc
fastqc 0input_files/* -o 1fastqc_output/

#Sposto i file nel mio OS Windows
for file in './1fastqc_output/*.html'; do
	cp $file /mnt/c/Users/Davide/Downloads/
done 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 2 pulizia file fastq

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Identifico adattatore contaminante
for forward in '0input_files/*_1.fastq'; do
	reverse=${forward/_1/_2}
	echo "Data analysis with fastp for file $forward and $reverse"
	fastp -i $forward -I $reverse --detect_adapter_for_pe
	rm fastp.json fastp.html
done

#Sono state identificate le stesse sequence adattatore contaminanti:

#Illumina TruSeq Adapter
#AGATCGGAAGAGCACACGTCTGAACTCCAGTCA for
#AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT rev

#Rimuovo adattatore dalla sequenza
conda deactivate Gp
conda activate cutadapt
mkdir 2cut
for forward in '0input_files/*_1.fastq'; do 
	reverse=${forward/_1/_2}
	basename_f=$(basename $forward)
	basename_r=$(basename $reverse)
	out_for=2cut/$basename_f
	out_rev=2cut/$basename_r

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o $out_for -p $out_rev  \
	--cores 10 \
	$forward $reverse
done
conda deactivate cutadapt
conda activate Gp

#Pulisco read che non soddisfano parametri:
mkdir 3trimmed
#trim_bin=/home/davide/miniconda3/envs/Gp/share/trimmomatic-0.39-2/trimmomatic.jar
out_dir="3trimmed"
for file_forward in `ls 2cut/*_1.fastq`; do
	file_reverse=${file_forward/_1./_2.}
	basename_f=$(basename $file_forward)
	basename_r=$(basename $file_reverse)	
	for_pai=$out_dir/${basename_f/_1/_1_PAI}
	rev_pai=$out_dir/${basename_r/_2/_2_PAI}
	for_unp=$out_dir/${basename_f/_1/_1_UNP}
	rev_unp=$out_dir/${basename_r/_2/_2_UNP}
	echo $for_pai
	trimmomatic PE -phred33 -threads 10 $file_forward $file_reverse $for_pai $for_unp $rev_pai $rev_unp \
	TRAILING:30 LEADING:30 HEADCROP:20 MINLEN:75 AVGQUAL:30
done

#verifico rimozione dell adattatore e stato dei file dopo trimming:
mkdir 4fastqc_output/
fastqc 3trimmed/*  -o 4fastqc_output/    

#Sposto i file nel mio OS Windows
for file in './4fastqc_output/*.html'; do
	cp $file /mnt/c/Users/Davide/Downloads/
done 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 3 allineamento: HISAT 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#prendo gli indici
#wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
#gunzip gunzip grch38_genome.tar.gz
#tar xvf grch38_genome.tar

mkdir 5hisat

#lancio hisat2
out_dir="5hisat"
for file_forward in `ls 3trimmed/*_1_PAI.fastq`; do
	file_id=${file_forward/_1_PAI.fastq/}
	file_id_basename=$(basename $file_id)
	out_file_sam=$out_dir/$file_id_basename".sam"
	for_pai=$file_forward
	rev_pai=${file_forward/_1_PAI.fastq/_2_PAI.fastq}
	for_unp=${file_forward/_1_PAI.fastq/_1_UNP.fastq}
	rev_unp=${file_forward/_1_PAI.fastq/_2_UNP.fastq}
	hisat2 -x 5hisat/grch38/genome -1 $for_pai -2 $rev_pai -U $for_unp,$rev_unp -S $out_file_sam
done

#hisat2 -x 5hisat/grch38/genome -1 trimmed/SRR15_for_pai -2 trimmed/SRR15_rev_pai -U trimmed/SRR15_for_unp,trimmed/SRR15_rev_unp -S hisat/sam_output


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# STEP 4 matrice delle conte: feature counts

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#scarico file di annotazione dei geni
#wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
#gunzip Homo_sapiens.GRCh38.110.gtf.gz

mkdir 6feature_counts
mv Homo_sapiens.GRCh38.110.gtf

genes=6feature_counts/Homo_sapiens.GRCh38.110.gtf
#lancio featureCounts
featureCounts -T 5 -p -B -t exon -g gene_id -a $genes -o 6feature_counts/count_matrix.txt 5hisat/*.sam

#featureCounts -T 5 -p -B -t exon -g gene_id -a 6feature_counts/Homo_sapiens.GRCh38.110.gtf -o feature_counts/count_matrix.txt hisat/sam_output

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."