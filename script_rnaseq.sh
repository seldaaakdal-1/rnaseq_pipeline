#!/bin/bash

reads="/mnt/c/Users/User/Desktop/sc/reads"
fastqc="/mnt/c/Users/User/Desktop/sc/fastqc"
counts="/mnt/c/Users/User/Desktop/sc/counts"
mapped="/mnt/c/Users/User/Desktop/sc/mapped"
trimmed="/mnt/c/Users/User/Desktop/sc/trimmed"
multiqc="/mnt/c/Users/User/Desktop/sc/multiqc"
genome_index="/mnt/c/Users/User/Desktop/sc/genome_index"

echo "Step 1: Prefetch SRA Data"

prefetch -O ${reads}/ PRJNA1090488


echo "Step 2: fasterq-dump"

fasterq-dump ${reads}/*.sra -O ${reads}/


echo "Step 3: Fastqc Quality Control"

fastqc ${reads}/*.fastq -O ${fastqc}/


echo "Step 4: Multiqc"

multiqc ${fastqc}/ -o ${multiqc}/


echo "Step 5: Trimmomatic"

for infile in ${reads}/*_1.fastq
do 
	base=$(basename ${infile} _1.fastq)
	java -jar /mnt/c/Users/User/Desktop/software/trimmomatic-0.40.jar PE \
	${infile} ${reads}/${base}_2.fastq \
	${trimmed}/${base}_1.trim.fastq ${trimmed}/${base}_1.untrim.fastq \
	${trimmed}/${base}_2.trim.fastq ${trimmed}/${base}_2.untrim.fastq \
	SLIDINGWINDOW:4:15 MINLEN:36 \
	ILLUMINACLIP:/mnt/c/Users/User/Desktop/software/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 


done


echo "Step 6: STAR Genome Index"

STAR --runMode genomeGenerate \
     --genomeDir ${genome_index}/ \
     --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     --sjdbGTFfile Saccharomyces_cerevisiae.R64-1-1.115.gtf \
     --runThreadN 2


echo "Step 7: STAR Mapping"

for infile in ${trimmed}/*_1.trim.fastq
do
	base=$(basename ${infile} _1.trim.fastq)
	STAR --genomeDir ${genome_index}/ \
		 --runThreadN 2 \
		 --readFilesIn ${infile} ${trimmed}/${base}_2.trim.fastq \
		 --outFileNamePrefix ${mapped}/${base} \
		 --outSAMtype BAM SortedByCoordinate \
		 --quantMode GeneCounts \
		 --outSAMattributes Standard

done


echo "Step 8: Calculating counts with featureCounts"

featureCounts -a Saccharomyces_cerevisiae.R64-1-1.115.gtf -o ${counts}/gene_count.txt -T 2 -p ${mapped}/*.bam









