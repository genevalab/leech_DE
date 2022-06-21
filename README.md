# leech_DE
Differential expression analysis of leech developmental stages

## Overview
1. Trimmomatic was run in paired end mode, with seed mismatches set to 2, palindrome clipping threshold set to 30 matches, simple clip threshold to 10 matches, and minimum adapter length to 1bp with the “keepBothPairs” option enabled. Leading and Trailing qualities were set to a value of 3. Sliding window size was set to 4bp, with a required quality of 15, and minimum read length was thresholded at 30bp.


```
echo "#.......fastqc initial quality analysis on $sample.....#"
fastqc -t 20 \
${sample}_R1_001.fastq.gz \
${sample}_R2_001.fastq.gz \
-o /projects/ccib/shain/fastqc/out #output directory

echo ""
echo "#......trimmomatic.......#"
java -jar /projects/f_geneva_1/programs/trimmomatic/trimmomatic-0.39.jar PE \
-threads 20 -phred33 -trimlog ${sample}_trim.log \
${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz \
${sample}_filtered.R1.fq ${sample}_filtered.unpaired.R1.fq \
${sample}_filtered.R2.fq ${sample}_filtered.unpaired.R2.fq \
ILLUMINACLIP:/projects/f_geneva_1/programs/trimmomatic/adapters/TruSeq3-PE-2.fa\:2:30:10:1:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30

echo ""
echo "#.......fastqc trimmomatic quality analysis on $sample......#"
fastqc -t 20 \
${sample}_filtered.R1.fq \
${sample}_filtered.R2.fq \
-o /projects/f_geneva_1/jody/popgen/fastqc/out

```


2. STAR was run with default parameters. The output was sorted and indexed with samtools prior to further analysis.

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/326/865/GCF_000326865.1_Helobdella_robusta_v1.0/GCF_000326865.1_Helobdella_robusta_v1.0_genomic.gtf.gz
gunzip GCF_000326865.1_Helobdella_robusta_v1.0_genomic.gtf.gz
mv GCF_000326865.1_Helobdella_robusta_v1.0_genomic.gtf H_robusta_v1.gtf
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/326/865/GCF_000326865.1_Helobdella_robusta_v1.0/GCF_000326865.1_Helobdella_robusta_v1.0_genomic.fna.gz
gunzip GCF_000326865.1_Helobdella_robusta_v1.0_genomic.fna.gz
mv GCF_000326865.1_Helobdella_robusta_v1.0_genomic.fna H_robusta_v1.fa


STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir H_robusta \
--genomeFastaFiles H_robusta_v1.fa \
--sjdbGTFfile H_robusta_v1.gtf \
--sjdbOverhang 149
```

```
STAR --genomeDir /GENOMEDIR \
--runThreadN 6 \
--readFilesIn READFILE.fq \
--outFileNamePrefix  OUTFILENAME \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

samtools sort

samtools index
```


3. Feature counts were run in non-stranded fashion to collect counts from the provided reference annotation (passed to us in Dan’s emails) to summarize features across gene features based on the corresponding gene_name or gene_id annotation in the reference annotation files (depending on analyzed species). The largest overlap option was turned on to ensure that reads were assigned to genes that they overlapped with best in the annotation.


4. While the DESeq2 code is not something I can share, please refer to the DESeq2 vignette here, which can provide you with a suitable idea of the manner of the code that was run.
