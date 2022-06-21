# Leech Differential Expression Project
Differential expression analysis of leech developmental stages


## Pipeline Steps
1. Trimmomatic
2. STAR
3. Subread
4. DESeq2


## Detailed steps and code
1. Trimmomatic is first used to remove sequencing adaptors and low quality bases

Code below can be run via the command ```sbatch trim.sh SAMPLENAME```

<details><summary>trim.sh</summary>
<p>
  
  ```
                                                                                                                               
#!/bin/bash

#SBATCH --partition=p_ccib_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=fastqc                       # job name for listing in queue
#SBATCH --output=slurm-%j-%x.out
#SBATCH --mem=50G                               # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=04:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=YOUREMAIL@rutgers.edu       # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons

echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module load java                                # load any modules needed
module load FastQC

sample=$1 
  
  
echo "Bash commands for the analysis you are going to run" 
  
echo "#.......fastqc initial quality analysis on $sample.....#"
fastqc -t 20 \
${sample}_R1_001.fastq.gz \
${sample}_R2_001.fastq.gz \
-o /projects/ccib/shain/fastqc/out #output directory

echo ""
echo "#......trimmomatic.......#"
java -jar /projects/f_geneva_1/programs/trimmomatic/trimmomatic-0.39.jar PE \
-threads 20 -phred33 -trimlog ${sample}_trim.log \
${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
${sample}_filtered.R1.fq ${sample}_filtered.unpaired.R1.fq \
${sample}_filtered.R2.fq ${sample}_filtered.unpaired.R2.fq \
ILLUMINACLIP:/projects/f_geneva_1/programs/trimmomatic/adapters/TruSeq3-PE-2.fa\:2:30:10:1:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30

echo ""
echo "#.......fastqc trimmomatic quality analysis on $sample......#"
fastqc -t 20 \
${sample}_filtered.R1.fq \
${sample}_filtered.R2.fq \
-o /projects/ccib/shain/fastqc/out
  
echo "done"

```
</p>
</details>


2. STAR aligns reads from each sample to the reference genome

<details><summary>trim.sh</summary>
<p>
  
  ```
                                                                                                                               
#!/bin/bash

#SBATCH --partition=p_ccib_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=STAR                      # job name for listing in queue
#SBATCH --output=slurm-%j-%x.out
#SBATCH --mem=50G                               # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=04:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=YOUREMAIL@rutgers.edu       # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons

echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/326/865/GCF_000326865.1_Helobdella_robusta_v1.0/GCF_000326865.1_Helobdella_robusta_v1.0_genomic.gtf.gz
gunzip GCF_000326865.1_Helobdella_robusta_v1.0_genomic.gtf.gz
mv GCF_000326865.1_Helobdella_robusta_v1.0_genomic.gtf H_robusta_v1.gtf
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/326/865/GCF_000326865.1_Helobdella_robusta_v1.0/GCF_000326865.1_Helobdella_robusta_v1.0_genomic.fna.gz
gunzip GCF_000326865.1_Helobdella_robusta_v1.0_genomic.fna.gz
mv GCF_000326865.1_Helobdella_robusta_v1.0_genomic.fna H_robusta_v1.fa


STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir /projects/ccib/shain/H_robusta \
--genomeFastaFiles H_robusta_v1.fa \
--sjdbGTFfile H_robusta_v1.gtf \
--sjdbOverhang 149
```

```
STAR --genomeDir /projects/ccib/shain/H_robusta \
--runThreadN 20 \
--readFilesIn ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
--outFileNamePrefix  ${sample} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

samtools sort -@ 20 -o ${sample}_sorted.bam ${sample}.bam

samtools index ${sample}_sorted.bam

```
</p>
</details>

3. Feature counts were run in non-stranded fashion to collect counts from the provided reference annotation (passed to us in Dan’s emails) to summarize features across gene features based on the corresponding gene_name or gene_id annotation in the reference annotation files (depending on analyzed species).

```
featureCounts ${sample}_sorted.bam -a H_robusta_v1.gtf -F GTF \
-G H_robusta_v1.fa -p -T 20 --largestOverlap -s 0

4. DESeq2
