# Leech Differential Expression Project
Differential expression analysis of leech developmental stages

## Pipeline Steps
1. Trimmomatic
2. STAR
3. Subread
4. DESeq2

![Pipeline Steps](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/RNAseqWorkflow.png)

Image from [HBC Training tutorial](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)


## Detailed steps and code
1. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is first used to remove sequencing adaptors and low quality bases

Code below can be run via the command ```sbatch trim.sh SAMPLENAME```

<summary>trim.sh</summary>
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

2. Create Conda environment and install needed software: [STAR](https://github.com/alexdobin/STAR) and [subread](https://sourceforge.net/projects/subread/)
```
conda create --name STAR
conda activate STAR
conda install -c bioconda star 
conda install -c bioconda subread

```

3. Download and set up STAR genome files (only needs to be run once ever)

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

4. Align reads from each sample to the reference genome. Needs to be run for each sample. Execute by running ```sbatch run_star.sh SAMPLE_NAME```

<summary>run_star.sh</summary>
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
eval "$(conda shell.bash hook)"
conda activate STAR

sample=$1 

STAR --genomeDir /projects/ccib/shain/H_robusta \
--runThreadN 20 \
--readFilesIn ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
--outFileNamePrefix  ${sample} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

#samtools sort -@ 20 -o ${sample}_sorted.bam ${sample}.bam

#samtools index ${sample}_sorted.bam

```
</p>

4. [Featurecounts](https://sourceforge.net/projects/subread/) counts the number of reads per gene

```
#!/bin/bash

#SBATCH --partition=p_ccib_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=featCount                      # job name for listing in queue
#SBATCH --output=slurm-%j-%x.out
#SBATCH --mem=50G                               # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=04:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=YOUREMAIL@rutgers.edu       # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons

module purge
eval "$(conda shell.bash hook)"
conda activate STAR

sample=$1 

featureCounts ${sample}Aligned.sortedByCoord.out.bam -a H_robusta_v1.gtf -F GTF \
-G H_robusta_v1.fa -p -T 20 --largestOverlap -s 0
```

5. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
*working on this based on [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

Once all featureCounts has been run for all samples we need to combine the counts into a single count matrix before inputting them to DESeq2. This can be accomplished via text manipulation on the command line:
```
for i in *counts.txt; do cut -f7 ${i} > ${i}_clean.txt ; done
ls -1  *counts.txt | head -1 | xargs cut -f1 > genes.txt
paste genes.txt *counts.txt_clean.txt > All_sample_count_matrix.txt
```

