# Structural Variants (SV) and Copy Number Variants (CNV) Detection

This tutorial provides a detailed, step-by-step guide for detecting and genotyping structural variants (SVs) and copy number variations (CNVs) using the Delly tool on the High Performance Computing (HPC) environment. The workflow efficiently handles large sample sizes using multi-threading with GNU Parallel, accelerating analysis and reducing processing time. It is best suited for fewer than 100 samples, ensuring optimal resource utilization without the need for Slurm array jobs.

The tutorial covers the following steps:
1.	SV Calling: Using Delly to call structural variants (SVs) from BAM files.
2.	Merging BCF Files: Merging individual BCF files from different samples into a single file for further analysis.
3.	SV Genotyping: Genotyping the identified SVs for each sample.
4.	Merging Genotyped BCF Files: Merging the genotyped BCF files to create a unified dataset.
5.	CNV Detection: Detecting copy number variations using Delly and the previously genotyped BCF files.
6.	Merging CNV BCF Files: Combining CNV BCF files into a single merged file for further downstream analysis.
7.	CNV Genotype Calling: Genotyping CNVs for each sample.
8.	Merging CNV Genotype Files: Combining all CNV genotype files into a single file for final analysis.


Required Files:
Alignment file: sample.bam (BAM format, containing the alignment results)  

Reference genome file: reference.fasta (FASTA format, the reference genome for SV calling)  

Mappability map file: map.fa.gz 

## Step 1: SV Calling with Delly
This step uses the Delly tool to call structural variants (SVs) from the BAM files, generating individual BCF files for each sample.

```
#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**          
#SBATCH --mem-per-cpu=**   
#SBATCH --time=**
export OMP_NUM_THREADS=** # Set number of threads for parallel processing
# Define input folder path for BAM files
export BAM_DIR=/home/BAM
# Define reference genome file path
export REF_FILE=/home/reference.fasta
# Define output directory for BCF files
export BCF_DIR=/home/SV/BCF
# Define path and name for the log file
export LOG_FILE=/home/SV/log.svcall.txt
module load gcc/9.3.0                        # Load required GCC module
module load delly/1.1.6                      # Load Delly module for SV calling
# Use GNU parallel to process all BAM files in parallel
find $BAM_DIR -name "*.bam" | parallel -j ** '
  BAM_FILE={}
  echo -e "$BAM_FILE"
  name=${BAM_FILE##*/}
  base=${name%.bam}
  echo $name
  echo $base
  BCF_FILE=$BCF_DIR/$base.bcf              # Define output BCF file path
  delly call -g $REF_FILE -o $BCF_FILE $BAM_FILE 2>> $LOG_FILE ' 

#-g $REF_FILE: Specifies the reference genome file (in .fasta format).
#-o $BCF_FILE: Specifies the output file path and name where the results will be saved in BCF format.
#$BAM_FILE: The input BAM file that contains the sample alignment data.
#2>> $LOG_FILE: Appends any error messages to the specified log file $LOG_FILE. 
```

## Step 2: Merging BCF Files
This step merges the individual BCF files generated in Step 1 into a single merged BCF file.

```
# List all BCF files generated from Delly call
ls /home/SV/BCF/*.bcf > list_bcf.txt

#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**             
#SBATCH --mem-per-cpu=**    
#SBATCH --time=**
module load gcc/9.3.0
module load delly/1.1.6
delly merge -o /home/SV/sites.sv.bcf list_bcf.txt 1> log_merge.txt 2>&1  
# Merge BCF files into one
# delly merge: This command merges multiple BCF files into one. 
#-o /home/SV/sites.sv.bcf: Specifies the output file path for the merged BCF file.
#1> log_merge.txt: Redirects the standard output (normal execution logs) to log_merge.txt.
#2>&1: Redirects the standard error output to the same log file (log_merge.txt).
```
## Step 3: SV Genotyping
This step uses Delly to genotype structural variants (SVs) for each sample, using the merged sites file from Step 2.
```
#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**           
#SBATCH --mem-per-cpu=**    
#SBATCH --time=**
export OMP_NUM_THREADS=**
export BAM_DIR=/home/BAM
export REF_FILE=/home/reference.fasta
export SITES_FILE=/home/SV/sites.sv.bcf
export BCF_DIR=/home/SV/BCF_GENO
export LOG_FILE=/home/SV/log_sv_geno.txt
module load gcc/9.3.0
module load delly/1.1.6
find $BAM_DIR -name "*.bam" | parallel -j ** '
  BAM_FILE={}
  echo -e "$BAM_FILE"
  name=${BAM_FILE##*/}
  base=${name%.bam}
  echo $name
  echo $base
  BCF_FILE=$BCF_DIR/$base.geno.bcf
  delly call -g $REF_FILE -v $SITES_FILE -o $BCF_FILE $BAM_FILE 2>> $LOG_FILE ' 
#-v $SITES_FILE: Specifies the merged structural variant sites file.
```
## Step 4: Merging Genotyped BCF Files
This step merges the individual genotyped BCF files into one final merged genotyped BCF file.
```
# List all genotyped BCF files generated in Step 3
ls /home/SV/BCF_GENO/*.geno.bcf > list_geno.bcf.txt

#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**           
#SBATCH --mem-per-cpu=**   
#SBATCH --time=**
module load gcc/9.3.0
module load bcftools/1.16
bcftools merge -m id -O b -o /home/SV/merged.sv.bcf -l list_geno.bcf.txt 1>log_bcftoolsmerge.txt 2>&1         
#-l list_geno.bcf.txt: A text file containing paths to all the BCF files to be merged.
```
## Step 5: CNV Detection with Delly
This step performs Copy Number Variation (CNV) detection for each BAM file, leveraging the genotyped BCF files generated in Step 3.
```
#!/bin/bash
#SBATCH --account=**                     
#SBATCH --cpus-per-task=**                  
#SBATCH --mem-per-cpu=**                  
#SBATCH --time=**                           
export OMP_NUM_THREADS=**                
# Define the input folder path for BAM files
export BAM_DIR=/home/BAM
# Define the input folder path for genotyped BCF files
export GENO_DIR=/home/SV/BCF_GENO
# Define the reference genome and mapping file paths
export REF_FILE=/home/reference.fasta
export MAP_FILE=/home/map.fa.gz
# Define the output folder path for CNV BCF files
export BCF_DIR=/home/CNV/BCF
# Define the path and name for the log file
export LOG_FILE=/home/CNV/log_cnv.txt
module load gcc/9.3.0                     
module load delly/1.1.6                     
# Iterate through all BAM files in the input directory and process them in parallel
find $BAM_DIR -name "*.bam" | parallel -j ** '
  BAM_FILE={}
  echo -e "$BAM_FILE"
  name=${BAM_FILE##*/}
  base=${name%.bam}
  echo $name
  echo $base
  # Construct the corresponding geno file path based on the BAM file prefix
  GENO_FILE="$GENO_DIR/$base.geno.bcf"
  # Define the output file path and name with a .bcf suffix
  BCF_FILE=$BCF_DIR/$base.cnv.bcf
  delly cnv -g $REF_FILE -o $BCF_FILE -m $MAP_FILE -l $GENO_FILE $BAM_FILE 2>> $LOG_FILE ' 
#-l $GENO_FILE: Specifies the genotype file, typically the genotyped BCF file generated in Step 3.
```
## Step 6: Merging CNV BCF Files
This step merges the individual BCF files generated in Step5 into a single merged BCF file.
```
# List all CNV BCF files
ls /home/CNV/BCF/*.bcf > list_cnvbcf.txt

#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**              
#SBATCH --mem-per-cpu=**    
#SBATCH --time=**
module load gcc/9.3.0
module load delly/1.1.6
# Merge CNV BCF files into one
delly merge -e -p -o /home/CNV/sites.cnv.bcf -m 1000 -n 100000 list_cnvbcf.txt 1> log_merge.cnv.txt 2>&1
```
## Step 7: CNV Genotype Calling
This step performs CNV genotype calling using Delly, based on the previously merged CNV sites from Step 6.

```
#!/bin/bash
#SBATCH --account=**
#SBATCH --cpus-per-task=**
#SBATCH --mem-per-cpu=**
#SBATCH --time=**
export OMP_NUM_THREADS=**
# Define input file paths
export BAM_DIR=/home/BAM
export REF_FILE=/home/ reference.fasta
export MAP_FILE=/home/map.fa.gz
export SITES_FILE=/home/CNV/sites.cnv.bcf
# Define output directory for genotyped BCF files
export BCF_DIR=/home/CNV/BCF_GENO
export LOG_FILE=/home/CNV/log_geno.cnv.txt
# Load required modules
module load gcc/9.3.0
module load delly/1.1.6
find $BAM_DIR -name "*.bam" | parallel -j ** '
  BAM_FILE={}
  echo -e "$BAM_FILE"
  name=${BAM_FILE##*/}
  base=${name%.bam}
  echo $name
  echo $base
  BCF_FILE=$BCF_DIR/$base.cnv.geno.bcf
# Run Delly CNV genotype calling
  delly cnv -u -v $SITES_FILE -g $REF_FILE -m $MAP_FILE -o $BCF_FILE $BAM_FILE 2>> $LOG_FILE '
```
## Step 8: Merging CNV Genotype Files
In this step, we merge the individual CNV genotype BCF files from Step 7 into a single merged BCF file using bcftools merge.
```
# List all the CNV genotype BCF files
ls /home/CNV/BCF_GENO/*.cnv.geno.bcf > list_geno.cnvbcf.txt

#!/bin/bash
#SBATCH --account=** 
#SBATCH --cpus-per-task=** 
#SBATCH --mem-per-cpu=** 
#SBATCH --time=** 
module load gcc/9.3.0
module load bcftools/1.16
bcftools merge -m id -O b -o /home/CNV/merged.cnv.bcf -l list_geno.cnvbcf.txt 1>log_bcftoolsmerge.cnv.txt 2>&1
```
