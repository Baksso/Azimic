#!/bin/bash
#
#SBATCH --job-name="trimmomatic"
#SBATCH --error="trimmomatic.slurm.err"
#SBATCH --output="trimmomatic.slurm.out"
#SBATCH --time=2:00:00
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --partition=standard

module load trimmomatic/0.39

for R1 in *R1*
do
   R2=${R1//R1_001.fastq/R2_001.fastq}
   R1paired=${R1//.fastq/_paired.fastq.gz}
   R1unpaired=${R1//.fastq/_unpaired.fastq.gz}	
   R2paired=${R2//.fastq/_paired.fastq.gz}
   R2unpaired=${R2//.fastq/_unpaired.fastq.gz}	
   java -jar /shared/software/trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 -phred33 $R1 $R2 $R1paired $R1unpaired $R2paired $R2unpaired ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 HEADCROP:20 LEADING:30 TRAILING:25
 MINLEN:36
done
