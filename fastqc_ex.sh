#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=50G

#Arguments from the command line input by user 
Out=$1
Job=$2

#Module loading and setting up constants
module load FastQC/0.11.8-Java-1.8.0_144 
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.12-GCC-10.2.0

fastq_path="/nobackup/proj/ghwtcmr/Programs_DD/FastQ-Screen-0.14.1"


#Define variables
sample_list=${Out}/sample_list.txt
fastq_folder=${Out}/${Job}_results/${Job}_fastqc

sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $sample_list |  awk '{print $1}'` 

##Run fastqc
fastqc -t 8 ${sample} --outdir=${fastq_folder} 2>> ${Out}/${Job}_logs/${Job}_${name}_fastqc.log

echo "Fastqc analysis done"

##Run fastq_screen
${fastq_path}/fastq_screen ${sample} --aligner bowtie2 --force \
--threads 8 --outdir ${fastq_folder} 2>> ${Out}/${Job}_logs/${Job}_fastq_screen.log

echo "Fastq_screen analysis done"