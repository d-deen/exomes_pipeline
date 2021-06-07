#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=10G


Out=$1
Job=$2

export LC_ALL=en_GB.utf-8
export LANG=en_GB.utf-8
module load MultiQC/1.7-foss-2018b-Python-3.6.6

fastq_folder=${Out}/${Job}_results/${Job}_fastqc
bwa_folder=${Out}/${Job}_results/${Job}_bam

multiqc ${fastq_folder} ${bwa_folder} -o ${Out}/${Job}_report -n ${Job}_QC_report

#multiqc ${bwa_folder}/mito_coverage_exons -o ${Out}/${Job}_report -n ${Job}_mito_coverage_exons_report
#multiqc ${bwa_folder}/mito_coverage_CDS -o ${Out}/${Job}_report -n ${Job}_mito_coverage_CDS_report

