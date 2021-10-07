#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=80G
#SBATCH -c 8


Out=$1
Job=$2

##Module loading
module load Java/11.0.2

##Setting constants for variant calling:
gatk_path='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/gatk-4.1.9.0'

fasta_ref="/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/genomes_DD/GCA_000001405.15_GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

bwa_folder=${Out}/${Job}_results/${Job}_bam
vcf_folder=${Out}/${Job}_results/${Job}_vcf

samplesheet=${Out}/samplesheet.txt
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'` 

##Running GATk to get gVCF for every sample
${gatk_path}/gatk --java-options "-Xmx80g" HaplotypeCaller  \
   -R ${fasta_ref} \
   -I ${bwa_folder}/${name}_bwa_dedup.bam  \
   -O "${Out}/${Job}_results/${Job}_vcf/${name}.g.vcf.gz" \
   -ERC GVCF \
   --do-not-run-physical-phasing true

