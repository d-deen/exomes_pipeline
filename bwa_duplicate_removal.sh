#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=250G
#SBATCH -c 8


Out=$1
Job=$2

##Module loading
module load SAMtools/1.12-GCC-10.2.0
module load BWA/0.7.17-foss-2017b
module load picard/2.2.4-intel-2017.03-GCC-6.3-Java-1.8.0_144

##Setting constants
#For running the alignment
threads=$SLURM_JOB_CPUS_PER_NODE
bwa_folder=${Out}/${Job}_results/${Job}_bam

index="/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/genomes_DD/GCA_000001405.15_GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
fasta_ref="/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/genomes_DD/GCA_000001405.15_GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

#For coverages calculations
mito_exons_bed='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/custom_ann_DD/mitocarta3.interval_list'
mito_cds_bed='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/custom_ann_DD/mitocarta3_CDS.interval_list'
twist_bed='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/custom_ann_DD/twist_compr.interval_list'


##Getting the name for the sample, r1 and r2
cd ${Out}/${Job}_rawfiles

samplesheet=${Out}/samplesheet.txt
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'` 
r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $2}'` 
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $3}'`


##Generate sorted and indexed bam outputs 
bwa mem ${index} $r1 $r2 -M -t ${threads} -R "@RG\tID:FlowCell.${name}\tSM:${name}\tPL:illumina\tLB:${Job}.${name}"   | \
samtools sort - -O bam | tee ${bwa_folder}/${name}_bwa_output.bam | \
samtools index - ${bwa_folder}/${name}_bwa_output.bam.bai 2>> ${Out}/${Job}_logs/${Job}_${name}_samtools.log

##Generating statistics for bam files (aligned reads and distribution by chromosomes)
samtools flagstat ${bwa_folder}/${name}_bwa_output.bam > ${bwa_folder}/${name}_stat_bwa_output.txt 
samtools idxstats ${bwa_folder}/${name}_bwa_output.bam > ${bwa_folder}/${name}_idxstat_bwa_output.txt

echo "The reads have been aligned"


##Removing duplicates

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=${bwa_folder}/${name}_bwa_output.bam \
      O=${bwa_folder}/${name}_bwa_dedup.bam \
      REMOVE_DUPLICATES=true \
      M=${bwa_folder}/${name}_dup_metrics.txt

samtools index ${bwa_folder}/${name}_bwa_dedup.bam 

echo "The duplicates reads have been removed"

##Calculating coverages

#Coverage for TWIST
#java -jar $EBROOTPICARD/picard.jar CollectWgsMetrics \
      #I=${bwa_folder}/${name}_bwa_dedup.bam \
      #O=${bwa_folder}/${name}_twist_metrics.txt \
      #R=${fasta_ref} \
      #COUNT_UNPAIRED=true \
      #INTERVALS=${twist_bed}


#Coverage for MitoCarta exons (coding and non-coding)
#java -jar $EBROOTPICARD/picard.jar CollectWgsMetrics \
      #I=${bwa_folder}/${name}_bwa_dedup.bam \
      #O=${bwa_folder}/mito_coverage_exons/${name}_exons_MTC_metrics.txt \
      #R=${fasta_ref} \
      #COUNT_UNPAIRED=true \
      #INTERVALS=${mito_exons_bed} 

#Coverage for MitoCarta CODING exons
java -jar $EBROOTPICARD/picard.jar CollectWgsMetrics \
      I=${bwa_folder}/${name}_bwa_dedup.bam \
      O=${bwa_folder}/mito_coverage_CDS/${name}_cds_MTC_metrics.txt \
      R=${fasta_ref} \
      COUNT_UNPAIRED=true \
      INTERVALS=${mito_cds_bed} 
     

echo "The coverage has been calculated"


##Checking bam files for contamination
verifybamid_path=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD
verifybamid_prefix=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/VerifyBamID/resource/hgdp.100k.b38.vcf.gz.dat

${verifybamid_path}/VerifyBamID.Linux.x86-64 \
--BamFile ${bwa_folder}/${name}_bwa_dedup.bam \
--SVDPrefix ${verifybamid_prefix}  \
--Reference ${fasta_ref} \
--Output ${bwa_folder}/${name}

echo "VerifyBamId verified a bam file"

#Exit 0 was inserted as sometimes verifybamid fails as not enough coverage to estimate the contamination

exit 0

