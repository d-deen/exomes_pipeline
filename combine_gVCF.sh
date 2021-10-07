#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=80G
#SBATCH -c 8


Out=$1
Job=$2


##Module loading
module load Java/11.0.2
module load picard/2.2.4-intel-2017.03-GCC-6.3-Java-1.8.0_144

##Setting constants for variant calling:
gatk_path='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/gatk-4.1.9.0'

fasta_ref="/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/genomes_DD/GCA_000001405.15_GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


##Removing bam files with the duplicates

bwa_folder=${Out}/${Job}_results/${Job}_bam
rm ${bwa_folder}/*_output.bam*


##Combinign gVCFs into a single vcf and calling variants:

#Getting a list of gVCF files
ls -d ${Out}/${Job}_results/${Job}_vcf/*.g.vcf.gz | awk '{print "--variant" OFS $0}'  > ${Out}/${Job}_results/${Job}_vcf/gvcf_list.txt

gvcf_list=${Out}/${Job}_results/${Job}_vcf/gvcf_list.txt



#Checking how many samples are there: if more than one, then combine them in a file, if one - rename it with a cohort. 

if test "$(cat $gvcf_list | wc -l)" -gt 1 ; then
    echo "proceed to combining the files"
    ${gatk_path}/gatk CombineGVCFs \
     -R ${fasta_ref} \
    --arguments_file ${gvcf_list} \
     -O ${Out}/${Job}_results/${Job}_vcf/${Job}_cohort.g.vcf.gz
else
    mv ${Out}/${Job}_results/${Job}_vcf/*.g.vcf.gz ${Out}/${Job}_results/${Job}_vcf/${Job}_cohort.g.vcf.gz
    ${gatk_path}/gatk IndexFeatureFile \
    -I ${Out}/${Job}_results/${Job}_vcf/${Job}_cohort.g.vcf.gz
fi




#Calling genotypes

${gatk_path}/gatk --java-options "-Xmx80g" GenotypeGVCFs \
   -R ${fasta_ref} \
   -V ${Out}/${Job}_results/${Job}_vcf/${Job}_cohort.g.vcf.gz \
   -O ${Out}/${Job}_results/${Job}_vcf/${Job}_cohort.vcf.gz 


##Filtering the variants by quality
#First, separating them into SNPs and indels

${gatk_path}/gatk SelectVariants \
-V ${Out}/${Job}_results/${Job}_vcf/${Job}_cohort.vcf.gz \
    -select-type SNP \
    -O  ${Out}/${Job}_results/${Job}_vcf/${Job}_snps.vcf.gz 


${gatk_path}/gatk SelectVariants \
-V ${Out}/${Job}_results/${Job}_vcf/${Job}_cohort.vcf.gz \
    -select-type INDEL \
    -O  ${Out}/${Job}_results/${Job}_vcf/${Job}_indels.vcf.gz


#Hard filtering of the variants 

${gatk_path}/gatk VariantFiltration \
    -V ${Out}/${Job}_results/${Job}_vcf/${Job}_snps.vcf.gz  \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${Out}/${Job}_results/${Job}_vcf/${Job}_snps_filtered.vcf.gz 


${gatk_path}/gatk VariantFiltration \
    -V ${Out}/${Job}_results/${Job}_vcf/${Job}_indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"  \
    -O ${Out}/${Job}_results/${Job}_vcf/${Job}_indels_filtered.vcf.gz

#Merging the variants
java -jar $EBROOTPICARD/picard.jar MergeVcfs \
          I=${Out}/${Job}_results/${Job}_vcf/${Job}_snps_filtered.vcf.gz \
          I=${Out}/${Job}_results/${Job}_vcf/${Job}_indels_filtered.vcf.gz \
          O=${Out}/${Job}_results/${Job}_vcf/${Job}_filtered_comb_cohort.vcf.gz


##Normalisation of the variants

#Loading bcftools module separately as there were conflicts when loading them all together
module load BCFtools/1.10.2-foss-2019b

bcftools norm -m - -Oz ${Out}/${Job}_results/${Job}_vcf/${Job}_filtered_comb_cohort.vcf.gz > ${Out}/${Job}_report/${Job}_filtered_cohort.vcf.gz


##This is the filtering that was done before April 2021 (selecting the variants with at least 30 reads)
#bcftools filter -i'FMT/DP>30'  | \
#bcftools view -f PASS > ${Out}/${Job}_results/${Job}_vcf/${Job}_filtered_cohort.vcf.gz


#Plotting statistics for vcf
module load texlive/20200406-GCCcore-10.2.0
bcftools stats -s - ${Out}/${Job}_report/${Job}_filtered_cohort.vcf.gz > ${Out}/${Job}_report/${Job}_vcf_stats.vchk

#Creating pdf report generates errors 
#plot-vcfstats -p ${Out}/${Job}_report/${Job}_vcf_stats -s -T $Job ${Out}/${Job}_report/${Job}_vcf_stats.vchk



##Creating the separate pivoted vcf files (pivoting: unique varaint calls for the cohort)

#Separating the combined vcf to separate samples, to get unique variant calls:
if test "$(cat $gvcf_list | wc -l)" -gt 1 ; then
    for sample in `bcftools query -l ${Out}/${Job}_report/${Job}_filtered_cohort.vcf.gz`; do
    bcftools view -c1 -Oz -s $sample -o ${Out}/${Job}_results/${Job}_vcf/vcf_pivot/$sample.vcf.gz ${Out}/${Job}_report/${Job}_filtered_cohort.vcf.gz
    bcftools index ${Out}/${Job}_results/${Job}_vcf/vcf_pivot/$sample.vcf.gz
    echo $sample.vcf.gz >> ${Out}/${Job}_results/${Job}_vcf/vcf_pivot/samplevcf_list.txt
    done



#Getting unique variant calls
    cd ${Out}/${Job}_results/${Job}_vcf/vcf_pivot

    no_of_samples=$(wc -l < samplevcf_list.txt)
    samplevcf_list=samplevcf_list.txt

    for i in $(seq 1 $no_of_samples)
    do
        sample=`sed -n "$i"p $samplevcf_list |  awk '{print $1}'` 
        samplevcf_list1=$(grep -v $sample $samplevcf_list)
        sample_order=$(echo $sample $samplevcf_list1)
        sample_name=$(echo $sample | sed 's/.vcf.gz//') 
        bcftools isec -C $sample_order -w 1 >  ${sample_name}_unique.vcf
    done

fi


