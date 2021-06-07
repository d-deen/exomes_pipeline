#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=250G
#SBATCH -c 4


Out=$1
Job=$2

module load Java/11.0.2
module load BCFtools/1.10.2-foss-2019b
module load SAMtools/1.12-GCC-10.2.0 


mutserve_path=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD
haplogrep_path=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD
haplocheck_path=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD


mito_ref="/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/mitoSequence_DD/rCRS.fasta"
fasta_ref="/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/genomes_DD/GCA_000001405.15_GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
mito_size="chrM:1-16569"
chrM="chrM"

cd ${Out}/${Job}_results/${Job}_bam

samtools depth -a -r ${mito_size} -H -q 30  -d 0 -G UNMAP -G SECONDARY -G QCFAIL -G DUP \
--reference  ${fasta_ref}  \
-o ${Out}/${Job}_report/mito/${Job}_mito_coverage.csv \
$(ls *bwa_dedup.bam)


#Calculating coverage

ls -1 *bwa_dedup.bam > bam_list.txt
bam_list=bam_list.txt

for f in $(ls *bwa_dedup.bam)
do 
    samtools coverage --ff UNMAP,SECONDARY,QCFAIL,DUP -r ${chrM} -H $f >> ${Job}_coverage_stats.txt
done

paste -d '\t' bam_list.txt ${Job}_coverage_stats.txt > ${Job}_cov_stats.txt
echo -e "Sample_Name\trName\tstartpos\tendpos\tnumreads\tcoveredbases\tcoverage\tmeandepth\tmeanbasequality\tmeanmapquality" | cat - ${Job}_cov_stats.txt \
 > ${Out}/${Job}_report/mito/${Job}_mito_coverage_stats.txt


${mutserve_path}/mutserve call --reference ${mito_ref} \
--level 0.01 --baseQ=30 --write-raw \
--output ${Out}/${Job}_results/${Job}_vcf/${Job}_mutserve.vcf \
$(ls ${Out}/${Job}_results/${Job}_bam/*bwa_dedup.bam) 



    #Indel calling
bcftools mpileup -C 50 -d 20000 -f ${fasta_ref} --threads 4 -L 20000 -m 10 -O z \
-Q 30 -r ${chrM} --per-sample-mF -F 0.1  -b ${bam_list} -o ${Out}/${Job}_results/${Job}_vcf/${Job}.mpileup

bcftools index  ${Out}/${Job}_results/${Job}_vcf/${Job}.mpileup

bcftools call -O v -r ${chrM} -m -v -V snps --threads 4 -o ${Out}/${Job}_results/${Job}_vcf/${Job}_indels.vcf  ${Out}/${Job}_results/${Job}_vcf/${Job}.mpileup

#Removing the variants in the blacklisted positions
module load BEDTools/2.27.1-GCCcore-6.4.0

blacklist=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/mitoSequence_DD/blacklisted_regions_chrM.bed
bedtools intersect -v -header -a ${Out}/${Job}_results/${Job}_vcf/${Job}_indels.vcf \
-b ${blacklist} > ${Out}/${Job}_report/mito/${Job}_indels_mito_filtered.vcf
bedtools intersect -v -header -a ${Out}/${Job}_results/${Job}_vcf/${Job}_mutserve.vcf \
-b ${blacklist} > ${Out}/${Job}_report/mito/${Job}_mutserve_filtered.vcf


#Annotation for haplogroup

${haplocheck_path}/haplocheck --out ${Out}/${Job}_report/${Job}_contamination.txt \
${Out}/${Job}_results/${Job}_vcf/${Job}_mutserve.vcf
    
${haplogrep_path}/haplogrep classify --format vcf --phylotree 17 --extend-report \
--in ${Out}/${Job}_results/${Job}_vcf/${Job}_mutserve.vcf \
--out ${Out}/${Job}_report/mito/${Job}_haplogroup_classifcation.txt 

