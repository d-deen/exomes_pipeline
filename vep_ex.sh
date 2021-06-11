#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=100G
#SBATCH -c 4


Out=$1
Job=$2
script_path=$3


module load Perl/5.30.2-GCCcore-9.3.0
module load tabix/0.2.6-foss-2017b
module load Bio-DB-HTS/3.01-GCC-8.2.0-2.31.1-Perl-5.28.1
module load DBD-mysql/4.050-foss-2019a-Perl-5.28.1

vep_path='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep'
dir=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep
dir_cache=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep
fasta=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/homo_sapiens/103_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

export PERL5LIB=$PERL5LIB:/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins

${vep_path}/vep --cache --dir $dir \
--dir_cache $dir_cache \
--offline \
--fasta $fasta \
--species homo_sapiens \
--input_file ${Out}/${Job}_report/${Job}_filtered_cohort.vcf.gz   \
--output_file ${Out}/${Job}_report/${Job}_filtered_cohort_annotated.vcf  \
--format vcf \
--force_overwrite  \
--vcf \
--no_check_variants_order \
--check_existing \
--freq_pop gnomAD \
--assembly GRCh38 \
--stats_file ${Out}/${Job}_report/${Job}_vep_stat.html \
--warning_file ${Out}/${Job}_logs/${Job}_vep_warning.txt \
--hgvs \
--variant_class \
--keep_csq \
--af_gnomad \
--polyphen p \
--sift p \
--symbol \
--total_length \
--max_af \
--plugin CADD,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/whole_genome_SNVs.tsv.gz,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/InDels.tsv.gz \
--plugin dbNSFP,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/dbNSFP4.1a_grch38.gz,MutationTaster_pred \
--plugin ExACpLI,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/ExACpLI_values.txt \
--plugin LoFtool \
--plugin DisGeNET,file=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/all_variant_disease_pmid_associations_final.tsv.gz,disease=1 \
--plugin REVEL,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/new_tabbed_revel_grch38.tsv.gz \
--plugin Mastermind,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/mastermind_cited_variants_reference-2021.04.02-grch38.vcf.gz,0,0,1 \
--fields "Gene,Feature,SYMBOL,Existing_variation,VARIANT_CLASS,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,HGVSc,HGVSp,BIOTYPE,IMPACT,CLIN_SIG,PolyPhen,SIFT,CADD_PHRED,CADD_RAW,MutationTaster_pred,REVEL,gnomAD_AF,MAX_AF,ExACpLI,LoFtool,DisGeNET_PMID,DisGeNET_SCORE,DisGeNET_disease,Mastermind_URL" \
--pick \
--pick_order rank,canonical,tsl \
--buffer_size 20000 \
--fork 4 


module load Python/3.7.0-foss-2018b
python  ${script_path}/final_output.py ${Out}/${Job}_report ${Out}/${Job}_report ${Job} 

cd  ${Out}/${Job}_results/${Job}_vcf/vcf_pivot

if ls  ${Out}/${Job}_results/${Job}_vcf/vcf_pivot/*unique* 2> /dev/null

then

for sample in $(ls -1 ${Out}/${Job}_results/${Job}_vcf/vcf_pivot/*unique*) 
do
    ${vep_path}/vep --cache --dir $dir \
    --dir_cache $dir_cache \
    --offline \
    --no_stats \
    --fasta $fasta \
    --species homo_sapiens \
    --input_file ${sample} \
    --output_file ${sample}_coding.vcf  \
    --format vcf \
    --force_overwrite  \
    --assembly GRCh38 \
    --vcf \
    --no_check_variants_order \
    --check_existing \
    --freq_pop gnomAD \
    --hgvs \
    --variant_class \
    --keep_csq \
    --af_gnomad \
    --polyphen p \
    --sift p \
    --symbol \
    --total_length \
    --coding_only \
    --plugin CADD,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/whole_genome_SNVs.tsv.gz,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/InDels.tsv.gz \
    --plugin dbNSFP,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/dbNSFP4.1a_grch38.gz,MutationTaster_pred \
    --plugin ExACpLI,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/ExACpLI_values.txt \
    --plugin LoFtool \
    --plugin DisGeNET,file=/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/all_variant_disease_pmid_associations_final.tsv.gz,disease=1 \
    --plugin REVEL,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/new_tabbed_revel_grch38.tsv.gz \
    --plugin Mastermind,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/mastermind_cited_variants_reference-2021.04.02-grch38.vcf.gz,0,0,1 \
    --fields "Gene,Feature,SYMBOL,Existing_variation,VARIANT_CLASS,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,HGVSc,HGVSp,BIOTYPE,IMPACT,CLIN_SIG,PolyPhen,SIFT,CADD_PHRED,CADD_RAW,MutationTaster_pred,REVEL,gnomAD_AF,MAX_AF,ExACpLI,LoFtool,DisGeNET_PMID,DisGeNET_SCORE,DisGeNET_disease,Mastermind_URL" \
    --pick \
    --pick_order rank,canonical,tsl \
    --buffer_size 20000 \
    --fork 4 
done

python  ${script_path}/top_candidates.py ${Out}/${Job}_results/${Job}_vcf/vcf_pivot ${Out}/${Job}_report/${Job}_priority ${Job} 

fi
