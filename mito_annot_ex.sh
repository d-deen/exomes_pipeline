#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=50G
#SBATCH -c 4

#this script annotates mito SNPs and indels

Out=$1
Job=$2
script_path=$3


module load Perl/5.30.2-GCCcore-9.3.0
module load tabix/0.2.6-foss-2017b
module load Bio-DB-HTS/3.01-GCC-8.2.0-2.31.1-Perl-5.28.1
module load DBD-mysql/4.050-foss-2019a-Perl-5.28.1

vep_path='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep'

export PERL5LIB=$PERL5LIB:/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins

fasta="/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/homo_sapiens/103_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
species="homo_sapiens"

${vep_path}/vep --cache --dir /nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep \
    --dir_cache /nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep \
    --offline \
    --fasta ${fasta} \
    --species ${species} \
    --input_file ${Out}/${Job}_report/mito/${Job}_mutserve_filtered.vcf   \
    --output_file ${Out}/${Job}_report/mito/${Job}_annotated_snp_mito.vcf  \
    --format vcf \
    --no_stats \
    --force_overwrite  \
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
    --plugin CADD,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/whole_genome_SNVs.tsv.gz,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/InDels.tsv.gz \
    --plugin dbNSFP,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/dbNSFP4.1a_grch38.gz,MutationTaster_pred \
    --fields "Gene,Feature,SYMBOL,Existing_variation,VARIANT_CLASS,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,HGVSc,HGVSp,BIOTYPE,IMPACT,CLIN_SIG,PolyPhen,SIFT,gnomAD_AF,CADD_PHRED,CADD_RAW,MutationTaster_pred" \
    --pick \
    --pick_order rank,canonical,tsl \
    --buffer_size 20000 \
    --fork 4 

    module load Python/3.7.0-foss-2018b
    python ${script_path}/format_mito_ex.py ${Out}/${Job}_report/mito ${Out}/${Job}_report/mito ${Job} 


${vep_path}/vep --cache --dir /nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep \
    --dir_cache /nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep \
    --offline \
    --fasta ${fasta} \
    --species ${species} \
    --input_file ${Out}/${Job}_report/mito/${Job}_indels_mito_filtered.vcf   \
    --output_file ${Out}/${Job}_report/mito/${Job}_annotated_indels_mito.vcf  \
    --format vcf \
    --no_stats \
    --force_overwrite  \
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
    --plugin CADD,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/whole_genome_SNVs.tsv.gz,/nobackup/proj/rtmngs/Mt_Exome_pipeline/genomes/Plugins/CADD/InDels.tsv.gz \
    --plugin dbNSFP,/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/programs_DD/ensembl-vep/Plugins/dbNSFP4.1a_grch38.gz,MutationTaster_pred \
    --fields "Gene,Feature,SYMBOL,Existing_variation,VARIANT_CLASS,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,HGVSc,HGVSp,BIOTYPE,IMPACT,CLIN_SIG,PolyPhen,SIFT,gnomAD_AF,CADD_PHRED,CADD_RAW,MutationTaster_pred" \
    --pick \
    --pick_order rank,canonical,tsl \
    --buffer_size 20000 \
    --fork 4   
    
python  ${script_path}/format_mito_indels_ex.py ${Out}/${Job}_report/mito ${Out}/${Job}_report/mito ${Job} 