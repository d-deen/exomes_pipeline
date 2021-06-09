import glob
import pandas as pd
import numpy as np
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("vcf_path", type=Path)
parser.add_argument("output_path", type=Path)
parser.add_argument("job")
args = parser.parse_args()

mitocarta3='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/custom_ann_DD/Human.MitoCarta3.0.csv'
omim2gene='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/custom_ann_DD/mim2gene.txt'
OMIMdisease='/nobackup/proj/rtmngs/Mt_Exome_pipeline/DD/custom_ann_DD/Gene_OMIM_disease.txt'

vcf_output=glob.glob(os.path.join(args.vcf_path,'*_filtered_cohort_annotated.vcf'))
csv_name_comments=args.job+'_descr_annotated.csv'
csv_name_annot=args.job+'_annotated.csv'
csv_output_comments=os.path.join(args.output_path,csv_name_comments)
csv_output_annot=os.path.join(args.output_path,csv_name_annot)

with open(vcf_output[0],'r') as f:
    lines = f.readlines()

#Preseving comments from the vcf file

comments = []
header=[]
for i in range(len(lines)):
    if lines[i].startswith('##'):
        comments.append(lines[i])
    elif lines[i].startswith('#'):
        header.append(lines[i])

#Opening .vcf file as a dataframe
in_df=pd.read_csv(vcf_output[0],  delimiter='\t', quotechar='"', quoting=2, comment='#', header=None)

#Dividing the dataframe into the VEP annotation part and the sample part, with vcf formnat data
in_df.columns=header[0].split()
df_A= in_df[['#CHROM','POS','REF','ALT','INFO','FORMAT']].copy()
df_B=in_df[in_df.columns[9:]]

#Parsing the VEP annotation

df_A[['INFO','VEP']]=in_df['INFO'].str.split('CSQ=', expand=True)
df_A[['Gene','Feature','SYMBOL','Existing_variation','VARIANT_CLASS','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','HGVSc','HGVSp','BIOTYPE','IMPACT','CLIN_SIG','PolyPhen','SIFT','CADD_PHRED','CADD_RAW','MutationTaster_pred','REVEL','gnomAD_AF','MAX_AF','ExACpLI','LoFtool','DisGeNET_PMID','DisGeNET_SCORE','DisGeNET_disease']]=df_A['VEP'].str.split('|',expand=True)


#Making gnomAD frequencies numeric to filter on them
df_A['gnomAD_AF'] = pd.to_numeric(df_A['gnomAD_AF'])

#Adding OMIM ID annotation

omim=pd.read_csv(omim2gene, delimiter='\t', quotechar='"', quoting=2, comment='#')
omim['MIM Number']=omim['MIM Number'].astype('str')
omim_agg=omim.groupby(['Gene'])['MIM Number'].apply('&'.join).reset_index()
df_A=df_A.merge(omim_agg, on='Gene',how='left')

#Adding OMIM disease annotation
OMIMdisease=pd.read_csv(OMIMdisease)
df_A=df_A.merge(OMIMdisease, on='SYMBOL', how='left')

#Adding MitoCarta annotation
mitocarta=pd.read_csv(mitocarta3)
df_A["In_MitoCarta3"]=df_A.Gene.isin(mitocarta.EnsemblGeneID) 


#Assigning Hom and Het to the varaints in the samples
def get_outcome(s):
    if '0/1' in s:
        return 'HET'
    elif '1/0' in s:
        return 'HET'
    elif '1/2' in s: 
        return 'HET'  
    elif '2/2' in s: 
        return 'HET' 
    elif '1/1'  in s:
        return 'HOM'
    elif '1|1'  in s:
        return 'HOM'
    elif '0/0' in s:
        return 'REF'
    elif '0|0' in s:
        return 'REF'
    else:
        return 'unknown'

flexcols = df_B.columns.tolist()
new_cols = []

for col in flexcols:
    new_cols.append(in_df[col].apply(get_outcome).rename(col+'_ZYG'))
    new_cols.append(in_df[col])
    

##Final pivoting of the files:
#remove the variants with frequencies higher than 0.05, replace the empty cells with '.'
combined=pd.concat([df_A]+new_cols,axis=1)
combined=combined.loc[ ~(combined['gnomAD_AF'] >= 0.05)]
combined['gnomAD_AF']=combined['gnomAD_AF'].fillna('.')


combined=combined.loc[(combined['In_MitoCarta3'] == True )]
combined=combined.loc[(combined['Consequence'].notna())]

#Remove the consequences we are not interested in
combined=combined.loc[~(combined['Consequence']=='intron_variant')]
combined=combined.loc[~(combined['Consequence']=='upstream_gene_variant')]
combined=combined.loc[~(combined['Consequence']=='downstream_gene_variant')]
combined=combined.loc[~(combined['Consequence']=='intergenic_variant')]
combined=combined.loc[~(combined['Consequence']=='5_prime_UTR_variant')]
combined=combined.loc[~(combined['Consequence']=='3_prime_UTR_variant')]
combined=combined.loc[~(combined['Consequence']=='synonymous_variant')]
combined=combined.loc[~(combined['Consequence']=='non_coding_transcript_variant')]
combined=combined.loc[~(combined['Consequence']=='intron_variant&non_coding_transcript_variant')]
combined=combined.loc[~(combined['Consequence']=='non_coding_transcript_exon_variant')]


#Write the output file into output directory  
with open(csv_output_comments, 'w') as f:
    f.writelines(comments)

with open(csv_output_annot, 'w') as f:
    combined.to_csv(f, index=False)


