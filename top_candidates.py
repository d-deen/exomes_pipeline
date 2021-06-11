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


vcf_output=glob.glob(os.path.join(args.vcf_path,'*_coding.vcf'))


for f in range(len(vcf_output)):
    csv_name_pivot=Path(vcf_output[f]).stem+'_pivoted.csv'
    csv_output_pivot=os.path.join(args.output_path,csv_name_pivot)
    with open(vcf_output[f],'r') as fil:
        lines = fil.readlines()

#Preseving comments from vcf file
    comments = []
    header=[]
    for i in range(len(lines)):
        if lines[i].startswith('##'):
            comments.append(lines[i])
        elif lines[i].startswith('#'):
            header.append(lines[i])

#Opening .vcf file as a dataframe
    in_df=pd.read_csv(vcf_output[f],  delimiter='\t', quotechar='"', quoting=2, comment='#', header=None)

#Dividing the dataframe into the VEP annotation part and the sample part, with vcf formnat data
    in_df.columns=header[0].split()
    df_A= in_df[['#CHROM','POS','REF','ALT','INFO','FORMAT']].copy()
    df_B=in_df[in_df.columns[9:]]

#Parsing the VEP annotation

    df_A[['INFO','VEP']]=in_df['INFO'].str.split('CSQ=', expand=True)
    df_A[['Gene','Feature','SYMBOL','Existing_variation','VARIANT_CLASS','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','HGVSc','HGVSp','BIOTYPE','IMPACT','CLIN_SIG','PolyPhen','SIFT','CADD_PHRED','CADD_RAW','MutationTaster_pred','REVEL','gnomAD_AF','MAX_AF','ExACpLI','LoFtool','DisGeNET_PMID','DisGeNET_SCORE','DisGeNET_disease','Mastermind_URL']]=df_A['VEP'].str.split('|',expand=True)




#Making GnomAD numeric and filtering by GnomAD

    df_A['gnomAD_AF'] = pd.to_numeric(df_A['gnomAD_AF'])
    

#Adding OMIM ID annotation

    omim=pd.read_csv(omim2gene, delimiter='\t', quotechar='"', quoting=2, comment='#')
    omim['MIM Number']=omim['MIM Number'].astype('str')
    omim_agg=omim.groupby(['Gene'])['MIM Number'].apply('&'.join).reset_index()
    df_A=df_A.merge(omim_agg, on='Gene',how='left')


#Adding MitoCarta annotation
    mitocarta=pd.read_csv(mitocarta3)
    df_A["In_MitoCarta3"]=df_A.Gene.isin(mitocarta.EnsemblGeneID) 
    

#Assigning Hom and Het to the samples
    def get_outcome(s):
        if '0/1' in s:
            return 'HET'
        elif '1/0' in s:
            return 'HET'
        elif '1/1'  in s:
            return 'HOM'
        elif '0/0' in s:
            return 'REF'
        elif '1/2' in s: 
            return 'HET'  
        elif '2/2' in s: 
            return 'HET' 
        else:
            return 'unknown'

    flexcols = df_B.columns.tolist()
    new_cols = []

    for col in flexcols:
        new_cols.append(in_df[col].apply(get_outcome).rename(col+'_ZYG'))
        new_cols.append(in_df[col])
    

 #Combining the annotation and Zyg information
    combined=pd.concat([df_A]+new_cols,axis=1)
   

#Removing the variants without annotation, leaving the ones in mitocarta3; pivoting by gnomAD
    combined=combined.loc[(combined['gnomAD_AF'].isnull()) | (combined['gnomAD_AF']<= 0.01)]
    combined=combined.loc[(combined['Consequence'].notna())]
    combined=combined.loc[(combined['In_MitoCarta3'] == True )]
    combined=combined.loc[(combined['Consequence'].isin(['synonymous_variant']) == False )]
    

#Write the output file into output directory 
    with open(csv_output_pivot, 'w') as out:
        combined.to_csv(out, index=False)





