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


vcf_output=glob.glob(os.path.join(args.vcf_path,'*_annotated_indels_mito.vcf'))

csv_name_annot=args.job+'_annotated_indels_mito.csv'

csv_output_annot=os.path.join(args.output_path,csv_name_annot)

with open(vcf_output[0],'r') as f:
    lines = f.readlines()

#Preseving comments from vcf file

comments = []
header=[]
for i in range(len(lines)):
    if lines[i].startswith('##'):
        comments.append(lines[i])
    elif lines[i].startswith('#'):
        header.append(lines[i])

#Opening .vcf file as a dataframe and checking  whether there are any indels to annotate:

try: 
    in_df=pd.read_csv(vcf_output[0],  delimiter='\t', quotechar='"', quoting=2, comment='#', header=None)
except:
    print('No indels in vcf file')
    exit()


#Dividing the dataframe into the VEP annotation part and the sample part, with vcf formnat data
in_df.columns=header[0].split()
df_A= in_df[['#CHROM','POS','REF','ALT','INFO','FORMAT']].copy()
df_B=in_df[in_df.columns[9:]]

#Parsing the VEP annotation

df_A[['INFO','VEP']]=in_df['INFO'].str.split('CSQ=', expand=True)
df_A[['Gene','Feature','SYMBOL','Existing_variation','VARIANT_CLASS','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','HGVSc','HGVSp','BIOTYPE','IMPACT','CLIN_SIG','PolyPhen','SIFT','gnomAD_AF','CADD_PHRED','CADD_RAW','MutationTaster_pred']]=df_A['VEP'].str.split('|',expand=True)


#Replacing GnomAD empty cells with '.'

df_A['gnomAD_AF']=df_A['gnomAD_AF'].replace(r'\s+',np.nan,regex=True).replace('','.')


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
    else:
        return 'unknown'

flexcols = df_B.columns.tolist()
new_cols = []

for col in flexcols:
    new_cols.append(in_df[col].apply(get_outcome).rename(col+'_ZYG'))
    new_cols.append(in_df[col])

#Write the output file into output directory  
combined=pd.concat([df_A]+new_cols,axis=1)


with open(csv_output_annot, 'w') as f:
    combined.to_csv(f, index=False)