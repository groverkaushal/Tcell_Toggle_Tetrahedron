########## this code takes the GSE id as input and downloads the soft file, gpl file. 
########## It then generates the matrix file for that GSE id including replacing the gene id with gene names.

import GEOparse
import pandas as pd
import numpy as np
import sys
import os
import subprocess


## This .py file was not used as the matrix file was downloaded from the supplementary material of GSE135390





input_gseid = "GSE210222"

##### make the matrix file with the gene id's
##### sometimes the soft files do not contain the matrix file so this code will not work.

gse = GEOparse.get_GEO(geo=input_gseid, destdir="./geodata")
# gse = GEOparse.get_GEO(filepath="./geodata/" + input_gseid + "_family.soft.gz")
gse_table = gse.pivot_samples('VALUE')



###### replacing the column names of GSM id's with the original name of the samples
sample_names = [gse.gsms[gsm_id].metadata['title'][0] for gsm_id in gse.gsms]
gse_table.columns = sample_names 


##### download the platform from GEO
gpl_id = gse.metadata['platform_id'][0]
gpl = GEOparse.get_GEO(geo=gpl_id, destdir='./geodata')
annot_table = gpl.table
annot_table = annot_table.dropna(subset=['GENE SYMBOL'])
# print(annot_table)



annot_table = annot_table.rename(columns={'GENE SYMBOL': 'gene_symbol'})

#### in case a gff file format is given, follow this:
# annot_table['GENE_SYMBOL'] = annot_table['GENE_SYMBOL'].str.extract(r'gene:([A-Z0-9]+)', expand=False)
# annot_table = annot_table.dropna(subset=['mrna_assignment'])
# annot_table['SPOT_ID'] = annot_table['SPOT_ID'].astype('int64')
# print(annot_table[['ID', 'gene_symbol']])
# annot_table[['ID', 'mrna_assignment']].to_csv("./geodata/annotationtable.csv",sep="\t",index=False)



# merge the two tables based on the 'ID_REF' column in the gsm_table and the 'ID' column in the annot_table
merged_table = gse_table.merge(annot_table[['ID', 'gene_symbol']], left_index=True, right_on='ID')



# move last 2 columns ['ID', 'gene_symbol'] to the front.
cols = merged_table.columns.tolist()
cols = cols[-2:] + cols[:-2]
merged_table = merged_table[cols]


# replace nan and none values with 0.
merged_table = pd.DataFrame(merged_table)
merged_table.fillna(0, inplace=True)



merged_table = merged_table.drop(columns = ['ID'])


###### if you want to convert log normalized values back to count
# import numpy as np
# merged_table[merged_table.columns[1:]] = np.exp2(merged_table[merged_table.columns[1:]])



# merged_table.iloc[:,0] = merged_table.iloc[:,0].str.upper()

temp = "./geodata/" + input_gseid + "_tabseparated.csv"
merged_table.to_csv(temp, index = False, sep="\t")
