###### this code converts counts data matrix file to tpm normalized


import pandas as pd
import numpy as np


input_gseid = "GSE135390"

data = pd.read_csv("./geodata/GSE135390_raw_counts.csv")
data.iloc[:,0] = data.iloc[:,0].str.upper()


#### annotation file is read. The below codes are used to add length of the genes to the matrix file.
lengths = pd.read_csv("./geodata/hg38_gene_length_kb.txt",sep="\t",header=None, names = ['geneid','symbol','length'])
lengths.iloc[:,1] = lengths.iloc[:,1].str.upper()
data = data.merge(lengths, left_on='Gene_ID', right_on='symbol')
# # length_saved = merged_table.iloc[:,merged_table.shape[1]-1]
data = data.drop(columns = ['geneid','symbol'])
data.iloc[:, 1:] = data.iloc[:, 1:].astype(float)


#### doing tpm normalization on the matrix file
scaling_factor = data.iloc[:, -1] 
data.iloc[:, 1:-1] = data.iloc[:, 1:-1].div(scaling_factor, axis=0)
total_counts = data.iloc[:, 1:-1].sum()
total_counts_divided = total_counts / 10**6
data.iloc[:, 1:-1] = data.iloc[:, 1:-1].div(total_counts_divided, axis=1)


data = data.drop(columns = ['length'])
print(data)
data[data.columns[1:]] = data[data.columns[1:]].add(1)
data.iloc[:, 1:] = np.log2(data.iloc[:, 1:])

# data['Gene_ID'] = data['Gene_ID'].str.upper()

data.to_csv("./geodata/GSE135390_tpm_log2.csv",sep="\t",index=False)