import pandas as pd
import gseapy as gp





input_gseid = "GSE135390"

######## ssgsea generator...the input files should be tab separated

ss = gp.ssgsea(data="./geodata/GSE135390_tpm_log2.csv",
       gene_sets="./geodata/genesets.gmt",
       outdir='./ssGSEAresults',
       sample_norm_method='rank',
       permutation_num=0,
       no_plot=True,
       processes=4, format='png', seed=9)





######## code to plot a scatter plot of all the samples with PCA values on the xaxis and yaxis
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os



rna_data = pd.read_csv("./geodata/"+ input_gseid +"_tpm_log2.csv", sep="\t",index_col=0)
rna_data = rna_data.T


rna_data = rna_data.reset_index()
rna_data.iloc[:, 0] = rna_data.iloc[:, 0].str.split('_').str[1]
rna_data = rna_data.sort_values(by=['index'])



#### the below list's length is the number of unique cells. and the values in the list is the values of how many samples in same type of cell.
nested_lengths= [3,3,3,3,3,3,3,3,3,3,3]
labels = ['Naive', 'Th1', 'Th1/17', 'Th17', 'Th2', 'Th22', 'Treg1', 'Treg1/17', 'Treg17', 'Treg2', 'Treg22']
# colors = ['green', 'blue', 'purple', 'red', 'yellow', 'orange', 'pink']



# colors = [color for color, count in zip(colors, nested_lengths) for _ in range(count)]
labels = [color for color, count in zip(labels, nested_lengths) for _ in range(count)]


gene_names = rna_data.iloc[:, 0].values
expression_data = rna_data.iloc[:, 1:].values


# Perform PCA on the expression data
pca = PCA(n_components=2)  # Set the number of components (2 for 2D visualization)
pca_result = pca.fit_transform(expression_data)


# Create a DataFrame to store PCA results and include gene names for reference
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
pca_df['Gene'] = gene_names



# Visualize the PCA results in a scatter plot
plt.figure(figsize=(15, 15))
plt.scatter(pca_df['PC1'], pca_df['PC2'], marker='o', alpha=0.5)
plt.xlabel('Principal Component 1 (PC1)')
plt.ylabel('Principal Component 2 (PC2)')
plt.title('PCA Visualization of RNA-seq Data')
plt.grid(True)

for i, name in enumerate(labels):
    plt.text(pca_df['PC1'][i], pca_df['PC2'][i], name, ha='center', va='bottom')

aspect_ratio = 1.5  # Width to height ratio
fig = plt.gcf()  # Get the current figure
fig.set_size_inches(10, 15 / aspect_ratio)  
fig.tight_layout(rect=(0, 0, 0.98, 1))

plt.savefig("./results/" + input_gseid + "PCA.png")
# plt.clf()
plt.show()










































# nested_lengths= [3,3,3,3,3,3,15]

# colors = [color for color, count in zip(colors, nested_lengths) for _ in range(count)]
# labels = [color for color, count in zip(labels, nested_lengths) for _ in range(count)]

# print(colors)



# gene_names = rna_data.iloc[:, 0].values
# expression_data = rna_data.iloc[:, 1:].values

# # Perform PCA on the expression data
# pca = PCA(n_components=2)  # Set the number of components (2 for 2D visualization)
# pca_result = pca.fit_transform(expression_data)

# # Create a DataFrame to store PCA results and include gene names for reference
# pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
# pca_df['Gene'] = gene_names



# # Visualize the PCA results in a scatter plot
# plt.figure(figsize=(15, 15))
# plt.scatter(pca_df['PC1'], pca_df['PC2'], marker='o', alpha=0.5,color = colors)
# plt.xlabel('Principal Component 1 (PC1)')
# plt.ylabel('Principal Component 2 (PC2)')
# plt.title('PCA Visualization of RNA-seq Data')
# plt.grid(True)

# import matplotlib.patches as mpatches
# patches = []
# for color, label in zip(colors, labels):
#     patch = mpatches.Patch(color=color, label=label)
#     patches.append(patch)

# plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=7)

# aspect_ratio = 1.5  # Width to height ratio
# fig = plt.gcf()  # Get the current figure
# fig.set_size_inches(10, 15 / aspect_ratio)  
# fig.tight_layout(rect=(0, 0, 0.98, 1))



# # Annotate some points (optional - to label specific genes)
# # Example: annotate the first 10 genes in the DataFrame
# # for i, gene in enumerate(pca_df['Gene'][:10]):
# #     plt.annotate(gene, (pca_df['PC1'][i], pca_df['PC2'][i]))

# # plt.show()
# plt.savefig("C:/Users/grover/Desktop/python_practice/iisc/results/" + input_gseid + "/1_1.PCA.png")
# # plt.clf()
# plt.show()









