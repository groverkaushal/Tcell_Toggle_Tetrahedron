Overview
Below is an overall description of our tool.

System Requirements:
python=3.11.4


Installation:
# python 3.12.1 -m venv venv_name
# source venv/bin/activate
# pip install -r requirements.txt



First run gseconverter.py
This code takes the GSE id as input and downloads the soft file, gpl file. 
Here you can either give the GSE id to the 
# input_gseid = "GSE210222"
If run successfully, It should create a tab separated .csv file of all the raw count values.
Many times for RNA-Seq data the .soft file does not contain the matrix values for the experiment, then you need to download the matrix file separately.



Now run tpm_converter.py
This code converts the counts matrix file to tpm normalised matrix.


Now we should have a tab separated tpm and log2 normalised .csv matrix file with first column as gene name and rest columns as samples.



Now we will run GSEA analysis on the normalised matrix file by running ssgsea.py
we will input the matrix.csv file, and the genesets.gmt file for the gsea function.
This code will also plot a scatter plot of all the samples with PCA values on the xaxis and yaxis.


At last we will run main.py 
This file generates all the plots showed in the paper. 









