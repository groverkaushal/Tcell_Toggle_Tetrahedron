## Overview

Below is an overall description of our tool. For details and citations about this work, please check:



## System Requirements

- Python 3.11.4

## Installation

1. Create a virtual environment:

   ```bash
   python 3.12.1 -m venv venv_name
   ```

2. Activate the virtual environment:

   ```bash
   source venv/bin/activate
   ```

3. Install required packages:

   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Step 1: Download Raw Count Values

Run `gseconverter.py` by providing the GSE id as input:

```python
input_gseid = "GSE210222"
```

This code will download soft and gpl files and create a tab-separated .csv file with raw count values.

### Step 2: Convert to TPM

Execute `tpm_converter.py` to convert the counts matrix file to a TPM-normalized matrix.

### Step 3: GSEA Analysis

Run `ssgsea.py` by providing the matrix.csv file and genesets.gmt file for GSEA analysis. The code will generate a scatter plot of samples with PCA values.

### Step 4: Generate Plots

Finally, execute `main.py` to generate plots showcased in the paper.

## Important Note

For RNA-Seq data where the .soft file lacks matrix values, download the matrix file separately.

## License

This project is licensed under the [MIT License](LICENSE).
