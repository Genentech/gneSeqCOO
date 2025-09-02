# gneSeqCOO

This package contains code for applying a Cell of Origin classifier to RNASeq data. The process consists of three steps: 
	1) normalizing
	2) calculating the LPS score
	3) splitting samples into GCB/ABC based on fixed cut-points.

A parallel process exists for performing Dark Zone signature (DZsig) classification. 

# Install

```
remotes::install_github("Genentech/gneSeqCOO")
```

# Basic usage

The input for the gneSeqCOO method is a raw count matrix, formatted as a DESeqDataSet. For a standard count matrix, this can be generated via the code below:

```{r,eval=FALSE}
## counts is a matrix, with rows of features and columns of samples.
cds = DESeqDataSetFromMatrix(counts,
  colData=data.frame(ID=colnames(counts)),
  rowData=data.frame(ID=rownames(counts)),
  design=~1)
```

Note that in order for the classifier to work, it is necessary for the row names of our CDS to be in the correct format. The classifier currently accepts Refseq “Gene IDs” (formatted as “GeneID:9294”) or Ensembl IDs (ENSG....). Genes in the dataset should not be filtered before running the algorithm, as the complete genome is used to normalize individual samples. 

The COO classifier can be applied to the cds object using the single command, coo_rnaseq(). An example is below:

```{r,eval=FALSE}

coo_pred = coo_rnaseq(cds)

```

The returned object consists of three columns:
  - Sample IDs, drawn from the column names of the cds object.
  - LPS - the Linear Predictor Score, used to split samples into GCB and ABC, and
  - COO - the Cell of Origin classification; either GCB, ABC, or Unclassified


The DZsig classifier can be applied to the same `cds` object: 

```{r,eval=FALSE}

dzsig_pred = dzsig_rnaseq(cds)

```

Thresholds for the DZsig classification are as published for the NanoString DLBCL90 assay (Ennishi et al., *J Clin Oncol* 2019; Data Supplement). 

It is also possible obtain both DZsig and COO calls and generate a combined "refined COO" classification: 

```{r,eval=FALSE}

refined_pred = refined_coo_rnaseq(cds)

```

The returned object of this last function has six columns: 
  - Sample IDs, drawn from the column names of the cds object.
  - LPS_COO - the linear predictor score for the COO classifiction
  - COO - the Cell of Origin classification; either GCB, ABC, or Unclassified
  - LPS_DZsig - the linear predictor score for the DZsig classification
  - DZsig - the Dark Zone signature classification; either DZsig-POS, DZsig-IND, or DZsig-NEG
  - refined_COO - the refined Cell of Origin classification, which is assigned heirarchically: 
    - If ABC, then ABC
    - If DZsig-POS, DZsig-POS
    - Else, COO (GCB or Unclassified)

For more details, please refer to the docmumentation in the package.
